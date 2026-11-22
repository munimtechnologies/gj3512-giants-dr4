#!/usr/bin/env python3
"""Amendment 2 calibration: the RV-amplitude-constrained inclination estimator
(gaiaorb/rvincl.py) and the shape-only estimator of rvdr4.py v0.2.0, side by side on
the same 300 Epochalypse systems as the v1 calibration (Yahalomi et al. 2026,
arXiv:2611.06915; Zenodo 10.5281/zenodo.23789041, CC BY 4.0).
Rules: PREREGISTRATION-amendment2.md (fixed 2026-11-22, before this job runs and
before Gaia DR4).

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
Reuses the frozen I/O, selection and v1 evaluation of epocha_job.py unchanged; the
shape-only numbers must reproduce the v1 run exactly (reported as a check).

Usage: python epocha2_job.py --data-dir /data --out-dir /out [--n-sys 300]
"""
import argparse
import json
import math
import os
import random
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import epocha_job as E  # noqa: E402
from gaiaorb import rvdr4 as R  # noqa: E402
from gaiaorb import rvincl as V  # noqa: E402

VERSION = "0.2.0"
T0 = time.time()
LOW_INC_DEG = 25.0
SIG_S_REL_PERT, SIG_OM_PERT = 0.03, 3.0     # RV amplitude and omega errors, perturbed mode
S_REL_ORACLE = 0.02                         # assumed (not applied) amplitude error, oracle mode
C_OVERALL_MIN, C_LOW_MIN, C_MASS_MIN = 0.85, 0.80, 0.85
MIN_N, MIN_LOW = 30, 15
V1_ORACLE_COVERAGE = 0.9049773755656109     # rvdr4-epochalypse-v1, oracle, P <= 11 yr


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


# ------------------------------------------------------------------ pure part
def ce_inputs(planet, star, mode, rng2):
    """RV amplitude, omega and m sin i handed to the constrained estimator.
    The perturbed mode scales s and m sin i by the same factor (both are
    proportional to K) and shifts omega; two draws are always consumed."""
    s_true = R.a0_min_mas(planet, star)
    if mode == "perturbed":
        f = 1.0 + SIG_S_REL_PERT * rng2.gauss(0.0, 1.0)
        dom = SIG_OM_PERT * rng2.gauss(0.0, 1.0)
        return dict(s=s_true * f, s_rel=SIG_S_REL_PERT, omega=planet["omega"] + dom, omega_sig=SIG_OM_PERT,
                    msini=planet["msini"] * f, msini_sig=planet["msini_sig"] * f)
    return dict(s=s_true, s_rel=S_REL_ORACLE, omega=planet["omega"], omega_sig=0.0,
                msini=planet["msini"], msini_sig=planet["msini_sig"])


def ce_row(post, incf_true, mass_true):
    if post is None:
        return dict(ce_status="not_run")
    if post.get("status") != "ok":
        return dict(ce_status=post.get("status", "error"))
    s, m = post["inc_fold_deg"], post["mass_true_mjup"]
    ev = post["log_evidence_branch"]
    return dict(ce_status="ok", ce_inc_med=s["median"], ce_inc_p02=s["p02"], ce_inc_p98=s["p98"],
                ce_covered_inc=V.covered(s, incf_true), ce_mass_med=m["median"],
                ce_covered_mass=V.covered(m, mass_true), ce_branch=post["preferred_branch"],
                ce_dlogev=ev["omega"] - ev["minus_omega"], ce_chi2min=post["chi2_min_constrained"])


def evaluate2(points, tr, mode, rng, rng2, ce_fn=V.posterior):
    """The v1 row of epocha_job.evaluate plus the constrained estimator."""
    planet, star = E.planet_star_from_truth(tr, mode, rng)
    R.STARS[star["key"]] = star
    res = R.run(points, planets=(planet,), star=star, log=E.quiet)
    p = res["planets"]["b"]
    v = R.verdict(p, True)
    post = p["posterior"]
    incf_true = R.fold_inc(tr["inc"])
    s, m = post["inc_fold_deg"], post["mass_true_mjup"]
    row = dict(sid=tr["sid"], mode=mode, P_yr=tr["P_yr"], ecc=tr["ecc"], alpha_mas=tr["alpha_mas"],
               inc_fold_true=incf_true, mass_true=tr["mpl"], n=res["n_transits"], jitter_mas=res["jitter_mas"],
               dchi2=p["dchi2"], lam=p["expected_dchi2_min_mass"], z=p["a0_sini_z"], verdict=v,
               a0_point=p["point"]["a0_mas"], inc_fold_med=s["median"], inc_fold_p02=s["p02"],
               inc_fold_p98=s["p98"], covered_inc=bool(s["p02"] <= incf_true <= s["p98"]),
               mass_med=m["median"], covered_mass=bool(m["p02"] <= tr["mpl"] <= m["p98"]))
    ci = ce_inputs(planet, star, mode, rng2)
    ce = None
    if v.startswith("TRUE_MASS_MEASURED"):
        ti = [p["ti"][k] for k in ("A", "B", "F", "G")]
        try:
            ce = ce_fn(ti, p["ti_cov"], ci["s"], ci["s_rel"], ci["omega"], ci["omega_sig"], ci["msini"],
                       ci["msini_sig"])
        except Exception as e:  # counted as not covering, never dropped
            ce = dict(status="error: %r" % (e,))
    row.update(ce_row(ce, incf_true, tr["mpl"]))
    return row


def write_csv(path, rows):
    """CSV with the union of all row keys (rows without a constrained fit lack ce_*)."""
    import csv
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _frac(xs):
    return (sum(1 for x in xs if x) / len(xs)) if xs else None


def _median(xs):
    s = sorted(xs)
    if not s:
        return None
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def compare(rows, p_max=E.P_GATE_YR, inc_lo=0.0, inc_hi=90.0001):
    """Shape-only vs constrained estimator on measured systems with P <= p_max and
    inc_lo <= true folded inclination < inc_hi. A constrained run that did not
    finish counts as not covering."""
    sel = [r for r in rows if r["P_yr"] <= p_max and r["verdict"].startswith("TRUE_MASS_MEASURED")
           and inc_lo <= r["inc_fold_true"] < inc_hi]
    ok = [r for r in sel if r.get("ce_status") == "ok"]
    return dict(
        n=len(sel), n_ce_ok=len(ok),
        shape_inc_coverage=_frac([r["covered_inc"] for r in sel]),
        ce_inc_coverage=_frac([r.get("ce_status") == "ok" and r["ce_covered_inc"] for r in sel]),
        shape_mass_coverage=_frac([r["covered_mass"] for r in sel]),
        ce_mass_coverage=_frac([r.get("ce_status") == "ok" and r["ce_covered_mass"] for r in sel]),
        shape_median_abs_dinc=_median([abs(r["inc_fold_med"] - r["inc_fold_true"]) for r in sel]),
        ce_median_abs_dinc=_median([abs(r["ce_inc_med"] - r["inc_fold_true"]) for r in ok]),
        shape_median_signed_dinc=_median([r["inc_fold_med"] - r["inc_fold_true"] for r in sel]),
        ce_median_signed_dinc=_median([r["ce_inc_med"] - r["inc_fold_true"] for r in ok]),
        shape_n_high=sum(1 for r in sel if r["inc_fold_med"] > r["inc_fold_true"]),
        ce_n_high=sum(1 for r in ok if r["ce_inc_med"] > r["inc_fold_true"]),
        ce_branch_omega=sum(1 for r in ok if r["ce_branch"] == "omega"))


def gates_ce(oracle_rows):
    allc = compare(oracle_rows)
    low = compare(oracle_rows, inc_hi=LOW_INC_DEG)
    g = dict(
        C1_overall_inclination=dict(n=allc["n"], coverage=allc["ce_inc_coverage"],
                                    pass_=bool(allc["n"] >= MIN_N and allc["ce_inc_coverage"] >= C_OVERALL_MIN),
                                    enough=allc["n"] >= MIN_N),
        C2_low_inclination=dict(n=low["n"], coverage=low["ce_inc_coverage"],
                                pass_=bool(low["n"] >= MIN_LOW and low["ce_inc_coverage"] >= C_LOW_MIN),
                                enough=low["n"] >= MIN_LOW),
        C3_mass=dict(n=allc["n"], coverage=allc["ce_mass_coverage"],
                     pass_=bool(allc["n"] >= MIN_N and allc["ce_mass_coverage"] >= C_MASS_MIN),
                     enough=allc["n"] >= MIN_N))
    if not all(v["enough"] for v in g.values()):
        outcome = "CE_UNDERPOWERED"
    elif all(v["pass_"] for v in g.values()):
        outcome = "CE_PASS"
    else:
        outcome = "CE_FAIL:" + ",".join(k for k, v in g.items() if not v["pass_"])
    return g, outcome


# ------------------------------------------------------------------ run
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--n-sys", type=int, default=300)
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    res = dict(project="gj3512-giants-dr4 / amendment 2 calibration", version=VERSION,
               rvdr4_version=R.VERSION, rvincl_version=V.VERSION)
    try:
        truth_path = E.fetch_to(E.TRUTH_FILE, os.path.join(args.data_dir, E.TRUTH_FILE))
        ep1, info1 = E.stream_first_epochs(E.TAR_1, args.data_dir)
        log("epochs (1 companion): %s %s" % (ep1, info1))
        sids1 = E.epoch_sids(ep1)
        tcols, tmap, tmiss, _tn = E.read_cols(truth_path, E.TRUTH_COLS, filter_sids=sids1)
        raw = [{k: tcols[c][i] for k, c in tmap.items()} for i in range(len(tcols[tmap["sid"]]))]
        truths = sorted((E.truth_from_raw(r) for r in raw), key=lambda tr: tr["sid"])
        cands = [tr for tr in truths if E.selectable(tr)]
        ecols, emap, _em, _en = E.read_cols(ep1, E.EPOCH_COLS, filter_sids={tr["sid"] for tr in cands})
        by = E.group_rows(ecols, emap)
        chosen = [tr for tr in cands if len(by.get(tr["sid"], ())) >= E.N_MIN_EPOCHS][:args.n_sys]
        log("selectable %d, chosen %d" % (len(cands), len(chosen)))
        rng = random.Random(20261122)       # the v1 ephemeris perturbations, same sequence
        rng2 = random.Random(20261123)      # amplitude and omega perturbations (new)
        rows = []
        for k, tr in enumerate(chosen):
            pts = E.points_from_rows(by[tr["sid"]])
            for mode in ("oracle", "perturbed"):
                try:
                    rows.append(evaluate2(pts, tr, mode, rng, rng2))
                except Exception as e:  # recorded, never silently dropped
                    rows.append(dict(sid=tr["sid"], mode=mode, P_yr=tr["P_yr"], verdict="ERROR", error=repr(e)))
            if k % 25 == 0:
                log("%d/%d %s" % (k, len(chosen), rows[-2]))
        ok = [r for r in rows if r["verdict"] != "ERROR"]
        oracle = [r for r in ok if r["mode"] == "oracle"]
        pert = [r for r in ok if r["mode"] == "perturbed"]
        write_csv(os.path.join(args.out_dir, "systems_ce.csv"), ok)
        g, outcome = gates_ce(oracle)
        v1 = E.summarize_rows(oracle)
        bins = [(0.0, 15.0), (15.0, 25.0), (25.0, 45.0), (45.0, 70.0), (70.0, 90.0001)]
        res.update(
            n_chosen=len(chosen), n_errors=len(rows) - len(ok),
            errors_sample=[r for r in rows if r["verdict"] == "ERROR"][:5],
            reproduces_v1=dict(shape_oracle_coverage=v1["inc_coverage_95"],
                               v1_value=V1_ORACLE_COVERAGE,
                               same=bool(v1["inc_coverage_95"] is not None
                                         and abs(v1["inc_coverage_95"] - V1_ORACLE_COVERAGE) < 1e-12)),
            oracle=dict(all=compare(oracle), low=compare(oracle, inc_hi=LOW_INC_DEG),
                        high=compare(oracle, inc_lo=LOW_INC_DEG),
                        bins={"%g-%g" % (lo, min(hi, 90.0)): compare(oracle, inc_lo=lo, inc_hi=hi)
                              for lo, hi in bins}),
            perturbed=dict(all=compare(pert), low=compare(pert, inc_hi=LOW_INC_DEG),
                           high=compare(pert, inc_lo=LOW_INC_DEG)),
            oracle_all_periods=compare(oracle, p_max=1e9),
            gates=g, outcome=outcome, status="ok")
        log("GATES %s" % json.dumps(g))
        log("OUTCOME %s" % outcome)
    except Exception:
        res.update(status="error", outcome="ERROR", error=traceback.format_exc())
        log(res["error"])
    res["runtime_s"] = time.time() - T0
    with open(os.path.join(args.out_dir, "results_ce.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
