#!/usr/bin/env python3
"""Calibration of the frozen RV-fixed Gaia DR4 fitter (gaiaorb/rvdr4.py) on the
Epochalypse catalogue: simulated Gaia DR4-like epoch astrometry with known injected
companions (Yahalomi et al. 2026, arXiv:2611.06915; Zenodo 10.5281/zenodo.23789041,
CC BY 4.0). Rules: VALIDATION-epochalypse.md (fixed 2026-11-22, before the run).

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
The pure functions are standard library only; pyarrow is imported lazily in I/O.

Usage: python epocha_job.py --data-dir /data --out-dir /out [--n-sys 300] [--n-null 300]
"""
import argparse
import csv
import glob
import json
import math
import os
import random
import sys
import tarfile
import time
import traceback
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from gaiaorb import rvdr4 as R  # noqa: E402

VERSION = "0.1.0"
T0 = time.time()
T0_JD = 2457936.875                 # Epochalypse t0 (= J2017.5, the gaiaorb origin)
MJ_PER_MSUN = 1047.5655
ZENODO = "https://zenodo.org/api/records/23789041/files/%s/content"
TRUTH_FILE = "injected_solutions_1_companion_high_snr.parquet"
TAR_1 = "simulated_astrometry_1_companion.tar"
TAR_0 = "simulated_astrometry_0_companion.tar"

# selection (VALIDATION-epochalypse.md, section 3)
P_MIN_YR, P_MAX_YR = 0.5, 20.0
ALPHA_MIN, ALPHA_MAX = 0.1, 3.0
E_MAX = 0.6
N_MIN_EPOCHS = 40
# gates (section 5)
P_GATE_YR = 11.0                    # two DR4 baselines
COVERAGE_MIN = 0.85
INCONSISTENT_MAX = 0.10
NULL_FP_MAX = 0.01
MIN_N = 30
MIN_NULL = 100
# ephemeris errors of the "perturbed" mode (section 4)
SIG_P_REL, SIG_E, SIG_TP_FRAC = 0.01, 0.02, 0.02

TRUTH_COLS = dict(sid=("gaia_source_id",), plx=("parallax_mas",), mstar=("mass_st_msun",),
                  sma=("sma_1", "sma"), ecc=("ecc_1", "ecc"), inc=("inc_1", "inc"),
                  M0=("M_anom_1", "M_anom", "M0_1", "M0"), mpl=("mass_pl_1", "mass_pl"),
                  omega=("omega_1", "omega"), Omega=("Omega_1", "Omega"))
EPOCH_COLS = dict(sid=("gaia_source_id",), t=("obs_time_tcb",), psi=("scan_pos_angle",),
                  f=("parallax_factor_al",), w=("centroid_pos_al",), sig=("centroid_pos_error_al",))


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def quiet(*_a):
    pass


# ------------------------------------------------------------------ pure part
def resolve(names, wanted):
    """Map logical names to the first candidate column present in `names`."""
    out, missing = {}, []
    for k, cands in wanted.items():
        hit = next((c for c in cands if c in names), None)
        if hit is None:
            missing.append(k)
        else:
            out[k] = hit
    return out, missing


def norm_sid(x):
    s = str(x).strip()
    if s.endswith(".0"):
        s = s[:-2]
    return s


def truth_from_raw(raw):
    """Normalise one truth row (logical keys) and derive P and alpha."""
    tr = dict(sid=norm_sid(raw["sid"]))
    for k in ("plx", "mstar", "sma", "ecc", "inc", "M0", "mpl", "omega", "Omega"):
        tr[k] = float(raw[k]) if raw.get(k) is not None else float("nan")
    mp = tr["mpl"] / MJ_PER_MSUN
    tr["P_yr"] = math.sqrt(tr["sma"] ** 3 / (tr["mstar"] + mp))
    tr["alpha_mas"] = mp / (tr["mstar"] + mp) * tr["sma"] * tr["plx"]
    return tr


def selectable(tr):
    return (P_MIN_YR <= tr["P_yr"] <= P_MAX_YR and ALPHA_MIN <= tr["alpha_mas"] <= ALPHA_MAX
            and 0.0 <= tr["ecc"] <= E_MAX and all(math.isfinite(tr[k]) for k in ("inc", "M0", "plx", "mstar")))


def tc_from_tp(Tp, P, e, omega_deg):
    """Inverse of rvdr4.tp_from_tc (conjunction at true anomaly 90 deg - omega)."""
    f = math.pi / 2.0 - math.radians(omega_deg)
    E = 2.0 * math.atan2(math.sqrt(1.0 - e) * math.sin(f / 2.0), math.sqrt(1.0 + e) * math.cos(f / 2.0))
    M = E - e * math.sin(E)
    return Tp + P * M / (2.0 * math.pi)


def planet_star_from_truth(tr, mode="oracle", rng=None, key="sim"):
    """An RV-like planet dict (P, e, Tc, K, m sin i) and star dict for rvdr4.
    K is set so that a0_min = alpha sin i exactly (for the mode's P and e)."""
    P_true = tr["P_yr"] * 365.25
    Tp = T0_JD - math.radians(tr["M0"]) * P_true / (2.0 * math.pi)
    P, e = P_true, tr["ecc"]
    if mode == "perturbed":
        P = P_true * (1.0 + SIG_P_REL * rng.gauss(0.0, 1.0))
        e = min(0.95, max(0.0, e + SIG_E * rng.gauss(0.0, 1.0)))
        Tp = Tp + SIG_TP_FRAC * P_true * rng.gauss(0.0, 1.0)
    om = (tr["omega"] + 180.0) % 360.0 if math.isfinite(tr["omega"]) else 90.0
    si = math.sin(math.radians(tr["inc"]))
    asini_au = tr["alpha_mas"] / tr["plx"] * si
    K = asini_au * R.AU_M * 2.0 * math.pi / (P * R.DAY_S * math.sqrt(1.0 - e * e))
    msini = tr["mpl"] * si
    planet = dict(key="b", star=key, P=P, e=e, omega=om, Tc=tc_from_tp(Tp, P, e, om), K=K,
                  msini=msini, msini_sig=0.02 * msini)
    star = dict(key=key, name=key, dr3_plx=tr["plx"], mass_msun=tr["mstar"], mass_sig=0.02 * tr["mstar"],
                rv_kms=0.0)
    return planet, star


def points_from_rows(rows):
    """Epochalypse epoch rows (logical keys t, psi, f, w, sig) -> gaiaorb points."""
    pts = []
    for r in sorted(rows, key=lambda r: r["t"]):
        psi = float(r["psi"])
        pts.append(dict(t=(float(r["t"]) - T0_JD) / 365.25, s=math.sin(psi), c=math.cos(psi), f=float(r["f"]),
                        fr=0.0, sig=float(r["sig"]), w=float(r["w"])))
    return pts


def evaluate(points, tr, mode="oracle", rng=None):
    planet, star = planet_star_from_truth(tr, mode, rng)
    R.STARS[star["key"]] = star          # rvdr4.expected_lambda looks the star up by key
    res = R.run(points, planets=(planet,), star=star, log=quiet)
    p = res["planets"]["b"]
    v = R.verdict(p, True)
    post = p["posterior"]
    incf_true = R.fold_inc(tr["inc"])
    s, m = post["inc_fold_deg"], post["mass_true_mjup"]
    return dict(sid=tr["sid"], mode=mode, P_yr=tr["P_yr"], ecc=tr["ecc"], alpha_mas=tr["alpha_mas"],
                inc_fold_true=incf_true, mass_true=tr["mpl"], n=res["n_transits"], jitter_mas=res["jitter_mas"],
                dchi2=p["dchi2"], lam=p["expected_dchi2_min_mass"], z=p["a0_sini_z"], verdict=v,
                a0_point=p["point"]["a0_mas"], inc_fold_med=s["median"], inc_fold_p02=s["p02"],
                inc_fold_p98=s["p98"], covered_inc=bool(s["p02"] <= incf_true <= s["p98"]),
                mass_med=m["median"], covered_mass=bool(m["p02"] <= tr["mpl"] <= m["p98"]))


def evaluate_null(points, tr_fake):
    planet, star = planet_star_from_truth(tr_fake, "oracle")
    R.STARS[star["key"]] = star
    res = R.run(points, planets=(planet,), star=star, log=quiet)
    p = res["planets"]["b"]
    return dict(fake_from=tr_fake["sid"], P_yr=tr_fake["P_yr"], n=res["n_transits"], dchi2=p["dchi2"],
                verdict=R.verdict(p, True), detected=bool(p["dchi2"] >= R.DET_DCHI2))


def frac(xs):
    return (sum(1 for x in xs if x) / len(xs)) if xs else None


def summarize_rows(rows, p_max=P_GATE_YR):
    sel = [r for r in rows if r["P_yr"] <= p_max]
    meas = [r for r in sel if r["verdict"].startswith("TRUE_MASS_MEASURED")]
    det = [r for r in sel if r["dchi2"] >= R.DET_DCHI2]
    verdicts = {}
    for r in sel:
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
    dinc = sorted(abs(r["inc_fold_med"] - r["inc_fold_true"]) for r in meas)
    return dict(n=len(sel), n_measured=len(meas), n_detected=len(det), verdicts=verdicts,
                inc_coverage_95=frac([r["covered_inc"] for r in meas]),
                mass_coverage_95=frac([r["covered_mass"] for r in meas]),
                inconsistent_rate=frac([r["verdict"] == "ORBIT_INCONSISTENT_WITH_RV" for r in det]),
                median_abs_dinc_deg=dinc[len(dinc) // 2] if dinc else None)


def gates(oracle_rows, null_rows):
    s = summarize_rows(oracle_rows)
    fp = frac([r["detected"] for r in null_rows])
    g = dict(
        G1_null=dict(n=len(null_rows), false_positive_rate=fp,
                     pass_=bool(len(null_rows) >= MIN_NULL and fp is not None and fp <= NULL_FP_MAX),
                     enough=len(null_rows) >= MIN_NULL),
        G2_inclination_coverage=dict(n=s["n_measured"], coverage=s["inc_coverage_95"],
                                     pass_=bool(s["n_measured"] >= MIN_N and s["inc_coverage_95"] >= COVERAGE_MIN),
                                     enough=s["n_measured"] >= MIN_N),
        G3_inconsistent_rate=dict(n=s["n_detected"], rate=s["inconsistent_rate"],
                                  pass_=bool(s["n_detected"] >= MIN_N and s["inconsistent_rate"] <= INCONSISTENT_MAX),
                                  enough=s["n_detected"] >= MIN_N))
    if not all(v["enough"] for v in g.values()):
        outcome = "EPOCHALYPSE_UNDERPOWERED"
    elif all(v["pass_"] for v in g.values()):
        outcome = "EPOCHALYPSE_PASS"
    else:
        outcome = "EPOCHALYPSE_FAIL:" + ",".join(k for k, v in g.items() if not v["pass_"])
    return g, outcome


def group_rows(cols, keymap):
    """Column dict (actual names) -> {sid: [row dicts with logical keys]}."""
    by = {}
    n = len(cols[keymap["sid"]])
    for i in range(n):
        r = {k: cols[c][i] for k, c in keymap.items()}
        by.setdefault(norm_sid(r["sid"]), []).append(r)
    return by


# ------------------------------------------------------------------ I/O
def fetch_to(name, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    req = urllib.request.Request(ZENODO % name, headers={"User-Agent": "mv2-research-job/0.1"})
    tmp = dest + ".part"
    with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as fh:
        while True:
            b = resp.read(1 << 20)
            if not b:
                break
            fh.write(b)
    os.replace(tmp, dest)
    return dest


def stream_first_epochs(tar_name, dest_dir, max_members=800):
    """Stream the tarball from Zenodo and keep only its first epochs_rank*.parquet."""
    pop = tar_name.replace("simulated_astrometry_", "").replace(".tar", "")
    have = sorted(glob.glob(os.path.join(dest_dir, "%s-epochs_rank*.parquet" % pop)))
    if have:
        return have[0], dict(cached=True)
    req = urllib.request.Request(ZENODO % tar_name, headers={"User-Agent": "mv2-research-job/0.1"})
    seen = []
    with urllib.request.urlopen(req, timeout=600) as resp:
        tf = tarfile.open(fileobj=resp, mode="r|*")
        for k, m in enumerate(tf):
            if len(seen) < 20:
                seen.append([m.name, m.size])
            if m.isfile() and "epochs_rank" in m.name and m.name.endswith(".parquet"):
                dest = os.path.join(dest_dir, "%s-%s" % (pop, os.path.basename(m.name)))
                src = tf.extractfile(m)
                with open(dest + ".part", "wb") as fh:
                    while True:
                        b = src.read(1 << 20)
                        if not b:
                            break
                        fh.write(b)
                os.replace(dest + ".part", dest)
                return dest, dict(cached=False, member=m.name, size=m.size, first_members=seen)
            if k >= max_members:
                break
    raise RuntimeError("no epochs_rank parquet among the first members: %r" % seen)


def read_cols(path, wanted, filter_sids=None):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    names = pq.read_schema(path).names
    keymap, missing = resolve(names, wanted)
    if "sid" in missing:
        raise KeyError("no source id column in %s: %r" % (path, names))
    t = pq.read_table(path, columns=sorted(set(keymap.values())))
    if filter_sids is not None:
        sid = pc.cast(t.column(keymap["sid"]), pa.string())
        t = t.filter(pc.is_in(sid, value_set=pa.array(sorted(filter_sids), type=pa.string())))
    cols = {c: t.column(c).to_pylist() for c in set(keymap.values())}
    return cols, keymap, missing, names


def epoch_sids(path):
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    t = pq.read_table(path, columns=["gaia_source_id"])
    return set(norm_sid(x) for x in pc.unique(pc.cast(t.column(0), pa.string())).to_pylist())


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--n-sys", type=int, default=300)
    ap.add_argument("--n-null", type=int, default=300)
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    res = dict(project="gj3512-giants-dr4 / epochalypse calibration", version=VERSION, rvdr4_version=R.VERSION)
    try:
        truth_path = fetch_to(TRUTH_FILE, os.path.join(args.data_dir, TRUTH_FILE))
        ep1, info1 = stream_first_epochs(TAR_1, args.data_dir)
        log("epochs (1 companion): %s %s" % (ep1, info1))
        res["epochs_1"] = dict(path=os.path.basename(ep1), **{k: v for k, v in info1.items()})
        sids1 = epoch_sids(ep1)
        tcols, tmap, tmiss, tnames = read_cols(truth_path, TRUTH_COLS, filter_sids=sids1)
        res["truth_columns"] = dict(resolved=tmap, missing=tmiss, all=tnames)
        if any(k in tmiss for k in ("plx", "mstar", "sma", "ecc", "inc", "M0", "mpl")):
            raise KeyError("truth columns missing: %r" % tmiss)
        raw = [{k: tcols[c][i] for k, c in tmap.items()} for i in range(len(tcols[tmap["sid"]]))]
        truths = sorted((truth_from_raw(r) for r in raw), key=lambda tr: tr["sid"])
        cands = [tr for tr in truths if selectable(tr)]
        log("shard systems %d, high-SNR truths in shard %d, selectable %d" % (len(sids1), len(truths), len(cands)))
        res["counts_1"] = dict(shard_systems=len(sids1), high_snr_in_shard=len(truths), selectable=len(cands))
        ecols, emap, emiss, enames = read_cols(ep1, EPOCH_COLS, filter_sids={tr["sid"] for tr in cands})
        res["epoch_columns"] = dict(resolved=emap, missing=emiss, all=enames)
        by = group_rows(ecols, emap)
        chosen = [tr for tr in cands if len(by.get(tr["sid"], ())) >= N_MIN_EPOCHS][:args.n_sys]
        log("chosen %d systems" % len(chosen))
        rng = random.Random(20261122)
        rows = []
        for k, tr in enumerate(chosen):
            pts = points_from_rows(by[tr["sid"]])
            for mode in ("oracle", "perturbed"):
                try:
                    rows.append(evaluate(pts, tr, mode, rng))
                except Exception as e:  # recorded, never silently dropped
                    rows.append(dict(sid=tr["sid"], mode=mode, P_yr=tr["P_yr"], verdict="ERROR", error=repr(e),
                                     dchi2=float("nan")))
            if k % 25 == 0:
                log("%d/%d %s" % (k, len(chosen), rows[-2]))
        ok = [r for r in rows if r["verdict"] != "ERROR"]
        oracle = [r for r in ok if r["mode"] == "oracle"]
        pert = [r for r in ok if r["mode"] == "perturbed"]
        write_csv(os.path.join(args.out_dir, "systems_1_companion.csv"), ok)
        res["n_errors"] = len(rows) - len(ok)
        res["errors_sample"] = [r for r in rows if r["verdict"] == "ERROR"][:5]

        ep0, info0 = stream_first_epochs(TAR_0, args.data_dir)
        log("epochs (0 companion): %s %s" % (ep0, info0))
        res["epochs_0"] = dict(path=os.path.basename(ep0), **{k: v for k, v in info0.items()})
        sids0 = sorted(epoch_sids(ep0))[:3 * args.n_null]
        c0, m0, _miss0, _n0 = read_cols(ep0, EPOCH_COLS, filter_sids=set(sids0))
        by0 = group_rows(c0, m0)
        null_ids = [s for s in sids0 if len(by0.get(s, ())) >= N_MIN_EPOCHS][:args.n_null]
        nulls = []
        for k, sid in enumerate(null_ids):
            if not chosen:
                break
            try:
                r = evaluate_null(points_from_rows(by0[sid]), chosen[k % len(chosen)])
                r["sid"] = sid
                nulls.append(r)
            except Exception as e:
                res.setdefault("null_errors", []).append([sid, repr(e)])
        write_csv(os.path.join(args.out_dir, "nulls_0_companion.csv"), nulls)

        g, outcome = gates(oracle, nulls)
        res.update(selection=dict(P_yr=[P_MIN_YR, P_MAX_YR], alpha_mas=[ALPHA_MIN, ALPHA_MAX], e_max=E_MAX,
                                  n_min_epochs=N_MIN_EPOCHS, n_chosen=len(chosen)),
                   oracle=summarize_rows(oracle), perturbed=summarize_rows(pert),
                   oracle_all_periods=summarize_rows(oracle, p_max=1e9),
                   perturbed_all_periods=summarize_rows(pert, p_max=1e9),
                   gates=g, outcome=outcome, status="ok")
        log("GATES %s" % json.dumps(g))
        log("OUTCOME %s" % outcome)
    except Exception:
        res.update(status="error", outcome="ERROR", error=traceback.format_exc())
        log(res["error"])
    res["runtime_s"] = time.time() - T0
    with open(os.path.join(args.out_dir, "results_epochalypse.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
