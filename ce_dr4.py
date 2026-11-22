#!/usr/bin/env python3
"""Amendment 2: the RV-amplitude-constrained inclinations and true masses of
GJ 3512 b/c, GJ 317 b/c, GJ 463 b and GJ 9773 b from their Gaia DR4 epochs, and the
published DR3-based inclination predictions checked against them.
Run AFTER rv_job_a1.py, on the normal points it wrote (points_<system>.csv).

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
Rules: PREREGISTRATION-amendment2.md. The verdicts (detected / consistent /
true mass measured) are exactly those of rv_job_a1.py (recomputed here by the same
code); this script only adds the constrained posteriors and their checks.

Usage:
  python ce_dr4.py --points-dir /out --out-dir /out --ce-outcome <outcome of epocha2_job.py>
        [--validated] [--systems gj317,gj463,gj9773,gj3512]
"""
import argparse
import csv
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import rv_job_a1 as A1  # noqa: E402
from gaiaorb import rvdr4 as R  # noqa: E402
from gaiaorb import rvincl as V  # noqa: E402

VERSION = "0.1.0"
T0 = time.time()


def quiet(*_a):
    pass


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def read_points(path):
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            d = {}
            for k, v in r.items():
                try:
                    d[k] = float(v)
                except (TypeError, ValueError):
                    d[k] = v
            out.append(d)
    return out


def ce_inputs_real(pl, star):
    return dict(s=R.a0_min_mas(pl, star), s_rel=V.s_rel_sigma(pl), omega=pl["omega"],
                omega_sig=pl.get("omega_sig", 0.0), msini=pl["msini"], msini_sig=pl["msini_sig"])


def combine_checks(c1, c2):
    """GJ 9773 b: a check counts only if the printed and the alternative omega agree."""
    if c1 is None or c2 is None:
        return None
    out = []
    for r1, r2 in zip(c1, c2):
        chk = r1["check"] if r1["check"] == r2["check"] else "INCONCLUSIVE_OMEGA_AMBIGUOUS"
        out.append(dict(ref=r1["ref"], inc_fold_pred=r1["inc_fold_pred"], mass_pred_mjup=r1["mass_pred_mjup"],
                        check=chk, printed_omega=r1, alternative_omega=r2))
    return out


def analyse_ce(points, key, validated, ce_outcome, star=None, ce_fn=V.posterior, log_fn=quiet):
    """`star` overrides the catalogue star only in tests (synthetic epochs)."""
    planets = A1.system_planets(key)
    st = star or R.STARS[key]
    res = A1.analyse(points, key, validated, log=log_fn, star=star)
    alts, changed = A1.alt_planets(planets)
    res_alt = R.run(points, planets=alts, star=star, log=log_fn) if changed else None
    primary = "constrained" if ce_outcome == "CE_PASS" else "shape"
    out = dict(system=key, version=VERSION, rvincl_version=V.VERSION, ce_outcome=ce_outcome,
               primary_estimator=primary, validated=bool(validated), verdicts=res["verdicts"], planets={})
    for pl, pla in zip(planets, alts):
        k = pl["key"]
        v = res["verdicts"][k]
        preds = R.PREDICTIONS.get((key, k), ())
        entry = dict(verdict=v, shape_posterior=res["planets"][k]["posterior"],
                     shape_checks=res["prediction_checks"][k])
        variants = [("printed_omega", pl, res["planets"][k])]
        if pla is not pl:
            variants.append(("alternative_omega", pla, res_alt["planets"][k]))
        checks = []
        for name, plv, pv in variants:
            if pv["dchi2"] < R.DET_DCHI2:
                entry[name] = dict(status="not_detected")
                checks.append(None)
                continue
            ci = ce_inputs_real(plv, st)
            post = ce_fn([pv["ti"][q] for q in ("A", "B", "F", "G")], pv["ti_cov"], ci["s"], ci["s_rel"],
                         ci["omega"], ci["omega_sig"], ci["msini"], ci["msini_sig"])
            entry[name] = post
            checks.append(R.prediction_checks(dict(posterior=post), v, preds) if post.get("status") == "ok"
                          else None)
        ce_checks = checks[0] if len(checks) == 1 else combine_checks(checks[0], checks[1])
        entry["ce_checks"] = ce_checks
        entry["checks"] = ce_checks if primary == "constrained" else entry["shape_checks"]
        out["planets"][k] = entry
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points-dir", default="/out")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--ce-outcome", required=True, help="outcome string of epocha2_job.py (results_ce.json)")
    ap.add_argument("--validated", action="store_true",
                    help="set only if the proxima-dr4-validate verdict was VALIDATED")
    ap.add_argument("--systems", default=",".join(A1.ALL_SYSTEMS))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    summary = dict(project="gj3512-giants-dr4 amendment 2", version=VERSION, ce_outcome=args.ce_outcome,
                   validated_flag=bool(args.validated), systems={})
    for key in [s.strip() for s in args.systems.split(",") if s.strip()]:
        path = os.path.join(args.points_dir, "points_%s.csv" % key)
        try:
            if not os.path.exists(path):
                res = dict(system=key, status="no_points")
            else:
                res = analyse_ce(read_points(path), key, args.validated, args.ce_outcome, log_fn=log)
                res["status"] = "ok"
        except Exception:
            res = dict(system=key, status="error", error=traceback.format_exc())
        with open(os.path.join(args.out_dir, "results_ce_dr4_%s.json" % key), "w") as f:
            json.dump(A1.etv_job.jsonable(res), f, indent=1)
        summary["systems"][key] = dict(status=res.get("status"), verdicts=res.get("verdicts"),
                                       checks={k: e.get("checks") for k, e in res.get("planets", {}).items()})
        log("SYSTEM %s %s" % (key, res.get("status")))
    summary["runtime_s"] = time.time() - T0
    with open(os.path.join(args.out_dir, "results_ce_dr4_summary.json"), "w") as f:
        json.dump(A1.etv_job.jsonable(summary), f, indent=1)


if __name__ == "__main__":
    main()
