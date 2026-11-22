#!/usr/bin/env python3
"""Amendment 1 (v0.2.0): the pre-registered Gaia DR4 true masses of GJ 317 b and c,
GJ 463 b and GJ 9773 b, the test of the published Gaia DR3-based inclination
predictions for these planets and for GJ 3512 b and c, and the mutual inclinations
of the two-planet systems. To be run on or after the Gaia DR4 release (2026-12-02).

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
Rules: PREREGISTRATION-amendment1.md (and PREREGISTRATION.md for GJ 3512). Only the
I/O (source lookup, download URLs and file format, shared with pceb-planets-dr4
through etv_job.load_points) may be adapted after the release; the analysis may not.
rv_job.py stays the primary GJ 3512 run; its verdicts are recomputed here by the same
code, only to attach the prediction checks.

Usage:
  python rv_job_a1.py --data-dir /data --out-dir /out [--validated]
        [--systems gj317,gj463,gj9773,gj3512] [--source-id STAR=N] [--epochs-file STAR=PATH]
"""
import argparse
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import etv_job  # noqa: E402  (shared DR4 source lookup and epoch download)
from gaiaorb import rvdr4  # noqa: E402

VERSION = "0.2.0"
T0 = time.time()
ALL_SYSTEMS = ("gj317", "gj463", "gj9773", "gj3512")


def system_planets(key):
    return rvdr4.PLANETS if key == "gj3512" else rvdr4.SYSTEMS_A1[key]


def alt_planets(planets):
    """The same planets with the alternative omega where Amendment 1 names one."""
    out, changed = [], False
    for pl in planets:
        alt = rvdr4.OMEGA_ALTERNATIVE.get((pl["star"], pl["key"]))
        if alt is not None:
            pl = dict(pl, omega=alt)
            changed = True
        out.append(pl)
    return tuple(out), changed


def analyse(points, key, validated, log=print, star=None):
    """All pre-registered numbers for one system from its DR4 normal points.
    `star` overrides the catalogue star only in tests (synthetic epochs)."""
    planets = system_planets(key)
    res = rvdr4.run(points, planets=planets, star=star, log=log)
    verdicts = {pl["key"]: rvdr4.verdict(res["planets"][pl["key"]], validated, res["status"]) for pl in planets}
    alts, changed = alt_planets(planets)
    if changed:
        res_alt = rvdr4.run(points, planets=alts, star=star, log=log)
        res["omega_alternative"] = dict(jitter_mas=res_alt["jitter_mas"], chi2_full=res_alt["chi2_full"],
                                        planets={k: dict(dchi2=v["dchi2"], a0_sini_z=v["a0_sini_z"],
                                                         point=v["point"], posterior=v["posterior"])
                                                 for k, v in res_alt["planets"].items()})
        for pl in alts:
            if (pl["star"], pl["key"]) in rvdr4.OMEGA_ALTERNATIVE:
                v_alt = rvdr4.verdict(res_alt["planets"][pl["key"]], validated, res_alt["status"])
                res["omega_alternative"]["verdict_" + pl["key"]] = v_alt
                verdicts[pl["key"]] = rvdr4.combine_omega_verdicts(verdicts[pl["key"]], v_alt)
    res["verdicts"] = verdicts
    res["prediction_checks"] = {
        pl["key"]: rvdr4.prediction_checks(res["planets"][pl["key"]], verdicts[pl["key"]],
                                           rvdr4.PREDICTIONS.get((key, pl["key"]), ()))
        for pl in planets}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--validated", action="store_true",
                    help="set only if the proxima-dr4-validate verdict was VALIDATED")
    ap.add_argument("--systems", default=",".join(ALL_SYSTEMS))
    ap.add_argument("--source-id", action="append", help="STAR=N")
    ap.add_argument("--epochs-file", action="append", help="STAR=PATH")
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    sids, files = etv_job.kv(args.source_id), etv_job.kv(args.epochs_file)
    summary = dict(project="gj3512-giants-dr4 amendment 1", version=VERSION, validated_flag=bool(args.validated),
                   systems={})
    for key in [s.strip() for s in args.systems.split(",") if s.strip()]:
        star = rvdr4.STARS[key]
        planets = system_planets(key)
        res = dict(system=key, a0_min_mas={p["key"]: rvdr4.a0_min_mas(p) for p in planets})
        try:
            sid, src, points, counts = etv_job.load_points(star, args, sids, files)
            res.update(source_id=sid, source=src, counts=counts)
            if not points:
                res["status"] = "no_epochs"
            else:
                etv_job.write_csv(os.path.join(args.out_dir, "points_%s.csv" % key), points)
                res.update(analyse(points, key, args.validated, log=etv_job.log))
        except Exception:
            res.update(status="error", error=traceback.format_exc())
        if "verdicts" not in res:
            res["verdicts"] = {p["key"]: "INCONCLUSIVE" for p in planets}
        with open(os.path.join(args.out_dir, "results_a1_%s.json" % key), "w") as f:
            json.dump(etv_job.jsonable(res), f, indent=1)
        summary["systems"][key] = dict(status=res.get("status"), verdicts=res["verdicts"],
                                       prediction_checks=res.get("prediction_checks"),
                                       mutual_inclination=res.get("mutual_inclination"))
        for k, v in res["verdicts"].items():
            etv_job.log("VERDICT %s %s %s" % (star["name"], k, v))
    summary["runtime_s"] = time.time() - T0
    with open(os.path.join(args.out_dir, "results_a1_summary.json"), "w") as f:
        json.dump(etv_job.jsonable(summary), f, indent=1)


if __name__ == "__main__":
    main()
