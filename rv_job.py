#!/usr/bin/env python3
"""The pre-registered Gaia DR4 measurement of the true masses of GJ 3512 b and c.
To be run on or after the Gaia DR4 release (2026-12-02).

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
Rules: PREREGISTRATION.md. Only the I/O (source lookup, download URLs and file
format, shared with pceb-planets-dr4 through etv_job.load_points) may be adapted
after the release; the analysis may not.

Usage:
  python rv_job.py --data-dir /data --out-dir /out [--validated]
        [--source-id gj3512=N] [--epochs-file gj3512=PATH]
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

VERSION = "0.1.0"
T0 = time.time()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--validated", action="store_true",
                    help="set only if the proxima-dr4-validate verdict was VALIDATED")
    ap.add_argument("--source-id", action="append", help="gj3512=N")
    ap.add_argument("--epochs-file", action="append", help="gj3512=PATH")
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    sids, files = etv_job.kv(args.source_id), etv_job.kv(args.epochs_file)
    star = rvdr4.STARS["gj3512"]
    res = dict(project="gj3512-giants-dr4", version=VERSION, validated_flag=bool(args.validated),
               a0_min_mas={p["key"]: rvdr4.a0_min_mas(p) for p in rvdr4.PLANETS})
    try:
        sid, src, points, counts = etv_job.load_points(star, args, sids, files)
        res.update(source_id=sid, source=src, counts=counts)
        if not points:
            res["status"] = "no_epochs"
        else:
            etv_job.write_csv(os.path.join(args.out_dir, "points_gj3512.csv"), points)
            res.update(rvdr4.run(points, log=etv_job.log))
    except Exception:
        res.update(status="error", error=traceback.format_exc())
    res["verdicts"] = {p["key"]: (rvdr4.verdict(res["planets"][p["key"]], args.validated, res.get("status"))
                                  if p["key"] in res.get("planets", {}) else "INCONCLUSIVE")
                       for p in rvdr4.PLANETS}
    res["runtime_s"] = time.time() - T0
    with open(os.path.join(args.out_dir, "results_gj3512.json"), "w") as f:
        json.dump(etv_job.jsonable(res), f, indent=1)
    for k, v in res["verdicts"].items():
        etv_job.log("VERDICT GJ 3512 %s %s" % (k, v))


if __name__ == "__main__":
    main()
