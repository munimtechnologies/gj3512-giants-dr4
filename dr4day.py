#!/usr/bin/env python3
"""Gaia DR4 release day: run MV2's four pre-registered near-miss tests in ONE job.

MV2 / Munim, Inc. (mv2@munimtech.com). Apache-2.0. Made with AI agents working in MV2's name.

Order (fixed 2026-11-22, before DR4):
  0. preflight  - is gaiadr4.gaia_source queryable? (if not: stop, nothing is graded)
  1. validate   - proxima-c-dr4 validate_job.py on the DR4 pre-release (--skip-rehearsal);
                  its verdict VALIDATED is the ONLY way --validated is passed to steps 2-7,
                  exactly as every pre-registration requires
  2. proxima    - proxima-c-dr4      proxima_job.py
  3. radio      - radio-planets-dr4  radio_job.py
  4. pceb       - pceb-planets-dr4   etv_job.py
  5. gj3512     - gj3512-giants-dr4  rv_job.py
  6. gj3512a1   - gj3512-giants-dr4  rv_job_a1.py (Amendment 1; writes points_<system>.csv)
  7. gj3512ce   - gj3512-giants-dr4  ce_dr4.py --ce-outcome <outcome of rvdr4-epochalypse-ce>
                  (Amendment 2; any outcome other than CE_PASS keeps the shape-only estimator primary)
A failing step is recorded and the next one still runs. Nothing here changes any analysis;
each script is the one committed in its repository before DR4.

Usage: python dr4day.py --data-dir /data --out-dir /out --ce-outcome CE_PASS|CE_FAIL|NOT_RUN
"""
import argparse
import json
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request

VERSION = "0.1.0"
T0 = time.time()
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
STEP_TIMEOUT_S = {"validate": 3600, "proxima": 5400, "radio": 5400, "pceb": 7200,
                  "gj3512": 3600, "gj3512a1": 7200, "gj3512ce": 7200}
VERDICT_KEYS = ("verdict", "verdicts", "outcome", "status", "n_passed", "n_gates")


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def plan(data_dir, out_dir, validated, ce_outcome, py=sys.executable):
    """The steps after validation, in order: (name, argv)."""
    v = ["--validated"] if validated else []

    def dd(name):
        return ["--data-dir", os.path.join(data_dir, name), "--out-dir", os.path.join(out_dir, name)]
    return [
        ("proxima", [py, "proxima_job.py"] + dd("proxima") + v),
        ("radio", [py, "radio_job.py"] + dd("radio") + v),
        ("pceb", [py, "etv_job.py"] + dd("pceb") + v),
        ("gj3512", [py, "rv_job.py"] + dd("gj3512") + v),
        ("gj3512a1", [py, "rv_job_a1.py"] + dd("gj3512a1") + v),
        ("gj3512ce", [py, "ce_dr4.py", "--points-dir", os.path.join(out_dir, "gj3512a1"),
                      "--out-dir", os.path.join(out_dir, "gj3512ce"), "--ce-outcome", ce_outcome] + v),
    ]


def validate_step(data_dir, out_dir, py=sys.executable):
    return ("validate", [py, "validate_job.py"] + ["--data-dir", os.path.join(data_dir, "validate"),
                                                   "--out-dir", os.path.join(out_dir, "validate"),
                                                   "--skip-rehearsal"])


def read_validation(out_dir):
    """True only if validate/validation.json exists and says VALIDATED."""
    p = os.path.join(out_dir, "validate", "validation.json")
    try:
        with open(p) as f:
            return json.load(f).get("verdict") == "VALIDATED"
    except (OSError, ValueError, AttributeError):
        return False


def harvest(step_out):
    """Top-level verdict-like keys of every JSON file a step wrote."""
    found = {}
    if not os.path.isdir(step_out):
        return found
    for name in sorted(os.listdir(step_out)):
        if not name.endswith(".json"):
            continue
        try:
            with open(os.path.join(step_out, name)) as f:
                d = json.load(f)
        except (OSError, ValueError):
            found[name] = "unreadable"
            continue
        if isinstance(d, dict):
            found[name] = {k: d[k] for k in VERDICT_KEYS if k in d}
    return found


def run_step(name, argv, out_dir, timeout_s, runner=subprocess.run):
    log("STEP %s: %s" % (name, " ".join(argv)))
    t = time.time()
    try:
        done = runner(argv, capture_output=True, text=True, timeout=timeout_s)
        code, so, se = done.returncode, done.stdout or "", done.stderr or ""
    except subprocess.TimeoutExpired as e:
        code, so, se = "timeout", str(e.stdout or ""), str(e.stderr or "")
    except OSError as e:
        code, so, se = "oserror", "", str(e)
    os.makedirs(os.path.join(out_dir, "logs"), exist_ok=True)
    with open(os.path.join(out_dir, "logs", name + ".txt"), "w") as f:
        f.write(so + "\n----- stderr -----\n" + se)
    rec = dict(step=name, exit=code, seconds=round(time.time() - t, 1), tail=(so + se)[-1500:],
               results=harvest(os.path.join(out_dir, name)))
    log("STEP %s exit %s" % (name, code))
    return rec


def preflight(fetch=None):
    """Is DR4 live? Returns (ok, detail)."""
    q = "SELECT TOP 1 source_id FROM gaiadr4.gaia_source"
    url = TAP + "?" + urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q))
    try:
        if fetch is None:
            req = urllib.request.Request(url, headers={"User-Agent": "gaiaorb-dr4/0.2 (open research code)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                text = r.read().decode("utf-8", "replace")
        else:
            text = fetch(url)
    except Exception as e:  # archive down or table missing
        return False, "TAP error: %s" % e
    lines = [x for x in text.strip().splitlines() if x.strip()]
    ok = len(lines) >= 2 and lines[0].strip().strip('"') == "source_id" and lines[1].strip().strip('"').isdigit()
    return ok, text[:300]


def main(argv=None, runner=subprocess.run, fetch=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--ce-outcome", required=True,
                    help="outcome in results_ce.json of rvdr4-epochalypse-ce (NOT_RUN if it has not run)")
    ap.add_argument("--skip-preflight", action="store_true")
    ap.add_argument("--only", default="", help="comma list of steps to run after validation (default all)")
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    summary = dict(project="near-misses Gaia DR4 day", version=VERSION, ce_outcome=args.ce_outcome, steps=[])

    def save():
        summary["runtime_s"] = round(time.time() - T0, 1)
        with open(os.path.join(args.out_dir, "dr4day_summary.json"), "w") as f:
            json.dump(summary, f, indent=1)

    if not args.skip_preflight:
        ok, detail = preflight(fetch)
        summary["preflight"] = dict(ok=ok, detail=detail)
        log("PREFLIGHT %s %s" % (ok, detail[:120]))
        if not ok:
            summary["outcome"] = "DR4_NOT_AVAILABLE"
            save()
            log("OUTCOME DR4_NOT_AVAILABLE")
            return summary
    name, cmd = validate_step(args.data_dir, args.out_dir)
    summary["steps"].append(run_step(name, cmd, args.out_dir, STEP_TIMEOUT_S[name], runner))
    validated = read_validation(args.out_dir)
    summary["validated"] = validated
    log("VALIDATED %s" % validated)
    save()
    only = {s.strip() for s in args.only.split(",") if s.strip()}
    for name, cmd in plan(args.data_dir, args.out_dir, validated, args.ce_outcome):
        if only and name not in only:
            continue
        summary["steps"].append(run_step(name, cmd, args.out_dir, STEP_TIMEOUT_S[name], runner))
        save()
    summary["outcome"] = "RAN"
    save()
    log("OUTCOME RAN")
    return summary


if __name__ == "__main__":
    main()
