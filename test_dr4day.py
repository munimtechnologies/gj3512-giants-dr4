"""Tests for dr4day.py (Gaia DR4 release-day orchestrator). Standard library only; no network,
no real subprocesses. MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name."""
import json
import os
import subprocess
import tempfile
import unittest

import dr4day


class Done:
    def __init__(self, code=0):
        self.returncode, self.stdout, self.stderr = code, "ok", ""


class FakeRunner:
    """Pretends to be each job: writes a small JSON into its --out-dir."""

    def __init__(self, verdict="VALIDATED", fail=()):
        self.calls, self.verdict, self.fail = [], verdict, set(fail)

    def __call__(self, argv, capture_output=True, text=True, timeout=None):
        self.calls.append(list(argv))
        script = argv[1]
        if script in self.fail:
            raise subprocess.TimeoutExpired(argv, timeout)
        out = argv[argv.index("--out-dir") + 1]
        os.makedirs(out, exist_ok=True)
        name = "validation.json" if script == "validate_job.py" else "results.json"
        body = dict(verdict=self.verdict if script == "validate_job.py" else "FAKE_" + script, extra=1)
        with open(os.path.join(out, name), "w") as f:
            json.dump(body, f)
        return Done()


def live(_url):
    return "source_id\n4295806720\n"


class TestPreflight(unittest.TestCase):
    def test_live(self):
        self.assertEqual(dr4day.preflight(live)[0], True)

    def test_error_page(self):
        ok, _d = dr4day.preflight(lambda u: '<?xml version="1.0"?><VOTABLE><INFO value="ERROR"/></VOTABLE>')
        self.assertFalse(ok)

    def test_exception(self):
        def boom(_u):
            raise OSError("503")
        self.assertFalse(dr4day.preflight(boom)[0])


class TestMain(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = os.path.join(self.tmp.name, "out")
        self.data = os.path.join(self.tmp.name, "data")

    def tearDown(self):
        self.tmp.cleanup()

    def args(self, ce="CE_PASS", extra=()):
        return ["--data-dir", self.data, "--out-dir", self.out, "--ce-outcome", ce] + list(extra)

    def test_not_available_stops(self):
        r = FakeRunner()
        s = dr4day.main(self.args(), runner=r, fetch=lambda u: "Service Unavailable")
        self.assertEqual(s["outcome"], "DR4_NOT_AVAILABLE")
        self.assertEqual(r.calls, [])
        with open(os.path.join(self.out, "dr4day_summary.json")) as f:
            self.assertEqual(json.load(f)["outcome"], "DR4_NOT_AVAILABLE")

    def test_validated_order_and_flags(self):
        r = FakeRunner("VALIDATED")
        s = dr4day.main(self.args("CE_PASS"), runner=r, fetch=live)
        scripts = [c[1] for c in r.calls]
        self.assertEqual(scripts, ["validate_job.py", "proxima_job.py", "radio_job.py", "etv_job.py",
                                   "rv_job.py", "rv_job_a1.py", "ce_dr4.py"])
        self.assertIn("--skip-rehearsal", r.calls[0])
        self.assertNotIn("--validated", r.calls[0])
        for c in r.calls[1:]:
            self.assertIn("--validated", c)
        ce = r.calls[-1]
        self.assertEqual(ce[ce.index("--ce-outcome") + 1], "CE_PASS")
        self.assertEqual(ce[ce.index("--points-dir") + 1], os.path.join(self.out, "gj3512a1"))
        self.assertTrue(s["validated"])
        self.assertEqual(s["outcome"], "RAN")
        self.assertEqual(len(s["steps"]), 7)
        self.assertEqual(s["steps"][1]["results"], {"results.json": {"verdict": "FAKE_proxima_job.py"}})

    def test_failed_validation_never_passes_flag(self):
        r = FakeRunner("FAILED")
        s = dr4day.main(self.args("NOT_RUN"), runner=r, fetch=live)
        self.assertFalse(s["validated"])
        for c in r.calls:
            self.assertNotIn("--validated", c)

    def test_timeout_recorded_and_next_steps_run(self):
        r = FakeRunner("VALIDATED", fail={"radio_job.py"})
        s = dr4day.main(self.args(), runner=r, fetch=live)
        ex = {st["step"]: st["exit"] for st in s["steps"]}
        self.assertEqual(ex["radio"], "timeout")
        self.assertEqual(ex["pceb"], 0)
        self.assertEqual(ex["gj3512ce"], 0)
        self.assertTrue(os.path.exists(os.path.join(self.out, "logs", "radio.txt")))

    def test_only_and_skip_preflight(self):
        r = FakeRunner("VALIDATED")
        s = dr4day.main(self.args(extra=["--skip-preflight", "--only", "pceb"]), runner=r,
                        fetch=lambda u: "down")
        self.assertEqual([c[1] for c in r.calls], ["validate_job.py", "etv_job.py"])
        self.assertNotIn("preflight", s)

    def test_missing_validation_file_is_not_validated(self):
        self.assertFalse(dr4day.read_validation(os.path.join(self.tmp.name, "nowhere")))


if __name__ == "__main__":
    unittest.main()
