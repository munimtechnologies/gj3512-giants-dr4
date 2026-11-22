"""Unit tests for epocha2_job.py (Amendment 2 calibration). Noise-free synthetic
epochs only; standard library only.
MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name."""
import math
import random
import unittest

import epocha_job as E
import epocha2_job as E2
from gaiaorb import rvdr4 as R
from test_epocha import TR_RAW, make_rows


def row(cov_shape, cov_ce, inc_true=40.0, P=2.0, verdict="TRUE_MASS_MEASURED", ce_ok=True):
    r = dict(P_yr=P, verdict=verdict, inc_fold_true=inc_true, covered_inc=cov_shape, covered_mass=cov_shape,
             inc_fold_med=inc_true + 3.0, dchi2=100.0)
    if ce_ok:
        r.update(ce_status="ok", ce_covered_inc=cov_ce, ce_covered_mass=cov_ce, ce_inc_med=inc_true + 0.5,
                 ce_branch="omega")
    else:
        r.update(ce_status="bad_covariance")
    return r


class TestGatesCE(unittest.TestCase):
    def test_pass(self):
        rows = [row(True, True) for _ in range(40)] + [row(False, i < 17, inc_true=10.0) for i in range(20)]
        g, out = E2.gates_ce(rows)
        self.assertEqual(out, "CE_PASS")
        self.assertAlmostEqual(g["C2_low_inclination"]["coverage"], 0.85, places=12)
        self.assertEqual(g["C1_overall_inclination"]["n"], 60)

    def test_fail_low(self):
        rows = [row(True, True) for _ in range(40)] + [row(False, i < 15, inc_true=10.0) for i in range(20)]
        _g, out = E2.gates_ce(rows)
        self.assertEqual(out, "CE_FAIL:C2_low_inclination")

    def test_not_finished_counts_as_miss(self):
        rows = [row(True, True) for _ in range(40)] + [row(False, True, inc_true=10.0, ce_ok=(i >= 5))
                                                       for i in range(20)]
        g, out = E2.gates_ce(rows)
        self.assertAlmostEqual(g["C2_low_inclination"]["coverage"], 0.75, places=12)
        self.assertEqual(out, "CE_FAIL:C2_low_inclination")

    def test_underpowered(self):
        rows = [row(True, True) for _ in range(40)] + [row(False, True, inc_true=10.0) for _ in range(10)]
        _g, out = E2.gates_ce(rows)
        self.assertEqual(out, "CE_UNDERPOWERED")

    def test_long_periods_and_unmeasured_excluded(self):
        rows = [row(True, True) for _ in range(40)] + [row(False, False, P=15.0) for _ in range(40)] + \
               [row(False, False, verdict="NOT_DETECTED") for _ in range(40)]
        c = E2.compare(rows)
        self.assertEqual(c["n"], 40)
        self.assertEqual(c["ce_inc_coverage"], 1.0)
        self.assertEqual(c["shape_n_high"], 40)
        self.assertAlmostEqual(c["shape_median_signed_dinc"], 3.0, places=12)
        self.assertAlmostEqual(c["ce_median_signed_dinc"], 0.5, places=12)


class TestInputs(unittest.TestCase):
    def test_oracle_and_perturbed(self):
        tr = E.truth_from_raw(TR_RAW)
        planet, star = E.planet_star_from_truth(tr, "oracle")
        ci = E2.ce_inputs(planet, star, "oracle", random.Random(1))
        si = math.sin(math.radians(tr["inc"]))
        self.assertAlmostEqual(ci["s"], tr["alpha_mas"] * si, places=9)
        self.assertEqual(ci["omega"], planet["omega"])
        self.assertEqual(ci["omega_sig"], 0.0)
        g = random.Random(7)
        f, dom = 1.0 + 0.03 * g.gauss(0.0, 1.0), 3.0 * g.gauss(0.0, 1.0)
        cp = E2.ce_inputs(planet, star, "perturbed", random.Random(7))
        self.assertAlmostEqual(cp["s"], tr["alpha_mas"] * si * f, places=9)
        self.assertAlmostEqual(cp["msini"], planet["msini"] * f, places=12)
        self.assertAlmostEqual(cp["omega"], planet["omega"] + dom, places=12)


class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tr = E.truth_from_raw(TR_RAW)
        cls.pts = E.points_from_rows(make_rows(cls.tr))

    def test_v1_fields_unchanged_and_ce_recovers(self):
        r1 = E.evaluate(self.pts, self.tr, "oracle")
        r2 = E2.evaluate2(self.pts, self.tr, "oracle", random.Random(1), random.Random(2))
        for k, v in r1.items():
            self.assertEqual(r2[k], v, k)
        self.assertEqual(r2["ce_status"], "ok")
        self.assertTrue(r2["ce_covered_inc"])
        self.assertTrue(r2["ce_covered_mass"])
        self.assertLess(abs(r2["ce_inc_med"] - 50.0), 1.5)
        self.assertEqual(r2["ce_branch"], "omega")

    def test_perturbed_v1_fields_unchanged(self):
        stub = lambda *a, **k: dict(status="stub")   # noqa: E731  (skip the slow grid here)
        r1 = E.evaluate(self.pts, self.tr, "perturbed", random.Random(4))
        r2 = E2.evaluate2(self.pts, self.tr, "perturbed", random.Random(4), random.Random(5), ce_fn=stub)
        for k, v in r1.items():
            self.assertEqual(r2[k], v, k)
        if r2["verdict"].startswith("TRUE_MASS_MEASURED"):
            self.assertEqual(r2["ce_status"], "stub")

    def test_unmeasured_not_run(self):
        pts0 = E.points_from_rows(make_rows(self.tr, with_orbit=False))
        r = E2.evaluate2(pts0, self.tr, "oracle", random.Random(1), random.Random(2))
        self.assertFalse(r["verdict"].startswith("TRUE_MASS_MEASURED"))
        self.assertEqual(r["ce_status"], "not_run")


class TestCsv(unittest.TestCase):
    def test_union_of_keys(self):
        import os
        import tempfile
        d = tempfile.mkdtemp()
        p = os.path.join(d, "x.csv")
        E2.write_csv(p, [dict(a=1, ce_status="not_run"), dict(a=2, ce_status="ok", ce_inc_med=3.0)])
        with open(p) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[0], "a,ce_status,ce_inc_med")
        self.assertEqual(lines[1], "1,not_run,")


if __name__ == "__main__":
    unittest.main()
