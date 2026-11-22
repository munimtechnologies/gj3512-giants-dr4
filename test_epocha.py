"""Unit tests for epocha_job.py (Epochalypse calibration of gaiaorb/rvdr4.py).
Noise-free synthetic data only; standard library only.
MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name."""
import math
import random
import unittest

from gaiaorb import kepler, model, rvdr4 as R
import epocha_job as E

TR_RAW = dict(sid=123456789, plx=50.0, mstar=0.4, sma=2.0, ecc=0.2, inc=50.0, M0=30.0, mpl=3.0,
              omega=70.0, Omega=100.0)


def make_rows(tr, n=120, sigma=0.05, seed=3, with_orbit=True):
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        t = -2.4 + 5.5 * i / (n - 1)
        psi = rng.uniform(-math.pi, math.pi)
        s, c = math.sin(psi), math.cos(psi)
        f = 0.8 * math.sin(2.0 * math.pi * t) * s + 0.5 * math.cos(2.0 * math.pi * t) * c
        w = 0.3 * s - 0.7 * c + tr["plx"] * f + 120.0 * t * s - 300.0 * t * c
        rows.append(dict(t=E.T0_JD + 365.25 * t, psi=psi, f=f, w=w, sig=sigma))
    if with_orbit:
        planet, _star = E.planet_star_from_truth(tr, "oracle")
        pts = E.points_from_rows(rows)
        ti = kepler.campbell_to_ti(tr["alpha_mas"], tr["inc"], planet["omega"], tr["Omega"])
        sig = model.orbit_signal(pts, ti, planet["P"], planet["e"], R.tp_days(planet))
        for r, p, v in zip(sorted(rows, key=lambda r: r["t"]), pts, sig):
            r["w"] += v
    return rows


class TestPure(unittest.TestCase):
    def test_tc_tp_round_trip(self):
        for e, om in ((0.0, 90.0), (0.3, 40.0), (0.6, 250.0)):
            tc = E.tc_from_tp(2458000.0, 1000.0, e, om)
            self.assertAlmostEqual(R.tp_from_tc(tc, 1000.0, e, om), 2458000.0, places=6)

    def test_truth_derivations(self):
        tr = E.truth_from_raw(dict(sid=1, plx=50.0, mstar=0.5, sma=2.0, ecc=0.1, inc=30.0, M0=0.0,
                                   mpl=0.01 * E.MJ_PER_MSUN, omega=0.0, Omega=0.0))
        self.assertAlmostEqual(tr["P_yr"], math.sqrt(8.0 / 0.51), places=9)
        self.assertAlmostEqual(tr["alpha_mas"], 0.01 / 0.51 * 2.0 * 50.0, places=9)
        self.assertEqual(tr["sid"], "1")

    def test_planet_matches_truth(self):
        tr = E.truth_from_raw(TR_RAW)
        planet, star = E.planet_star_from_truth(tr, "oracle")
        si = math.sin(math.radians(tr["inc"]))
        self.assertAlmostEqual(R.a0_min_mas(planet, star), tr["alpha_mas"] * si, places=9)
        tp_want = E.T0_JD - math.radians(tr["M0"]) * tr["P_yr"] * 365.25 / (2.0 * math.pi)
        self.assertAlmostEqual(R.tp_days(planet) + R.J2017_5_JD, tp_want, places=5)
        self.assertAlmostEqual(planet["msini"], 3.0 * si, places=12)
        self.assertEqual(star["rv_kms"], 0.0)

    def test_resolve_and_sid(self):
        km, miss = E.resolve(["gaia_source_id", "sma_1", "ecc"], E.TRUTH_COLS)
        self.assertEqual(km["sma"], "sma_1")
        self.assertEqual(km["ecc"], "ecc")
        self.assertIn("inc", miss)
        self.assertEqual(E.norm_sid(123.0), "123")
        self.assertEqual(E.norm_sid(" 456 "), "456")

    def test_points(self):
        pts = E.points_from_rows([dict(t=E.T0_JD + 365.25, psi=math.pi / 2, f=0.3, w=2.0, sig=0.1),
                                  dict(t=E.T0_JD, psi=0.0, f=0.1, w=1.0, sig=0.1)])
        self.assertAlmostEqual(pts[0]["t"], 0.0, places=12)
        self.assertAlmostEqual(pts[1]["t"], 1.0, places=12)
        self.assertAlmostEqual(pts[1]["s"], 1.0, places=12)
        self.assertAlmostEqual(pts[0]["c"], 1.0, places=12)

    def test_selectable(self):
        tr = E.truth_from_raw(TR_RAW)
        self.assertTrue(E.selectable(tr))
        self.assertFalse(E.selectable(dict(tr, ecc=0.9)))
        self.assertFalse(E.selectable(dict(tr, alpha_mas=0.01)))


class TestGates(unittest.TestCase):
    def rows(self, n_cov, n=40):
        return [dict(P_yr=2.0, verdict="TRUE_MASS_MEASURED", dchi2=100.0, covered_inc=(i < n_cov),
                     covered_mass=True, inc_fold_med=30.0, inc_fold_true=31.0) for i in range(n)]

    def nulls(self, n, n_det):
        return [dict(detected=(i < n_det)) for i in range(n)]

    def test_pass(self):
        g, out = E.gates(self.rows(36), self.nulls(150, 1))
        self.assertEqual(out, "EPOCHALYPSE_PASS")
        self.assertAlmostEqual(g["G2_inclination_coverage"]["coverage"], 0.9, places=12)

    def test_fail_coverage(self):
        _g, out = E.gates(self.rows(30), self.nulls(150, 1))
        self.assertEqual(out, "EPOCHALYPSE_FAIL:G2_inclination_coverage")

    def test_fail_null(self):
        _g, out = E.gates(self.rows(36), self.nulls(150, 2))   # 1.33% > 1%
        self.assertEqual(out, "EPOCHALYPSE_FAIL:G1_null")

    def test_underpowered(self):
        _g, out = E.gates(self.rows(36), self.nulls(50, 0))
        self.assertEqual(out, "EPOCHALYPSE_UNDERPOWERED")

    def test_long_periods_excluded(self):
        rows = self.rows(36) + [dict(r, P_yr=15.0, covered_inc=False) for r in self.rows(0)]
        s = E.summarize_rows(rows)
        self.assertEqual(s["n"], 40)
        self.assertEqual(E.summarize_rows(rows, p_max=1e9)["n"], 80)


class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tr = E.truth_from_raw(TR_RAW)
        cls.pts = E.points_from_rows(make_rows(cls.tr))

    def test_oracle_recovers(self):
        r = E.evaluate(self.pts, self.tr, "oracle")
        self.assertEqual(r["verdict"], "TRUE_MASS_MEASURED")
        self.assertAlmostEqual(r["a0_point"], self.tr["alpha_mas"], delta=1e-3)
        self.assertTrue(r["covered_inc"])
        self.assertTrue(r["covered_mass"])
        self.assertAlmostEqual(r["inc_fold_true"], 50.0, places=9)

    def test_perturbed_runs(self):
        r = E.evaluate(self.pts, self.tr, "perturbed", random.Random(1))
        self.assertIn(r["verdict"], ("TRUE_MASS_MEASURED", "TRUE_MASS_MEASURED_BROWN_DWARF",
                                     "ORBIT_INCONSISTENT_WITH_RV", "RV_ORBIT_NOT_SEEN", "NOT_DETECTED"))

    def test_null(self):
        pts0 = E.points_from_rows(make_rows(self.tr, with_orbit=False))
        r = E.evaluate_null(pts0, self.tr)
        self.assertLess(r["dchi2"], 1e-6)
        self.assertFalse(r["detected"])

    def test_group_rows(self):
        cols = dict(gaia_source_id=["1", "2", "1"], obs_time_tcb=[1.0, 2.0, 3.0])
        by = E.group_rows(cols, dict(sid="gaia_source_id", t="obs_time_tcb"))
        self.assertEqual(sorted(by), ["1", "2"])
        self.assertEqual([r["t"] for r in by["1"]], [1.0, 3.0])


if __name__ == "__main__":
    unittest.main()
