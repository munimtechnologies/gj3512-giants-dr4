"""Unit tests for gaiaorb.rvdr4 (GJ 3512 giants in Gaia DR4). Standard library only.
MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name."""
import math
import random
import unittest

from gaiaorb import kepler, model, rvdr4 as R

B, C = R.PLANETS
STAR0 = dict(R.STARS["gj3512"], rv_kms=0.0)   # synthetic data have no perspective term
INC_TRUE, OMEGA_NODE_TRUE = 36.0, 40.0


def synth(n=320, sigma=0.05, seed=7, planets=R.PLANETS, inc=INC_TRUE, node=OMEGA_NODE_TRUE, noise=False):
    rng = random.Random(seed)
    pts = []
    for i in range(n):
        t = -2.4 + 5.5 * i / (n - 1)
        psi = rng.uniform(0.0, 2.0 * math.pi)
        s, c = math.sin(psi), math.cos(psi)
        f = 0.9 * math.sin(2.0 * math.pi * t + 0.4) * s + 0.6 * math.cos(2.0 * math.pi * t + 0.4) * c
        w = 1.0 * s - 2.0 * c + 105.0 * f - 260.0 * t * s - 1279.0 * t * c
        pts.append(dict(t=t, s=s, c=c, f=f, fr=0.0, sig=sigma, w=w))
    for pl in planets:
        a0 = R.a0_min_mas(pl) / math.sin(math.radians(inc))
        ti = kepler.campbell_to_ti(a0, inc, pl["omega"], node)
        sgn = model.orbit_signal(pts, ti, pl["P"], pl["e"], R.tp_days(pl))
        for p, v in zip(pts, sgn):
            p["w"] += v + (rng.gauss(0.0, sigma) if noise else 0.0)
    return pts


class TestPredictions(unittest.TestCase):
    def test_tp_circular(self):
        self.assertAlmostEqual(R.tp_from_tc(1000.0, 100.0, 0.0, 90.0), 1000.0, places=9)
        self.assertAlmostEqual(R.tp_from_tc(1000.0, 100.0, 0.0, 0.0), 975.0, places=9)

    def test_tp_round_trip_b(self):
        tp = R.tp_from_tc(B["Tc"], B["P"], B["e"], B["omega"])
        f = R.true_anomaly(B["Tc"], B["P"], B["e"], tp)
        want = math.pi / 2.0 - math.radians(B["omega"])
        d = (f - want + math.pi) % (2.0 * math.pi) - math.pi
        self.assertLess(abs(d), 1e-9)
        self.assertGreater(tp, B["Tc"])            # conjunction comes 7.8 d before periastron
        self.assertAlmostEqual(tp - B["Tc"], 7.84, delta=0.02)

    def test_minimum_wobbles(self):
        self.assertAlmostEqual(R.a0_min_mas(B), 0.1266, places=3)   # Morales+2019: ~120 microarcsec
        self.assertAlmostEqual(R.a0_min_mas(C), 0.6359, places=3)

    def test_min_wobble_gives_msini(self):
        self.assertAlmostEqual(R.mass_from_a0(R.a0_min_mas(B), B), 0.461, delta=0.01)
        self.assertAlmostEqual(R.mass_from_a0(R.a0_min_mas(C), C), 0.453, delta=0.01)

    def test_cholesky(self):
        M = [[4.0, 2.0], [2.0, 3.0]]
        L = R.cholesky(M)
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(sum(L[i][k] * L[j][k] for k in range(2)), M[i][j], places=12)


class TestRecovery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = R.run(synth(), star=STAR0, log=lambda *a: None)

    def test_status(self):
        self.assertEqual(self.res["status"], "ok")
        self.assertEqual(self.res["n_transits"], 320)
        self.assertEqual(self.res["n_clipped"], 0)

    def test_b_recovered(self):
        b = self.res["planets"]["b"]
        a0_true = R.a0_min_mas(B) / math.sin(math.radians(INC_TRUE))
        self.assertAlmostEqual(b["point"]["a0_mas"], a0_true, delta=1e-3)
        self.assertAlmostEqual(b["point"]["inc_deg"], INC_TRUE, delta=0.5)
        self.assertAlmostEqual(b["point"]["Omega_deg"], OMEGA_NODE_TRUE, delta=0.5)
        self.assertLess(b["omega_diff_mod180_deg"], 0.5)
        self.assertGreater(b["dchi2"], 100.0)
        self.assertGreater(b["expected_dchi2_min_mass"], R.MISSING_LAMBDA)
        self.assertLess(b["a0_sini_z"], 1.0)
        self.assertAlmostEqual(b["posterior"]["mass_true_mjup"]["median"], 0.461 / math.sin(math.radians(INC_TRUE)),
                               delta=0.1)
        self.assertEqual(R.verdict(b, True), "TRUE_MASS_MEASURED")
        self.assertEqual(R.verdict(b, False), "INCONCLUSIVE")

    def test_c_recovered(self):
        c = self.res["planets"]["c"]
        a0_true = R.a0_min_mas(C) / math.sin(math.radians(INC_TRUE))
        self.assertAlmostEqual(c["point"]["a0_mas"], a0_true, delta=5e-3)
        self.assertAlmostEqual(c["point"]["inc_deg"], INC_TRUE, delta=2.0)
        self.assertGreater(c["dchi2"], R.DET_DCHI2)


class TestNull(unittest.TestCase):
    def test_no_planets_not_detected(self):
        res = R.run(synth(planets=()), star=STAR0, log=lambda *a: None)
        b = res["planets"]["b"]
        self.assertLess(b["dchi2"], 1e-6)
        self.assertEqual(R.verdict(b, True), "RV_ORBIT_NOT_SEEN")   # min-mass orbit would have been obvious


class TestVerdict(unittest.TestCase):
    def mk(self, d, lam, z, p16):
        return dict(dchi2=d, expected_dchi2_min_mass=lam, a0_sini_z=z,
                    posterior=dict(mass_true_mjup=dict(p16=p16)))

    def test_branches(self):
        self.assertEqual(R.verdict(self.mk(10.0, 20.0, 0.0, 1.0), True), "NOT_DETECTED")
        self.assertEqual(R.verdict(self.mk(10.0, 80.0, 0.0, 1.0), True), "RV_ORBIT_NOT_SEEN")
        self.assertEqual(R.verdict(self.mk(90.0, 80.0, 4.0, 1.0), True), "ORBIT_INCONSISTENT_WITH_RV")
        self.assertEqual(R.verdict(self.mk(90.0, 80.0, 1.0, 20.0), True), "TRUE_MASS_MEASURED_BROWN_DWARF")
        self.assertEqual(R.verdict(self.mk(90.0, 80.0, 1.0, 0.7), True), "TRUE_MASS_MEASURED")
        self.assertEqual(R.verdict(self.mk(90.0, 80.0, 1.0, 0.7), True, status="error"), "INCONCLUSIVE")


if __name__ == "__main__":
    unittest.main()
