"""Unit tests for Amendment 1 of gj3512-giants-dr4 (GJ 317 b/c, GJ 463 b, GJ 9773 b).
Noise-free synthetic epochs only; standard library only.
MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name."""
import math
import random
import unittest

from gaiaorb import kepler, model, rvdr4 as R
import rv_job_a1 as J

B317, C317 = R.SYSTEMS_A1["gj317"]
B463, = R.SYSTEMS_A1["gj463"]
B9773, = R.SYSTEMS_A1["gj9773"]
NEW = (B317, C317, B463, B9773)
QUIET = dict(log=lambda *a: None)


def star0(key):
    return dict(R.STARS[key], rv_kms=0.0)   # synthetic data have no perspective term


def synth(planets, inc, node, n=320, sigma=0.05, seed=11):
    rng = random.Random(seed)
    pts = []
    for i in range(n):
        t = -2.4 + 5.5 * i / (n - 1)
        psi = rng.uniform(0.0, 2.0 * math.pi)
        s, c = math.sin(psi), math.cos(psi)
        f = 0.9 * math.sin(2.0 * math.pi * t + 0.4) * s + 0.6 * math.cos(2.0 * math.pi * t + 0.4) * c
        w = 1.0 * s - 2.0 * c + 60.0 * f - 460.0 * t * s + 805.0 * t * c
        pts.append(dict(t=t, s=s, c=c, f=f, fr=0.0, sig=sigma, w=w))
    for pl in planets:
        a0 = R.a0_min_mas(pl) / math.sin(math.radians(inc))
        ti = kepler.campbell_to_ti(a0, inc, pl["omega"], node)
        for p, v in zip(pts, model.orbit_signal(pts, ti, pl["P"], pl["e"], R.tp_days(pl))):
            p["w"] += v
    return pts


class TestInputs(unittest.TestCase):
    def test_minimum_wobbles(self):
        # hand values: K sqrt(1-e^2) P / 2pi / au x parallax
        self.assertAlmostEqual(R.a0_min_mas(B317), 0.2997, places=3)
        self.assertAlmostEqual(R.a0_min_mas(C317), 1.1269, places=3)
        self.assertAlmostEqual(R.a0_min_mas(B463), 0.6212, places=3)
        self.assertAlmostEqual(R.a0_min_mas(B9773), 0.2527, places=3)

    def test_min_wobble_gives_msini(self):
        for pl in NEW:
            self.assertAlmostEqual(R.mass_from_a0(R.a0_min_mas(pl), pl), pl["msini"], delta=0.03)

    def test_tp_round_trip(self):
        for pl in NEW:
            tp = R.tp_from_tc(pl["Tc"], pl["P"], pl["e"], pl["omega"])
            f = R.true_anomaly(pl["Tc"], pl["P"], pl["e"], tp)
            d = (f - (math.pi / 2.0 - math.radians(pl["omega"])) + math.pi) % (2.0 * math.pi) - math.pi
            self.assertLess(abs(d), 1e-9)

    def test_stars_and_predictions(self):
        for key in ("gj317", "gj463", "gj9773"):
            self.assertIn(key, R.STARS)
            for pl in R.SYSTEMS_A1[key]:
                self.assertEqual(pl["star"], key)
                self.assertIn((key, pl["key"]), R.PREDICTIONS)
        self.assertEqual(R.OMEGA_ALTERNATIVE[("gj9773", "b")], 81.0)
        self.assertEqual(len(R.PREDICTIONS[("gj463", "b")]), 5)

    def test_gj3512_unchanged(self):
        self.assertEqual(len(R.PLANETS), 2)
        self.assertAlmostEqual(R.a0_min_mas(R.PLANETS[0]), 0.1266, places=3)
        self.assertIs(J.system_planets("gj3512"), R.PLANETS)


class TestGeometry(unittest.TestCase):
    def test_mutual_inclination(self):
        self.assertAlmostEqual(R.mutual_inc_deg(30.0, 10.0, 30.0, 10.0), 0.0, places=5)
        self.assertAlmostEqual(R.mutual_inc_deg(90.0, 0.0, 90.0, 30.0), 30.0, places=9)
        self.assertAlmostEqual(R.mutual_inc_deg(30.0, 0.0, 30.0, 180.0), 60.0, places=9)
        self.assertAlmostEqual(R.mutual_inc_deg(0.0, 0.0, 90.0, 50.0), 90.0, places=9)

    def test_fold_and_circ(self):
        self.assertEqual(R.fold_inc(152.0), 28.0)
        self.assertEqual(R.fold_inc(21.8), 21.8)
        self.assertAlmostEqual(R.circ_dist_deg(350.0, 10.0), 20.0, places=9)
        self.assertAlmostEqual(R.circ_dist_deg(10.0, 200.0), 170.0, places=9)


class TestRules(unittest.TestCase):
    def test_combine_omega(self):
        self.assertEqual(R.combine_omega_verdicts("TRUE_MASS_MEASURED", "TRUE_MASS_MEASURED"), "TRUE_MASS_MEASURED")
        self.assertEqual(R.combine_omega_verdicts("TRUE_MASS_MEASURED", "ORBIT_INCONSISTENT_WITH_RV"),
                         "INCONCLUSIVE_OMEGA_AMBIGUOUS")

    def test_prediction_checks(self):
        pres = dict(posterior=dict(inc_fold_deg=dict(p02=20.0, p16=21.0, median=22.0, p84=23.0, p98=24.0),
                                   mass_true_mjup=dict(p02=4.0, p16=4.2, median=4.4, p84=4.6, p98=5.0)))
        preds = R.PREDICTIONS[("gj463", "b")]
        rows = R.prediction_checks(pres, "TRUE_MASS_MEASURED", preds)
        self.assertEqual(rows[0]["check"], "AGREES")                 # Morales: 21.8
        self.assertAlmostEqual(rows[0]["inc_offset_sigma"], -0.2, places=9)
        self.assertTrue(rows[0]["mass_in_dr4_95"])
        self.assertEqual([r["check"] for r in rows[1:]], ["DISAGREES"] * 4)
        self.assertFalse(rows[3]["mass_in_dr4_95"])                  # 3.3 < 4.0
        rows = R.prediction_checks(pres, "NOT_DETECTED", preds)
        self.assertEqual({r["check"] for r in rows}, {"NOT_EVALUATED"})


class TestGJ317(unittest.TestCase):
    INC, NODE = 34.0, 40.0

    @classmethod
    def setUpClass(cls):
        cls.res = R.run(synth((B317, C317), cls.INC, cls.NODE), planets=(B317, C317), star=star0("gj317"), **QUIET)

    def test_b(self):
        b = self.res["planets"]["b"]
        self.assertEqual(self.res["n_clipped"], 0)
        self.assertAlmostEqual(b["point"]["a0_mas"], R.a0_min_mas(B317) / math.sin(math.radians(self.INC)),
                               delta=1e-3)
        self.assertAlmostEqual(b["point"]["inc_deg"], self.INC, delta=0.5)
        self.assertAlmostEqual(b["point"]["Omega_deg"], self.NODE, delta=0.5)
        self.assertLess(b["omega_diff_mod180_deg"], 0.5)
        self.assertFalse(b["rv_branch_flipped"])
        self.assertGreater(b["dchi2"], 100.0)
        self.assertEqual(R.verdict(b, True), "TRUE_MASS_MEASURED")
        self.assertEqual(b["excess_noise_mass_mjup"], [])

    def test_c(self):
        c = self.res["planets"]["c"]
        a0 = R.a0_min_mas(C317) / math.sin(math.radians(self.INC))
        self.assertAlmostEqual(c["point"]["a0_mas"], a0, delta=0.01 * a0)
        self.assertAlmostEqual(c["point"]["inc_deg"], self.INC, delta=2.0)
        self.assertGreater(c["dchi2"], R.DET_DCHI2)

    def test_coplanar(self):
        m = self.res["mutual_inclination"]
        self.assertLess(m["point_same_deg"], 2.0)
        self.assertGreater(m["point_flipped_deg"], 60.0)          # 2 x 34 deg when Omega is flipped
        self.assertEqual(m["rv_branch"], "same")
        self.assertGreater(m["n"], 3000)


class TestGJ463(unittest.TestCase):
    INC, NODE = 21.8, 147.0

    @classmethod
    def setUpClass(cls):
        pts = synth((B463,), cls.INC, cls.NODE, sigma=0.005)
        cls.res = J.analyse(pts, "gj463", True, star=star0("gj463"), **QUIET)

    def test_recovered(self):
        b = self.res["planets"]["b"]
        self.assertAlmostEqual(b["point"]["inc_deg"], self.INC, delta=0.5)
        self.assertAlmostEqual(b["point"]["Omega_deg"], self.NODE, delta=0.5)
        self.assertAlmostEqual(b["posterior"]["mass_true_mjup"]["median"],
                               B463["msini"] / math.sin(math.radians(self.INC)), delta=0.3)
        self.assertEqual(self.res["verdicts"]["b"], "TRUE_MASS_MEASURED")
        self.assertNotIn("mutual_inclination", self.res)
        self.assertNotIn("omega_alternative", self.res)

    def test_discriminates_published_predictions(self):
        rows = self.res["prediction_checks"]["b"]
        self.assertEqual(rows[0]["check"], "AGREES")                 # truth = Morales+2025 here
        self.assertEqual([r["check"] for r in rows[1:]], ["DISAGREES"] * 4)


class TestGJ9773OmegaAmbiguity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pts = synth((B9773,), 22.3, 29.0)
        cls.res = J.analyse(pts, "gj9773", True, star=star0("gj9773"), **QUIET)

    def test_both_fits_and_combined_verdict(self):
        alt = self.res["omega_alternative"]
        prim = R.verdict(self.res["planets"]["b"], True)
        self.assertEqual(prim, "TRUE_MASS_MEASURED")
        self.assertIn("verdict_b", alt)
        want = prim if alt["verdict_b"] == prim else "INCONCLUSIVE_OMEGA_AMBIGUOUS"
        self.assertEqual(self.res["verdicts"]["b"], want)
        self.assertGreaterEqual(alt["jitter_mas"], 0.0)

    def test_unvalidated_is_inconclusive(self):
        pts = synth((B9773,), 22.3, 29.0, n=120)
        res = J.analyse(pts, "gj9773", False, star=star0("gj9773"), **QUIET)
        self.assertEqual(res["verdicts"]["b"], "INCONCLUSIVE")


if __name__ == "__main__":
    unittest.main()
