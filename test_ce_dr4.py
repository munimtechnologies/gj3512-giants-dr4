"""Unit tests for ce_dr4.py (Amendment 2 applied to the DR4 targets).
Noise-free synthetic epochs only; standard library only.
MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name."""
import math
import os
import tempfile
import unittest

import ce_dr4 as C
from gaiaorb import rvdr4 as R
from test_rvdr4_a1 import B463, B9773, star0, synth


def stub_post(*_a, **_k):
    return dict(status="ok", inc_fold_deg=dict(p02=20.0, p16=21.0, median=22.0, p84=23.0, p98=24.0),
                mass_true_mjup=dict(p02=4.0, p16=4.2, median=4.4, p84=4.6, p98=5.0))


class TestGJ463(unittest.TestCase):
    INC, NODE = 21.8, 147.0

    @classmethod
    def setUpClass(cls):
        pts = synth((B463,), cls.INC, cls.NODE, sigma=0.005)
        cls.res = C.analyse_ce(pts, "gj463", True, "CE_PASS", star=star0("gj463"))

    def test_constrained_primary_and_recovered(self):
        self.assertEqual(self.res["primary_estimator"], "constrained")
        b = self.res["planets"]["b"]
        self.assertEqual(b["verdict"], "TRUE_MASS_MEASURED")
        post = b["printed_omega"]
        self.assertEqual(post["status"], "ok")
        self.assertLess(abs(post["inc_fold_deg"]["median"] - self.INC), 1.0)
        self.assertTrue(post["inc_fold_deg"]["p02"] <= self.INC <= post["inc_fold_deg"]["p98"])
        m_true = B463["msini"] / math.sin(math.radians(self.INC))
        self.assertLess(abs(post["mass_true_mjup"]["median"] / m_true - 1.0), 0.05)
        self.assertNotIn("alternative_omega", b)

    def test_checks(self):
        rows = self.res["planets"]["b"]["checks"]
        self.assertEqual(rows[0]["check"], "AGREES")                 # truth = Morales+2025 here
        self.assertEqual([r["check"] for r in rows[1:]], ["DISAGREES"] * 4)
        self.assertEqual(rows, self.res["planets"]["b"]["ce_checks"])


class TestRules(unittest.TestCase):
    def test_shape_primary_when_calibration_fails(self):
        pts = synth((B463,), 21.8, 147.0, sigma=0.005, n=120)
        res = C.analyse_ce(pts, "gj463", True, "CE_FAIL:C2_low_inclination", star=star0("gj463"), ce_fn=stub_post)
        b = res["planets"]["b"]
        self.assertEqual(res["primary_estimator"], "shape")
        self.assertEqual(b["checks"], b["shape_checks"])
        self.assertEqual(b["ce_checks"][0]["check"], "AGREES")      # 21.8 in [20, 24] (stub)

    def test_gj9773_both_omegas(self):
        pts = synth((B9773,), 22.3, 29.0, n=120)
        res = C.analyse_ce(pts, "gj9773", True, "CE_PASS", star=star0("gj9773"), ce_fn=stub_post)
        b = res["planets"]["b"]
        self.assertIn("alternative_omega", b)
        cc = b["ce_checks"]
        if b["verdict"].startswith("TRUE_MASS_MEASURED"):
            self.assertEqual(cc[0]["check"], "AGREES")              # 22.3 in [20, 24] for both (stub)
        elif cc is not None:
            self.assertEqual({r["check"] for r in cc}, {"NOT_EVALUATED"})

    def test_combine_checks(self):
        a = [dict(ref="x", inc_fold_pred=1.0, mass_pred_mjup=2.0, check="AGREES")]
        d = [dict(ref="x", inc_fold_pred=1.0, mass_pred_mjup=2.0, check="DISAGREES")]
        self.assertEqual(C.combine_checks(a, a)[0]["check"], "AGREES")
        self.assertEqual(C.combine_checks(a, d)[0]["check"], "INCONCLUSIVE_OMEGA_AMBIGUOUS")
        self.assertIsNone(C.combine_checks(a, None))

    def test_inputs(self):
        ci = C.ce_inputs_real(B463, R.STARS["gj463"])
        self.assertAlmostEqual(ci["s"], 0.6212, places=3)
        self.assertAlmostEqual(ci["s_rel"], math.sqrt((0.71 / 36.27) ** 2 + (22.0 / 3459.0) ** 2 +
                                                      (0.146 * 0.015 / (1 - 0.146 ** 2)) ** 2), places=12)
        self.assertEqual(ci["omega_sig"], 5.9)

    def test_read_points_round_trip(self):
        pts = synth((B463,), 21.8, 147.0, n=10)
        d = tempfile.mkdtemp()
        p = os.path.join(d, "points_gj463.csv")
        C.A1.etv_job.write_csv(p, pts)
        back = C.read_points(p)
        self.assertEqual(len(back), 10)
        for a, b in zip(pts, back):
            for k in a:
                self.assertEqual(a[k], b[k])


if __name__ == "__main__":
    unittest.main()
