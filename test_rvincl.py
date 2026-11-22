"""Unit tests for gaiaorb/rvincl.py (Amendment 2: RV-amplitude-constrained inclination).
Noise-free inputs and exact integrals only; standard library only.
MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name."""
import math
import random
import unittest

from gaiaorb import kepler, rvincl as V


def i0_series(x, n=80):
    """Modified Bessel function I0 by its power series."""
    s, term = 0.0, 1.0
    for k in range(n):
        if k > 0:
            term *= (x / 2.0) ** 2 / (k * k)
        s += term
    return s


def diag_cov(sig):
    return [[sig * sig if r == c else 0.0 for c in range(4)] for r in range(4)]


def vm(k, shift=0.0):
    """Quadratic forms giving, with mu = 1 and eps = 0, chi2(W) = k (1 - cos(W - shift))."""
    return (0.0, 0.0, 0.0, 0.5 * k * math.cos(shift), 0.5 * k * math.sin(shift), k)


def rand_spd(rng):
    M = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(4)]
    return [[sum(M[r][k] * M[c][k] for k in range(4)) + (1.0 if r == c else 0.0) for c in range(4)]
            for r in range(4)]


class TestIntegral(unittest.TestCase):
    def test_flat(self):
        li, e = V.log_int_w((0.0,) * 6, 1.0, 0.0)
        self.assertAlmostEqual(li, math.log(2.0 * math.pi), places=12)
        self.assertEqual(e, 0.0)

    def test_von_mises_wide(self):
        k = 10.0
        li, _e = V.log_int_w(vm(k), 1.0, 0.0)
        exact = math.log(2.0 * math.pi) - k / 2.0 + math.log(i0_series(k / 2.0))
        self.assertAlmostEqual(li, exact, places=9)

    def test_von_mises_narrow_refined(self):
        k = 1e6                        # peak width ~1.4e-3 rad: fine local grid
        for shift in (0.0, 0.01):      # peak on and off the coarse grid
            li, e = V.log_int_w(vm(k, shift), 1.0, 0.0)
            exact = 0.5 * math.log(4.0 * math.pi / k) + math.log(1.0 + 1.0 / (4.0 * k))
            self.assertAlmostEqual(li, exact, places=7)
            self.assertLess(e, 0.01)

    def test_von_mises_very_narrow_laplace(self):
        k = 1e9                        # narrower than the fine grid: Laplace
        li, _e = V.log_int_w(vm(k, 0.01), 1.0, 0.0)
        self.assertAlmostEqual(li, 0.5 * math.log(4.0 * math.pi / k), places=6)

    def test_two_peaks(self):
        # chi2 = k (1 - cos 2W): two equal narrow peaks at 0 and pi
        k = 1e6
        li, _e = V.log_int_w((-k, k, 0.0, 0.0, 0.0, k), 1.0, 0.0)
        exact = math.log(2.0) + 0.5 * math.log(4.0 * math.pi / (4.0 * k)) + math.log(1.0 + 1.0 / (4.0 * k))
        self.assertAlmostEqual(li, exact, places=6)

    def test_a_integral_closed_form(self):
        A2, A1, bQb, mu = 2.0, 1.3, 1.5, 0.8
        for eps in (0.09, 0.0004):
            tau = math.sqrt(eps)
            n, lo, hi = 20001, mu - 12.0 * tau, mu + 12.0 * tau
            h = (hi - lo) / (n - 1)
            tot = 0.0
            for j in range(n):
                a = lo + j * h
                w = 0.5 if j in (0, n - 1) else 1.0
                tot += w * math.exp(-0.5 * (a - mu) ** 2 / eps) / math.sqrt(2.0 * math.pi * eps) * \
                    math.exp(-0.5 * (A2 * a * a - 2.0 * A1 * a + bQb))
            got, _E = V.a_integral(A2, A1, bQb, mu, eps)
            self.assertAlmostEqual(got, math.log(tot * h), places=9)
        got, E = V.a_integral(A2, A1, bQb, mu, 0.0)
        self.assertAlmostEqual(E, A2 * mu * mu - 2.0 * A1 * mu + bQb, places=12)
        self.assertAlmostEqual(got, -0.5 * E, places=12)


class TestQuadraticForms(unittest.TestCase):
    def test_chi2_matches_campbell(self):
        rng = random.Random(5)
        for _ in range(20):
            Q = rand_spd(rng)
            b = [rng.gauss(0, 1) for _ in range(4)]
            a, inc, om, W = rng.uniform(0.1, 3), rng.uniform(0, 180), rng.uniform(0, 360), rng.uniform(0, 360)
            ti = kepler.campbell_to_ti(a, inc, om, W)
            d = [ti[k] - b[k] for k in range(4)]
            direct = sum(d[r] * Q[r][c] * d[c] for r in range(4) for c in range(4))
            bQb = sum(b[r] * Q[r][c] * b[c] for r in range(4) for c in range(4))
            got = V.chi2_direct(V.qforms(Q, b, bQb, inc, om), a, math.radians(W))
            self.assertAlmostEqual(got, direct, delta=1e-9 * max(1.0, abs(direct)))

    def test_terms_match_a_integral(self):
        rng = random.Random(9)
        Q = rand_spd(rng)
        b = [rng.gauss(0, 1) for _ in range(4)]
        bQb = sum(b[r] * Q[r][c] * b[c] for r in range(4) for c in range(4))
        qs = V.qforms(Q, b, bQb, 37.0, 120.0)
        W = 1.1
        lt, es = V._terms(qs, 0.9, 0.04, V._trig([W]))
        c, s = math.cos(W), math.sin(W)
        A2 = c * c * qs[0] + 2.0 * c * s * qs[2] + s * s * qs[1]
        A1 = c * qs[3] + s * qs[4]
        want, E = V.a_integral(A2, A1, bQb, 0.9, 0.04)
        self.assertAlmostEqual(lt[0], want, places=12)
        self.assertAlmostEqual(es[0], E, places=12)


class TestPosterior(unittest.TestCase):
    def run_case(self, i0, s=0.1, om=50.0, W0=30.0, sig=0.005, s_rel=0.02, om_sig=0.0, gen_om=None,
                 msini=2.0, **kw):
        a = s / math.sin(math.radians(i0))
        ti = kepler.campbell_to_ti(a, i0, om if gen_om is None else gen_om, W0)
        return V.posterior(ti, diag_cov(sig), s, s_rel, om, om_sig, msini, 0.0, n_draw=2000, **kw)

    def test_face_on_recovered(self):
        r = self.run_case(8.0)
        self.assertEqual(r["status"], "ok")
        self.assertLess(abs(r["inc_fold_deg"]["median"] - 8.0), 1.0)
        self.assertTrue(V.covered(r["inc_fold_deg"], 8.0))

    def test_intermediate_recovered_and_mass(self):
        r = self.run_case(25.0)
        self.assertLess(abs(r["inc_fold_deg"]["median"] - 25.0), 1.0)
        self.assertTrue(V.covered(r["inc_fold_deg"], 25.0))
        m_true = 2.0 / math.sin(math.radians(25.0))
        self.assertLess(abs(r["mass_true_mjup"]["median"] / m_true - 1.0), 0.03)
        self.assertTrue(V.covered(r["mass_true_mjup"], m_true))

    def test_inclined_and_branch_evidence(self):
        r = self.run_case(60.0, s=0.2, s_rel=0.01)
        self.assertLess(abs(r["inc_fold_deg"]["median"] - 60.0), 1.5)
        self.assertTrue(V.covered(r["inc_fold_deg"], 60.0))
        ev = r["log_evidence_branch"]
        self.assertEqual(r["preferred_branch"], "omega")
        self.assertGreater(ev["omega"] - ev["minus_omega"], 10.0)
        self.assertLess(r["chi2_min_constrained"], 0.05)

    def test_flipped_convention_still_recovered(self):
        r = self.run_case(60.0, s=0.2, s_rel=0.01, gen_om=-50.0)
        self.assertEqual(r["preferred_branch"], "minus_omega")
        self.assertTrue(V.covered(r["inc_fold_deg"], 60.0))

    def test_omega_nodes_run(self):
        r = self.run_case(40.0, om_sig=3.0, n_coarse=60, n_fine=80)
        self.assertEqual(r["status"], "ok")
        self.assertTrue(V.covered(r["inc_fold_deg"], 40.0))

    def test_bad_covariance(self):
        r = V.posterior([0.1, 0.0, 0.0, 0.1], [[0.0] * 4 for _ in range(4)], 0.1, 0.02, 0.0, 0.0, 1.0, 0.0)
        self.assertEqual(r["status"], "bad_covariance")


class TestHelpers(unittest.TestCase):
    def test_omega_nodes(self):
        self.assertEqual(V.omega_nodes(0.0), ((0.0, 1.0),))
        nd = V.omega_nodes(3.0)
        self.assertEqual(len(nd), 13)                       # spacing 1.5 deg over +-9 deg
        self.assertAlmostEqual(nd[0][0], -9.0, places=12)
        self.assertAlmostEqual(sum(w for _x, w in nd), 1.0, places=12)
        self.assertAlmostEqual(sum(w * x for x, w in nd), 0.0, places=12)
        var = sum(w * x * x for x, w in nd)
        self.assertGreater(var, 0.95 * 9.0)
        self.assertLess(var, 1.01 * 9.0)
        self.assertEqual(len(V.omega_nodes(0.64)), 13)      # spacing 0.32 deg
        self.assertEqual(len(V.omega_nodes(10.0)), 41)

    def test_s_rel_sigma(self):
        self.assertAlmostEqual(V.s_rel_sigma(dict(K=10.0, K_sig=0.2, P=100.0, P_sig=0.0, e=0.0, e_sig=0.1)),
                               0.02, places=12)
        r = V.s_rel_sigma(dict(K=10.0, K_sig=0.0, P=100.0, P_sig=1.0, e=0.6, e_sig=0.064))
        self.assertAlmostEqual(r, math.hypot(0.01, 0.6 * 0.064 / 0.64), places=12)

    def test_quantile_uniform(self):
        xs, h = V._grid(0.0, 10.0, 10)
        p = [0.1] * 10
        self.assertAlmostEqual(V._quantile(xs, h, p, 0.5), 5.0, places=9)
        self.assertAlmostEqual(V._quantile(xs, h, p, 0.025), 0.25, places=9)

    def test_logsumexp(self):
        self.assertAlmostEqual(V.logsumexp([0.0, 0.0]), math.log(2.0), places=12)
        self.assertEqual(V.logsumexp([-math.inf]), -math.inf)


if __name__ == "__main__":
    unittest.main()
