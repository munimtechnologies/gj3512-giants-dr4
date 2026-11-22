"""RV-amplitude-constrained inclination posterior for an RV-fixed Gaia orbit
(Amendment 2, PREREGISTRATION-amendment2.md; fixed 2026-11-22, before Gaia DR4).

Part of gj3512-giants-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents
working in MV2's name. Pure standard library.

Why. rvdr4.py (v0.2.0) reads the inclination from the SHAPE of the fitted
Thiele-Innes ellipse alone. The Epochalypse calibration (2026-11-22) showed that
this passes overall (95% intervals cover the truth in 90.5% of cases) but is
biased toward edge-on for nearly face-on orbits: for true folded inclinations
below 25 deg the intervals covered the truth in 14 of 21 cases, and the median
was too high in 20 of 21. The radial velocities also measure the AMPLITUDE
s = a0 sin i (from K, P, e and the parallax) and the star's argument of
periastron omega. Given s and omega, the astrometric orbit has only the free
parameters i, Omega and a0 (with a0 sin i tied to s within the RV error).

Method. The linear fit gives the TI estimate b and its 4x4 covariance C (the
block of the full fit's covariance, so the other parameters are profiled
exactly); Q = C^-1. For TI = a t(W) with t = u cos W + v sin W (u, v from i
and omega, kepler.campbell_to_ti convention):
  chi2(a, W) = A2 a^2 - 2 A1 a + bQb,  A2 = t'Qt,  A1 = t'Qb.
Prior and RV likelihood: isotropic orientation, log-uniform true mass, and
s_obs ~ N(a sin i, sigma_s). In (i, a) this is, up to a constant,
  sin i * N(a; mu, tau^2),  mu = s_obs / sin i,  tau = sigma_s / sin i,
so the a-integral is exact:
  int N(a; mu, tau^2) exp(-chi2 / 2) da = exp(-E / 2) / sqrt(1 + A2 tau^2),
  E = bQb + (A2 mu^2 - 2 A1 mu - A1^2 tau^2) / (1 + A2 tau^2).
The W-integral uses a periodic trapezoid grid (spectrally accurate); a peak
narrower than the grid is re-integrated on a fine local grid, and if still
narrower, by Laplace's method from the local curvature. The posterior over the
folded inclination i_f in (0, 90] sums i_f and 180 - i_f, and two omega
branches, omega and -omega, with equal weight, because sign conventions for
omega differ between codes (omega + 180 and a mirrored sky are already covered
by the full i and Omega ranges). The RV error of omega is integrated on a
Gaussian-weighted grid of nodes (spacing <= 1.5 deg, +-3 sigma). A coarse grid finds the support; a fine grid over it gives
the quantiles. True mass = m sin i / sin i_f.
"""
import math
import random

from .linalg import cholesky, chol_inverse

VERSION = "0.1.0"
N_W = 180                  # Omega grid (2 deg)
N_REFINE = 400             # points of the local grid over +-1 coarse step
N_COARSE = 180             # folded-inclination grid, coarse pass (0.5 deg)
N_FINE = 240               # fine pass over the support
SUPPORT_TAIL = 1e-6        # mass left out at each end when choosing the fine range
NARROW_FRAC = 0.4          # largest-term share of a sum above which the peak counts as narrow
OMEGA_STEP_MAX = 1.5       # deg, largest spacing of the omega-error nodes

_H = 2.0 * math.pi / N_W
_WTAB = [(math.cos(k * _H), math.sin(k * _H), math.cos(2 * k * _H), math.sin(2 * k * _H)) for k in range(N_W)]


def logsumexp(xs):
    xs = [x for x in xs if x != -math.inf]
    if not xs:
        return -math.inf
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def a_integral(A2, A1, bQb, mu, eps):
    """(log of int N(a; mu, eps) exp(-chi2(a) / 2) da, E). eps = tau^2 >= 0."""
    den = 1.0 + A2 * eps
    E = bQb + (A2 * mu * mu - 2.0 * A1 * mu - A1 * A1 * eps) / den
    return -0.5 * E - 0.5 * math.log(den), E


def _terms(qs, mu, eps, trig):
    uQu, vQv, uQv, uQb, vQb, bQb = qs
    c0, c1 = 0.5 * (uQu + vQv), 0.5 * (uQu - vQv)
    lt, es = [], []
    for cw, sw, c2w, s2w in trig:
        A2 = c0 + c1 * c2w + uQv * s2w
        A1 = cw * uQb + sw * vQb
        den = 1.0 + A2 * eps
        E = bQb + (A2 * mu * mu - 2.0 * A1 * mu - A1 * A1 * eps) / den
        es.append(E)
        lt.append(-0.5 * E - 0.5 * math.log(den))
    return lt, es


def _trig(Ws):
    return [(math.cos(W), math.sin(W), math.cos(2 * W), math.sin(2 * W)) for W in Ws]


def _refine(qs, mu, eps, W0):
    """Integral over [W0 - H, W0 + H] on a fine grid; Laplace if still narrower."""
    step = 2.0 * _H / N_REFINE
    Ws = [W0 - _H + (j + 0.5) * step for j in range(N_REFINE)]
    lt, es = _terms(qs, mu, eps, _trig(Ws))
    m = max(lt)
    tot = sum(math.exp(x - m) for x in lt)
    if 1.0 / tot <= NARROW_FRAC:
        return m + math.log(tot * step), min(es)
    j = lt.index(m)
    if 0 < j < N_REFINE - 1:
        curv = -(lt[j + 1] - 2.0 * lt[j] + lt[j - 1]) / (step * step)
        if curv > 0.0:
            # vertex of the parabola through the three points, then Laplace; cells
            # more than 3 steps from the peak are added (broad components)
            dj = 0.5 * (lt[j - 1] - lt[j + 1]) / (lt[j - 1] - 2.0 * lt[j] + lt[j + 1])
            peak = m + 0.5 * curv * (dj * step) ** 2
            parts = [peak + 0.5 * math.log(2.0 * math.pi / curv)]
            parts += [lt[q] + math.log(step) for q in range(N_REFINE) if abs(q - j) > 3]
            return logsumexp(parts), min(es)
    return m + math.log(tot * step), min(es)


def log_int_w(qs, mu, eps):
    """log of int over W in [0, 2 pi) of int da N(a; mu, eps) exp(-chi2(a, W) / 2),
    and the smallest E met. qs = (uQu, vQv, uQv, uQb, vQb, bQb)."""
    lt, es = _terms(qs, mu, eps, _WTAB)
    m = max(lt)
    e_min = min(es)
    tot = sum(math.exp(x - m) for x in lt)
    if 1.0 / tot <= NARROW_FRAC:
        return m + math.log(tot * _H), e_min
    # Narrow peak(s): each local maximum's window [W_k - H, W_k + H] is integrated
    # finely; it replaces trapezoid cell k and half of cells k - 1 and k + 1.
    acc = []
    n = N_W
    left = [1.0] * n
    for k in range(n):
        x = lt[k]
        if x < m - 30.0 or not (x > lt[k - 1] and x >= lt[(k + 1) % n]):
            continue
        li, e = _refine(qs, mu, eps, k * _H)
        acc.append(li)
        e_min = min(e_min, e)
        left[k] = 0.0
        left[k - 1] = max(0.0, left[k - 1] - 0.5)
        left[(k + 1) % n] = max(0.0, left[(k + 1) % n] - 0.5)
    if not acc:
        return m + math.log(tot * _H), e_min
    acc += [lt[j] + math.log(left[j] * _H) for j in range(n) if left[j] > 0.0]
    return logsumexp(acc), e_min


def uv_unit(inc_deg, omega_deg):
    """TI = a (u cos W + v sin W) (kepler.campbell_to_ti convention)."""
    ci = math.cos(math.radians(inc_deg))
    cw, sw = math.cos(math.radians(omega_deg)), math.sin(math.radians(omega_deg))
    u = (cw, sw * ci, -sw, cw * ci)
    v = (-sw * ci, cw, -cw * ci, -sw)
    return u, v


def _qf(Q, x, y):
    return sum(x[r] * sum(Q[r][c] * y[c] for c in range(4)) for r in range(4))


def qforms(Q, b, bQb, inc_deg, omega_deg):
    u, v = uv_unit(inc_deg, omega_deg)
    return (_qf(Q, u, u), _qf(Q, v, v), _qf(Q, u, v), _qf(Q, u, b), _qf(Q, v, b), bQb)


def chi2_direct(qs, a, W):
    """chi2(a, W) from the quadratic forms (for tests and diagnostics)."""
    uQu, vQv, uQv, uQb, vQb, bQb = qs
    c, s = math.cos(W), math.sin(W)
    A2 = c * c * uQu + 2.0 * c * s * uQv + s * s * vQv
    A1 = c * uQb + s * vQb
    return A2 * a * a - 2.0 * A1 * a + bQb


def omega_nodes(sig_deg):
    """Nodes (offset in deg, weight) for the RV error of omega: a Gaussian-weighted
    grid over +-3 sigma with spacing min(1.5 deg, sigma / 2). Discrete nodes far
    apart would bias the fit when the astrometry pins omega better than the RVs."""
    if not sig_deg > 0.0:
        return ((0.0, 1.0),)
    d = min(OMEGA_STEP_MAX, 0.5 * sig_deg)
    n = int(math.ceil(3.0 * sig_deg / d - 1e-9))
    xs = [j * d for j in range(-n, n + 1)]
    ws = [math.exp(-0.5 * (x / sig_deg) ** 2) for x in xs]
    t = sum(ws)
    return tuple((x, w / t) for x, w in zip(xs, ws))


def log_post_fold(Q, b, s_mas, s_rel_sig, omega_deg, omega_sig_deg, inc_fold_deg, branches=(1.0, -1.0)):
    """Unnormalised log posterior density of the folded inclination, the
    per-branch log integrals (without the sin i prior), and the smallest E met."""
    bQb = _qf(Q, b, b)
    sin_i = math.sin(math.radians(inc_fold_deg))
    if sin_i <= 0.0:
        return -math.inf, [-math.inf] * len(branches), math.inf
    mu = s_mas / sin_i
    eps = (s_rel_sig * s_mas / sin_i) ** 2
    per_branch, e_min = [], math.inf
    for sgn in branches:
        terms = []
        for xo, wo in omega_nodes(omega_sig_deg):
            om = sgn * (omega_deg + xo)
            for inc in (inc_fold_deg, 180.0 - inc_fold_deg):
                li, e = log_int_w(qforms(Q, b, bQb, inc, om), mu, eps)
                e_min = min(e_min, e)
                terms.append(math.log(wo) + li)
        per_branch.append(logsumexp(terms) + math.log(1.0 / len(branches)))
    return logsumexp(per_branch) + math.log(sin_i), per_branch, e_min


def _grid(lo, hi, n):
    h = (hi - lo) / n
    return [lo + (k + 0.5) * h for k in range(n)], h


def _quantile(centers, h, probs, q):
    acc = 0.0
    for x, p in zip(centers, probs):
        if acc + p >= q:
            frac = (q - acc) / p if p > 0 else 0.5
            return x - 0.5 * h + frac * h
        acc += p
    return centers[-1] + 0.5 * h


def summarize(xs):
    s = sorted(xs)

    def pct(q):
        k = q * (len(s) - 1)
        lo = int(math.floor(k))
        hi = min(lo + 1, len(s) - 1)
        return s[lo] + (s[hi] - s[lo]) * (k - lo)
    return dict(median=pct(0.5), p16=pct(0.16), p84=pct(0.84), p02=pct(0.025), p98=pct(0.975))


def posterior(ti, ti_cov, s_mas, s_rel_sig, omega_deg, omega_sig_deg, msini, msini_sig,
              n_draw=4000, seed=20261202, n_coarse=N_COARSE, n_fine=N_FINE):
    """Constrained posterior of the folded inclination and the true mass."""
    try:
        Q = chol_inverse(cholesky([list(r) for r in ti_cov]))
    except ValueError as e:
        return dict(status="bad_covariance", error=str(e))
    b = list(ti)
    args = (Q, b, s_mas, s_rel_sig, omega_deg, omega_sig_deg)
    xc, hc = _grid(0.0, 90.0, n_coarse)
    lc = [log_post_fold(*args, inc_fold_deg=x)[0] for x in xc]
    m = max(lc)
    pc = [math.exp(v - m) for v in lc]
    tot = sum(pc)
    pc = [p / tot for p in pc]
    acc, k_lo = 0.0, 0
    for k, p in enumerate(pc):
        acc += p
        if acc > SUPPORT_TAIL:
            k_lo = k
            break
    acc, k_hi = 0.0, len(pc) - 1
    for k in range(len(pc) - 1, -1, -1):
        acc += pc[k]
        if acc > SUPPORT_TAIL:
            k_hi = k
            break
    lo = max(0.0, xc[k_lo] - 1.5 * hc)
    hi = min(90.0, xc[k_hi] + 1.5 * hc)
    xf, hf = _grid(lo, hi, n_fine)
    lf, branch_acc, e_min = [], [[], []], math.inf
    for x in xf:
        lp, pb, em = log_post_fold(*args, inc_fold_deg=x)
        lf.append(lp)
        e_min = min(e_min, em)
        ls = math.log(math.sin(math.radians(x)))
        for j in range(2):
            branch_acc[j].append(pb[j] + ls)
    m = max(lf)
    pf = [math.exp(v - m) for v in lf]
    tot = sum(pf)
    pf = [p / tot for p in pf]
    inc = dict(median=_quantile(xf, hf, pf, 0.5), p16=_quantile(xf, hf, pf, 0.16),
               p84=_quantile(xf, hf, pf, 0.84), p02=_quantile(xf, hf, pf, 0.025),
               p98=_quantile(xf, hf, pf, 0.975))
    rng = random.Random(seed)
    masses = []
    for _ in range(n_draw):
        i_f = _quantile(xf, hf, pf, rng.random())
        si = math.sin(math.radians(max(i_f, 1e-6)))
        masses.append((msini + rng.gauss(0.0, msini_sig)) / si)
    lb = [logsumexp(v) for v in branch_acc]
    return dict(status="ok", version=VERSION, inc_fold_deg=inc, mass_true_mjup=summarize(masses),
                log_evidence_branch=dict(omega=lb[0], minus_omega=lb[1]),
                preferred_branch="omega" if lb[0] >= lb[1] else "minus_omega",
                chi2_min_constrained=e_min, fine_range_deg=[lo, hi], s_mas=s_mas, s_rel_sig=s_rel_sig,
                omega_deg=omega_deg, omega_sig_deg=omega_sig_deg)


def s_rel_sigma(planet):
    """Relative error of a sin i = K sqrt(1 - e^2) P / 2 pi from the RV errors."""
    rk = planet.get("K_sig", 0.0) / planet["K"]
    rp = planet.get("P_sig", 0.0) / planet["P"]
    e = planet["e"]
    re = e * planet.get("e_sig", 0.0) / (1.0 - e * e)
    return math.sqrt(rk * rk + rp * rp + re * re)


def covered(summary, truth):
    return bool(summary["p02"] <= truth <= summary["p98"])
