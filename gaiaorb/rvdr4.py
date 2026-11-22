"""Pre-registered Gaia DR4 measurement of the true masses of RV giant planets
around M dwarfs, with the RV orbital elements fixed.

  v0.1.0 (frozen 2026-11-22): GJ 3512 b and c (PREREGISTRATION.md).
  v0.2.0 (frozen 2026-11-22, Amendment 1, PREREGISTRATION-amendment1.md):
         GJ 317 b and c, GJ 463 b, GJ 9773 b; tests of the published DR3-based
         inclination predictions; mutual inclination of two-planet systems.
         The GJ 3512 analysis and its verdicts are unchanged; v0.2.0 only adds
         reported quantities (folded-inclination posterior, mutual inclination,
         prediction checks).

Part of gj3512-giants-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents
working in MV2's name. Pure standard library.

Physics. The radial velocities (Morales et al. 2025, A&A, arXiv:2507.15516,
Table 4) fix each planet's period P, eccentricity e, argument of periastron
omega, time of conjunction and the star's a sin i = K sqrt(1 - e^2) P / 2 pi.
Gaia sees the star's reflex orbit, of angular semi-major axis
a0 = (a sin i) plx / sin i  >=  a0_min = (a sin i) plx.
With P, e and the periastron time fixed from the RVs, the along-scan signal of
each planet is LINEAR in its four Thiele-Innes constants (A, B, F, G), so the
primary fit is one exact weighted least-squares solve:
  5 single-star parameters + 4 TI constants per planet,
with one jitter term making the reduced chi^2 of that full model equal 1.
The TI constants give a0, i, omega, Omega for each planet; the true mass is
m sin i (RV) / sin i (Gaia), cross-checked by the mass from a0 alone.
"""
import math
import random

from . import kepler, model
from .linalg import wls

VERSION = "0.2.0"
J2017_5_JD = 2457936.875
M_SUN_IN_MJUP = 1047.5655
AU_M = 149597870700.0
DAY_S = 86400.0
DET_DCHI2 = 23.51          # chi^2 with 4 dof, p = 1e-4
MISSING_LAMBDA = 50.0      # expected dchi2 of the minimum-mass orbit above which a non-detection contradicts the RVs
CONSISTENCY_Z = 3.0
BD_MASS_MJUP = 13.0

# Gaia DR3 astrometry (epoch J2016.0), fetched 2026-11-22 from the Gaia TAP.
# GJ 3512: astrometric_excess_noise 0.264 mas, RUWE 1.25, non_single_star 0.
# GJ 317:  excess noise 0.358 mas, RUWE 1.111, non_single_star 0.
# GJ 463:  excess noise 0.183 mas, RUWE 1.407, non_single_star 0.
# GJ 9773: excess noise 0.133 mas, RUWE 1.195, non_single_star 0.
# (The excess-noise values match Morales et al. 2025, Table 1, confirming the
# identifications.) Stellar masses: Morales et al. 2025, Table 1.
STARS = dict(
    gj3512=dict(key="gj3512", name="GJ 3512", dr3_source_id=1040827088726734976,
                dr3_ra=130.331590389547, dr3_dec=59.49165865053633,
                dr3_pmra=-260.27618577042017, dr3_pmdec=-1279.5622335452344,
                dr3_plx=105.29350582131478, dr3_g=13.082282, dr3_ruwe=1.2506707,
                rv_kms=7.132652, mass_msun=0.123, mass_sig=0.009),
    gj317=dict(key="gj317", name="GJ 317", dr3_source_id=5701750715316821248,
               dr3_ra=130.24445565917404, dr3_dec=-23.45269701927076,
               dr3_pmra=-461.344454674079, dr3_pmdec=805.2206235882711,
               dr3_plx=65.87719202272015, dr3_g=10.765108, dr3_ruwe=1.1111134,
               rv_kms=87.10827, mass_msun=0.402, mass_sig=0.014),
    gj463=dict(key="gj463", name="GJ 463", dr3_source_id=1584088431025057792,
               dr3_ra=185.7438361631177, dr3_dec=64.03248298122433,
               dr3_pmra=-673.1157390294371, dr3_pmdec=373.31127424004893,
               dr3_plx=54.44735322749599, dr3_g=10.552364, dr3_ruwe=1.4072167,
               rv_kms=8.916901, mass_msun=0.486, mass_sig=0.014),
    gj9773=dict(key="gj9773", name="GJ 9773", dr3_source_id=2723557788822526720,
                dr3_ra=333.1502482522987, dr3_dec=8.550242428936258,
                dr3_pmra=106.38005183095764, dr3_pmdec=-671.1541035343739,
                dr3_plx=62.87845281368041, dr3_g=10.90822, dr3_ruwe=1.1952326,
                rv_kms=8.1909685, mass_msun=0.358, mass_sig=0.011),
)

# Morales et al. 2025 (A&A; arXiv:2507.15516), Table 4. Tc = time of conjunction
# (BJD); omega = argument of periastron of the star's orbit (RV convention);
# errors are the larger of the two published one-sided errors.
PLANETS = (
    dict(key="b", star="gj3512", P=203.109, P_sig=0.035, K=71.28, K_sig=0.32, e=0.4279, e_sig=0.0036,
         omega=126.80, omega_sig=0.64, Tc=2457739.51, Tc_sig=0.22, msini=0.461, msini_sig=0.023),
    dict(key="c", star="gj3512", P=2354.0, P_sig=36.0, K=28.02, K_sig=0.36, e=0.089, e_sig=0.012,
         omega=50.0, omega_sig=10.0, Tc=2457339.0, Tc_sig=37.0, msini=0.453, msini_sig=0.023),
)
# Morales et al. 2025, Section 6: true masses from the DR3 excess noise alone.
EXCESS_NOISE_MASS = dict(b=(0.78, 0.27, 0.04), c=(0.77, 0.26, 0.04))   # (value, minus, plus) M_Jup

# ---- Amendment 1 (v0.2.0): same paper, Table 4.
SYSTEMS_A1 = dict(
    gj317=(
        dict(key="b", star="gj317", P=695.69, P_sig=0.26, K=71.39, K_sig=0.56, e=0.0826, e_sig=0.0071,
             omega=11.3, omega_sig=4.8, Tc=2458089.4, Tc_sig=1.7, msini=1.688, msini_sig=0.041),
        dict(key="c", star="gj317", P=6500.0, P_sig=169.0, K=29.33, K_sig=0.74, e=0.217, e_sig=0.014,
             omega=106.7, omega_sig=5.4, Tc=2462735.0, Tc_sig=179.0, msini=1.430, msini_sig=0.058),
    ),
    gj463=(
        dict(key="b", star="gj463", P=3459.0, P_sig=22.0, K=36.27, K_sig=0.71, e=0.146, e_sig=0.015,
             omega=329.3, omega_sig=5.9, Tc=2454313.0, Tc_sig=30.0, msini=1.649, msini_sig=0.045),
    ),
    gj9773=(
        dict(key="b", star="gj9773", P=2965.0, P_sig=41.0, K=15.11, K_sig=0.32, e=0.218, e_sig=0.020,
             omega=99.0, omega_sig=5.2, Tc=2459362.0, Tc_sig=18.0, msini=0.543, msini_sig=0.016),
    ),
)
# Table 4 prints e cos(omega) = +0.034 and e sin(omega) = 0.214 for GJ 9773 b,
# which give omega = 81.0 deg, but omega = 99.0 deg. The primary fit uses the
# printed omega (99.0); the same fit with 81.0 is run as well and a verdict for
# GJ 9773 b is given only if both agree (Amendment 1, section 4).
OMEGA_ALTERNATIVE = {("gj9773", "b"): 81.0}

# Published inclination / true-mass predictions made from Gaia DR3 (excess noise,
# and for GJ 463 the Hipparcos-Gaia proper-motion anomaly). inc_fold = the
# inclination folded to [0, 90] deg (prograde and retrograde solutions coincide).
PREDICTIONS = {
    ("gj3512", "b"): (dict(ref="Morales+2025 Table 13", inc_fold=36.0, mass_mjup=0.78),),
    ("gj3512", "c"): (dict(ref="Morales+2025 Table 13", inc_fold=36.0, mass_mjup=0.77),),
    ("gj317", "b"): (dict(ref="Morales+2025 Table 13", inc_fold=33.9, mass_mjup=3.03),),
    ("gj317", "c"): (dict(ref="Morales+2025 Table 13", inc_fold=33.9, mass_mjup=2.56),),
    ("gj463", "b"): (dict(ref="Morales+2025 Table 13", inc_fold=21.8, mass_mjup=4.44),
                     dict(ref="Sozzetti 2023 (A&A 670, L17), prograde", inc_fold=27.0, mass_mjup=3.4),
                     dict(ref="Sozzetti 2023 (A&A 670, L17), retrograde", inc_fold=28.0, mass_mjup=3.6),
                     dict(ref="Piccinini+2026 (arXiv:2601.09401), prograde", inc_fold=27.7, mass_mjup=3.3),
                     dict(ref="Piccinini+2026 (arXiv:2601.09401), retrograde", inc_fold=27.7, mass_mjup=3.4)),
    ("gj9773", "b"): (dict(ref="Morales+2025 Table 13", inc_fold=22.3, mass_mjup=1.43),),
}


# ---------------------------------------------------------------- predictions
def tp_from_tc(Tc, P, e, omega_deg):
    """Periastron time from the time of (inferior) conjunction, at which the
    true anomaly is f = 90 deg - omega (radvel / juliet convention)."""
    f = math.pi / 2.0 - math.radians(omega_deg)
    E = 2.0 * math.atan2(math.sqrt(1.0 - e) * math.sin(f / 2.0), math.sqrt(1.0 + e) * math.cos(f / 2.0))
    M = E - e * math.sin(E)
    return Tc - P * M / (2.0 * math.pi)


def true_anomaly(t, P, e, Tp):
    E = kepler.ecc_anomaly(2.0 * math.pi * (t - Tp) / P, e)
    return 2.0 * math.atan2(math.sqrt(1.0 + e) * math.sin(E / 2.0), math.sqrt(1.0 - e) * math.cos(E / 2.0))


def a_sini_au(K, P, e):
    return K * math.sqrt(1.0 - e * e) * P * DAY_S / (2.0 * math.pi) / AU_M


def a0_min_mas(planet, star=None):
    star = star or STARS[planet["star"]]
    return a_sini_au(planet["K"], planet["P"], planet["e"]) * star["dr3_plx"]


def tp_days(planet):
    """Periastron time in days from J2017.5 (the gaiaorb time origin)."""
    return tp_from_tc(planet["Tc"], planet["P"], planet["e"], planet["omega"]) - J2017_5_JD


def mass_from_a0(a0_mas, planet, star=None, m1=None):
    star = star or STARS[planet["star"]]
    m1 = star["mass_msun"] if m1 is None else m1
    return kepler.companion_mass(a0_mas, star["dr3_plx"], planet["P"], m1) * M_SUN_IN_MJUP


def fold_inc(inc_deg):
    return min(inc_deg, 180.0 - inc_deg)


def circ_dist_deg(a, b):
    d = (a - b) % 360.0
    return min(d, 360.0 - d)


def mutual_inc_deg(i1, W1, i2, W2):
    """Angle between two orbit normals (degrees)."""
    i1, W1, i2, W2 = (math.radians(x) for x in (i1, W1, i2, W2))
    c = math.cos(i1) * math.cos(i2) + math.sin(i1) * math.sin(i2) * math.cos(W1 - W2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


# ---------------------------------------------------------------- the fit
def design(points, planets, mu_r, drop=None):
    cols, y = model.five_cols(points, mu_r)
    for pl in planets:
        if pl["key"] != drop:
            cols = cols + model.orbit_cols(points, pl["P"], pl["e"], tp_days(pl))
    return cols, y


def fit_all(points, planets, mu_r, clip=5.0):
    """Full model (5p + 4 TI per planet) with jitter and one pass of
    |r|/sigma > clip rejection against the FULL model."""
    pts = list(points)
    n_par = 5 + 4 * len(planets)
    n_clipped = 0
    for rnd in range(2):
        cols, y = design(pts, planets, mu_r)
        base = [p["sig"] for p in pts]
        eps = model.fit_jitter(cols, y, base, n_par)
        sig = [math.sqrt(s * s + eps * eps) for s in base]
        f = wls(cols, y, sig)
        if rnd == 1 or clip is None:
            break
        keep = [i for i in range(len(pts)) if abs(f["resid"][i]) <= clip * sig[i]]
        if len(keep) == len(pts):
            break
        n_clipped = len(pts) - len(keep)
        pts = [pts[i] for i in keep]
    return dict(points=pts, sig=sig, eps=eps, fit=f, n_clipped=n_clipped, n=len(pts))


def dchi2_without(pts, sig, planets, mu_r, key, chi2_full):
    cols, y = design(pts, planets, mu_r, drop=key)
    return wls(cols, y, sig)["chi2"] - chi2_full


def expected_lambda(pts, sig, planets, mu_r, planet, inc=90.0, n_node=12):
    """Median over Omega of the noise-free dchi2 that the minimum-mass orbit
    (a0_min, i = 90 deg, RV omega and phase) would give with these epochs."""
    cols_full, _ = design(pts, planets, mu_r)
    cols_drop, _ = design(pts, planets, mu_r, drop=planet["key"])
    lams = []
    for k in range(n_node):
        ti = kepler.campbell_to_ti(a0_min_mas(planet), inc, planet["omega"], 180.0 * k / n_node)
        y = model.orbit_signal(pts, ti, planet["P"], planet["e"], tp_days(planet))
        lams.append(wls(cols_drop, y, sig)["chi2"] - wls(cols_full, y, sig)["chi2"])
    lams.sort()
    return 0.5 * (lams[(n_node - 1) // 2] + lams[n_node // 2])


def cholesky(M):
    n = len(M)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = M[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                L[i][j] = math.sqrt(max(s, 1e-30))
            else:
                L[i][j] = s / L[j][j]
    return L


def pct(xs, q):
    s = sorted(xs)
    k = q * (len(s) - 1)
    lo = int(math.floor(k))
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def summarize(xs):
    return dict(median=pct(xs, 0.5), p16=pct(xs, 0.16), p84=pct(xs, 0.84), p02=pct(xs, 0.025), p98=pct(xs, 0.975))


def planet_posterior(beta, cov, planet, star=None, n_draw=4000, seed=20261202):
    """Monte Carlo over the TI covariance (and m sin i, M_star errors)."""
    star = star or STARS[planet["star"]]
    L = cholesky(cov)
    rng = random.Random(seed)
    a0s, incs, sinis, a0sini, mtrue, ma0, oms, incf = [], [], [], [], [], [], [], []
    for _ in range(n_draw):
        z = [rng.gauss(0.0, 1.0) for _ in range(4)]
        ti = [beta[i] + sum(L[i][k] * z[k] for k in range(i + 1)) for i in range(4)]
        a, inc, om, _Om = kepler.ti_to_campbell(*ti)
        if a <= 0.0 or not math.isfinite(inc):
            continue
        si = math.sin(math.radians(inc))
        a0s.append(a)
        incs.append(inc)
        incf.append(fold_inc(inc))
        sinis.append(si)
        a0sini.append(a * si)
        oms.append(om)
        msini = planet["msini"] + rng.gauss(0.0, planet["msini_sig"])
        mtrue.append(msini / max(si, 1e-6))
        m1 = max(0.05, star["mass_msun"] + rng.gauss(0.0, star["mass_sig"]))
        ma0.append(mass_from_a0(a, planet, star, m1))
    return dict(n=len(a0s), a0_mas=summarize(a0s), inc_deg=summarize(incs), sin_i=summarize(sinis),
                a0_sini_mas=summarize(a0sini), mass_true_mjup=summarize(mtrue),
                mass_from_a0_mjup=summarize(ma0), omega_deg=summarize(oms), inc_fold_deg=summarize(incf))


def mutual_posterior(beta8, cov8, n_draw=4000, seed=20261203):
    """Mutual inclination of two planets from joint draws of their 8 TI
    constants. Astrometry alone cannot tell Omega from Omega + 180 deg, so both
    branches are reported: 'same' (Omegas as returned, in [0, 180)) and
    'flipped' (one Omega + 180 deg)."""
    L = cholesky(cov8)
    rng = random.Random(seed)
    same, flip = [], []
    for _ in range(n_draw):
        z = [rng.gauss(0.0, 1.0) for _ in range(8)]
        ti = [beta8[i] + sum(L[i][k] * z[k] for k in range(i + 1)) for i in range(8)]
        a1, i1, _w1, W1 = kepler.ti_to_campbell(*ti[0:4])
        a2, i2, _w2, W2 = kepler.ti_to_campbell(*ti[4:8])
        if a1 <= 0.0 or a2 <= 0.0 or not (math.isfinite(i1) and math.isfinite(i2)):
            continue
        same.append(mutual_inc_deg(i1, W1, i2, W2))
        flip.append(mutual_inc_deg(i1, W1, i2, W2 + 180.0))
    return dict(n=len(same), same=summarize(same), flipped=summarize(flip))


def run(points, planets=PLANETS, star=None, clip=5.0, log=print):
    star = star or STARS[planets[0]["star"]]
    mu_r = model.mu_r_from_rv(star["rv_kms"], star["dr3_plx"])
    fa = fit_all(points, planets, mu_r, clip)
    pts, sig, f = fa["points"], fa["sig"], fa["fit"]
    res = dict(status="ok", n_transits=fa["n"], n_clipped=fa["n_clipped"], jitter_mas=fa["eps"],
               chi2_full=f["chi2"], median_sigma_mas=sorted(sig)[len(sig) // 2],
               span_years=pts[-1]["t"] - pts[0]["t"], planets={}, rvdr4_version=VERSION)
    for k, pl in enumerate(planets):
        j = 5 + 4 * k
        beta = f["beta"][j:j + 4]
        cov = [row[j:j + 4] for row in f["cov"][j:j + 4]]
        a, inc, om, Om = kepler.ti_to_campbell(*beta)
        d = dchi2_without(pts, sig, planets, mu_r, pl["key"], f["chi2"])
        lam = expected_lambda(pts, sig, planets, mu_r, pl)
        post = planet_posterior(beta, cov, pl, star)
        pred = a0_min_mas(pl, star)
        s = post["a0_sini_mas"]
        width = 0.5 * (s["p84"] - s["p16"])
        z = abs(s["median"] - pred) / width if width > 0 else float("inf")
        om_diff = (om - pl["omega"]) % 180.0
        om_diff = min(om_diff, 180.0 - om_diff)
        # Report only: the Omega branch whose omega lies nearer the RV omega.
        rv_flip = circ_dist_deg(om, pl["omega"]) > 90.0
        res["planets"][pl["key"]] = dict(
            ti=dict(A=beta[0], B=beta[1], F=beta[2], G=beta[3]), ti_cov=cov,
            point=dict(a0_mas=a, inc_deg=inc, omega_deg=om, Omega_deg=Om),
            Omega_rv_branch_deg=(Om + (180.0 if rv_flip else 0.0)) % 360.0, rv_branch_flipped=rv_flip,
            dchi2=d, expected_dchi2_min_mass=lam, a0_min_mas=pred, a0_sini_z=z,
            omega_diff_mod180_deg=om_diff, posterior=post,
            excess_noise_mass_mjup=list(EXCESS_NOISE_MASS.get(pl["key"], ())) if pl["star"] == "gj3512" else [])
        log("%s %s: dchi2=%.1f (expected at min mass %.1f) a0=%.4f mas i=%.1f m_true=%.3f [%.3f, %.3f] MJ z=%.2f" % (
            pl["star"], pl["key"], d, lam, a, inc, post["mass_true_mjup"]["median"],
            post["mass_true_mjup"]["p16"], post["mass_true_mjup"]["p84"], z))
    if len(planets) == 2:
        b8 = f["beta"][5:13]
        c8 = [row[5:13] for row in f["cov"][5:13]]
        p1, p2 = (res["planets"][pl["key"]] for pl in planets)
        mp = mutual_posterior(b8, c8)
        mp["point_same_deg"] = mutual_inc_deg(p1["point"]["inc_deg"], p1["point"]["Omega_deg"],
                                              p2["point"]["inc_deg"], p2["point"]["Omega_deg"])
        mp["point_flipped_deg"] = mutual_inc_deg(p1["point"]["inc_deg"], p1["point"]["Omega_deg"],
                                                 p2["point"]["inc_deg"], p2["point"]["Omega_deg"] + 180.0)
        mp["rv_branch"] = "flipped" if p1["rv_branch_flipped"] != p2["rv_branch_flipped"] else "same"
        res["mutual_inclination"] = mp
    return res


def verdict(pres, validation_ok, status="ok"):
    """Pre-registered verdict for one planet (PREREGISTRATION.md, section 5)."""
    if not validation_ok or status != "ok":
        return "INCONCLUSIVE"
    if pres["dchi2"] < DET_DCHI2:
        return "RV_ORBIT_NOT_SEEN" if pres["expected_dchi2_min_mass"] >= MISSING_LAMBDA else "NOT_DETECTED"
    if pres["a0_sini_z"] > CONSISTENCY_Z:
        return "ORBIT_INCONSISTENT_WITH_RV"
    m = pres["posterior"]["mass_true_mjup"]
    if m["p16"] > BD_MASS_MJUP:
        return "TRUE_MASS_MEASURED_BROWN_DWARF"
    return "TRUE_MASS_MEASURED"


# ---------------------------------------------------------------- Amendment 1
def combine_omega_verdicts(v_primary, v_alt):
    """GJ 9773 b: a verdict only if the printed and the alternative omega agree."""
    return v_primary if v_primary == v_alt else "INCONCLUSIVE_OMEGA_AMBIGUOUS"


def prediction_checks(pres, verdict_str, preds):
    """Amendment 1, section 5: a published DR3-based inclination AGREES if it
    lies inside the DR4 95% interval of the folded inclination, else DISAGREES.
    Evaluated only for TRUE_MASS_MEASURED(_BROWN_DWARF) planets."""
    out = []
    for p in preds:
        row = dict(ref=p["ref"], inc_fold_pred=p["inc_fold"], mass_pred_mjup=p["mass_mjup"])
        if not verdict_str.startswith("TRUE_MASS_MEASURED"):
            row["check"] = "NOT_EVALUATED"
        else:
            s = pres["posterior"]["inc_fold_deg"]
            m = pres["posterior"]["mass_true_mjup"]
            row["check"] = "AGREES" if s["p02"] <= p["inc_fold"] <= s["p98"] else "DISAGREES"
            half = 0.5 * (s["p84"] - s["p16"])
            row["inc_offset_sigma"] = (p["inc_fold"] - s["median"]) / half if half > 0 else None
            row["mass_in_dr4_95"] = bool(m["p02"] <= p["mass_mjup"] <= m["p98"])
        out.append(row)
    return out
