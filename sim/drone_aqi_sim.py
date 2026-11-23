#!/usr/bin/env python3
"""
Drone V AQI -- first-order physics simulation.

Three coupled models:
  1. Rotor airflow (actuator-disk momentum theory) -> how much air a hovering
     drone actually moves through its rotors.
  2. Fibrous / electret filter physics (single-fibre filtration theory, Davies
     pressure drop) -> capture efficiency vs particle size, pressure drop,
     and the extra rotor power / flight-time penalty that filter causes.
  3. Well-mixed box models -> what N drones do to PM2.5 in (a) all of Delhi and
     (b) a small "clean-air zone" such as a park or school ground.
  4. Mist / spray scavenging (Slinn below-cloud collection efficiency) ->
     how much water is needed to wash PM10 / PM2.5 out of the air.

Everything is order-of-magnitude engineering physics, not CFD. The point is to
find whether an idea is off by 10x or by 10,000x before building hardware.

Run:  python3 sim/drone_aqi_sim.py
"""
import math
import numpy as np

# ---------------------------------------------------------------- constants
RHO_AIR = 1.25        # kg/m3, cold (~12 C) winter air in Delhi
MU_AIR = 1.78e-5      # Pa s
T_AIR = 285.0         # K
KB = 1.380649e-23
E_CHG = 1.602e-19
EPS0 = 8.854e-12
G = 9.81
RHO_P = 1600.0        # kg/m3, typical urban/smoke PM density
RHO_W = 1000.0
MFP = 0.066e-6        # mean free path of air, m


def cunningham(dp):
    kn = 2 * MFP / dp
    return 1 + kn * (1.257 + 0.4 * np.exp(-1.1 / kn))


def diffusivity(dp):
    return KB * T_AIR * cunningham(dp) / (3 * math.pi * MU_AIR * dp)


# ================================================================ 1. ROTORS
class Drone:
    def __init__(self, name, mass_kg, n_rotors, prop_diam_m, batt_wh,
                 fm=0.6, eta_drive=0.85, filter_mass_kg=0.0):
        self.name = name
        self.mass = mass_kg
        self.n = n_rotors
        self.D = prop_diam_m
        self.batt_wh = batt_wh          # usable energy
        self.fm = fm                    # rotor figure of merit
        self.eta = eta_drive            # motor + ESC efficiency
        self.filter_mass = filter_mass_kg

    @property
    def disk_area(self):
        return self.n * math.pi * (self.D / 2) ** 2

    def hover(self, filt=None, area_ratio=1.0):
        """Return dict of hover state, optionally with a filter in the flow.

        A filter anywhere in the rotor stream (above = inflow, below = wake)
        adds a pressure loss dp. Momentum balance: the rotor must make up that
        dp over the full flow Q, so extra shaft power ~ Q*dp. Extra mass of
        the filter also raises thrust and induced velocity. We iterate.
        area_ratio = filter media area / rotor disk area (pleating > 1).
        """
        m = self.mass + (self.filter_mass if filt else 0.0)
        T = m * G
        A = self.disk_area
        vi = math.sqrt(T / (2 * RHO_AIR * A))       # induced velocity at disk
        Q = vi * A                                  # m3/s through rotors
        p_ideal = T * vi
        p_rotor = p_ideal / self.fm
        dp = 0.0
        eff = {}
        if filt:
            U = Q / (A * area_ratio)                # face velocity on media
            dp = filt.pressure_drop(U)
            # loss power is supplied by an off-design rotor: assume it is
            # added at the same aerodynamic efficiency as the rotor (FM)
            p_rotor += Q * dp / self.fm
            eff = filt.efficiency_table(U)
            eff["face_velocity"] = U
        p_elec = p_rotor / self.eta
        endurance_min = self.batt_wh / p_elec * 60
        return dict(mass=m, thrust_N=T, disk_loading_Pa=T / A, vi=vi, Q=Q,
                    filter_dp_Pa=dp, power_W=p_elec,
                    endurance_min=endurance_min, eff=eff)


# ================================================================ 2. FILTERS
class FibrousFilter:
    """Single-fibre filtration theory (Hinds, Aerosol Technology, ch.9).

    d_f: fibre diameter, alpha: solidity, L: thickness,
    sigma: electret surface charge on fibres (C/m2); 0 = uncharged media.
    """

    def __init__(self, name, d_f, alpha, L, sigma=0.0, eps_p=2.5):
        self.name, self.d_f, self.alpha, self.L = name, d_f, alpha, L
        self.sigma, self.eps_p = sigma, eps_p

    def pressure_drop(self, U):
        a = self.alpha
        return MU_AIR * U * self.L * 64 * a ** 1.5 * (1 + 56 * a ** 3) / self.d_f ** 2

    def single_fibre_eff(self, dp, U):
        a, df = self.alpha, self.d_f
        Ku = -math.log(a) / 2 - 0.75 + a - a * a / 4
        cc = cunningham(dp)
        # diffusion
        Pe = U * df / diffusivity(dp)
        eD = 2 * Pe ** (-2 / 3)
        # interception
        R = dp / df
        eR = (1 - a) * R ** 2 / (Ku * (1 + R))
        # inertial impaction (Stokes number based on fibre)
        Stk = RHO_P * dp ** 2 * cc * U / (18 * MU_AIR * df)
        J = (29.6 - 28 * a ** 0.62) * R ** 2 - 27.5 * R ** 2.8 if R < 0.4 else 2.0
        eI = Stk * max(J, 0) / (2 * Ku ** 2)
        # electrostatic (electret fibres)
        eE = 0.0
        if self.sigma > 0:
            qf = math.pi * df * self.sigma               # C per metre of fibre
            # polarisation force on neutral particles
            kin = ((self.eps_p - 1) / (self.eps_p + 2)) * cc * dp ** 2 * qf ** 2 / (
                3 * math.pi ** 2 * EPS0 * MU_AIR * df ** 3 * U)
            e_in = 0.06 * kin ** 0.4
            # Coulomb force on Boltzmann-charged particles (mean |n|)
            n_avg = 2.37 * math.sqrt(dp * 1e6)
            kc = cc * n_avg * E_CHG * qf / (3 * math.pi ** 2 * EPS0 * MU_AIR * dp * df * U)
            e_c = 0.59 * kc ** 0.83
            eE = 1 - (1 - min(e_in, 1)) * (1 - min(e_c, 1))
        tot = 1 - np.prod([1 - min(e, 1.0) for e in (eD, eR, eI, eE)])
        return tot

    def penetration(self, dp, U):
        e = self.single_fibre_eff(dp, U)
        a = self.alpha
        return math.exp(-4 * a * e * self.L / (math.pi * self.d_f * (1 - a)))

    def efficiency(self, dp, U):
        return 1 - self.penetration(dp, U)

    def mass_eff(self, U, mmd, gsd, dmax):
        """Mass-weighted efficiency over a lognormal mass distribution,
        truncated at dmax (2.5 um for PM2.5, 10 um for PM10)."""
        d = np.logspace(math.log10(0.02e-6), math.log10(dmax), 300)
        w = np.exp(-(np.log(d / mmd)) ** 2 / (2 * math.log(gsd) ** 2))
        e = np.array([self.efficiency(x, U) for x in d])
        return float((w * e).sum() / w.sum())

    def efficiency_table(self, U):
        out = {f"{s}um": self.efficiency(s * 1e-6, U) for s in (0.1, 0.3, 1.0, 2.5, 10.0)}
        # Delhi winter smoke-dominated PM2.5: accumulation mode MMD ~0.5 um
        out["PM2.5_mass"] = self.mass_eff(U, 0.5e-6, 2.0, 2.5e-6)
        # coarse dust mode for PM10
        out["PM10_mass"] = self.mass_eff(U, 4.0e-6, 2.2, 10e-6)
        return out


FILTERS = {
    "electret_coarse": FibrousFilter("Electret coarse (MERV-8-like, charged)", 20e-6, 0.03, 3e-3, sigma=1e-5),
    "uncharged_coarse": FibrousFilter("Same media, discharged", 20e-6, 0.03, 3e-3, sigma=0.0),
    "electret_fine": FibrousFilter("Electret fine (N95-like media)", 8e-6, 0.06, 1.5e-3, sigma=1e-5),
    "hepa_glass": FibrousFilter("HEPA glass-fibre media", 1.0e-6, 0.05, 0.4e-3, sigma=0.0),
}


# ================================================================ 3. BOX MODELS
def box_reduction(volume_m3, cadr_total_m3s, k_vent):
    """Steady-state fractional PM reduction in a ventilated well-mixed box.
    dC/dt = k_vent*(C_bg - C) - (CADR/V)*C  ->  C/C_bg = k_vent/(k_vent+k_d)."""
    k_d = cadr_total_m3s / volume_m3
    return k_d / (k_vent + k_d), 1 / k_d


# ================================================================ 4. MIST
def drop_terminal_velocity(D):
    """Terminal velocity of a water drop (m/s): Stokes below 80 um,
    Rogers & Yau linear law u = 8000 s^-1 * r up to 1.2 mm,
    Atlas et al. (1973) above."""
    if D < 80e-6:
        return (RHO_W - RHO_AIR) * G * D ** 2 / (18 * MU_AIR)
    if D < 1.2e-3:
        return 8000.0 * D / 2
    return 9.65 - 10.3 * math.exp(-0.6 * D * 1e3)


def slinn_efficiency(dp, D):
    """Slinn (1983) collection efficiency of a falling drop D for particle dp."""
    Ut = drop_terminal_velocity(D)
    Re = (D / 2) * Ut * RHO_AIR / MU_AIR
    Sc = MU_AIR / (RHO_AIR * diffusivity(dp))
    phi = dp / D
    omega = 1.0e-3 / MU_AIR
    tau = RHO_P * dp ** 2 * cunningham(dp) / (18 * MU_AIR)
    St = 2 * tau * Ut / D
    S_star = (1.2 + math.log(1 + Re) / 12) / (1 + math.log(1 + Re))
    e_b = 4 / (Re * Sc) * (1 + 0.4 * Re ** 0.5 * Sc ** (1 / 3) + 0.16 * Re ** 0.5 * Sc ** 0.5)
    e_int = 4 * phi * (1 / omega + (1 + 2 * Re ** 0.5) * phi)
    e_imp = 0.0
    if St > S_star:
        e_imp = ((St - S_star) / (St - S_star + 2 / 3)) ** 1.5 * (RHO_W / RHO_P) ** 0.5
    return min(e_b + e_int + e_imp, 1.0)


def mass_scav_eff(D, mmd, gsd, dmax):
    d = np.logspace(math.log10(0.02e-6), math.log10(dmax), 200)
    w = np.exp(-(np.log(d / mmd)) ** 2 / (2 * math.log(gsd) ** 2))
    e = np.array([slinn_efficiency(x, D) for x in d])
    return float((w * e).sum() / w.sum())


def evap_time(D, rh, T_c=15.0):
    """Lifetime of a pure water drop in sub-saturated air (diffusion-limited)."""
    Dv = 2.4e-5
    psat = 610.94 * math.exp(17.625 * T_c / (T_c + 243.04))
    rho_vs = psat * 0.018 / (8.314 * (T_c + 273.15))
    return RHO_W * D ** 2 / (8 * Dv * rho_vs * max(1 - rh, 1e-6))


# ================================================================ REPORT
def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def main():
    drones = [
        Drone("Small quad 2 kg, 4x10in props", 2.0, 4, 0.254, batt_wh=95, filter_mass_kg=0.15),
        Drone("Heavy octo 25 kg, 8x30in props (agri-spray class)", 25.0, 8, 0.762,
              batt_wh=1400, filter_mass_kg=1.5),
    ]

    hr("1. ROTOR AIRFLOW (momentum theory, hover, no filter)")
    for d in drones:
        h = d.hover()
        print(f"{d.name}\n  disk area {d.disk_area:.2f} m2 | disk loading {h['disk_loading_Pa']:.0f} Pa"
              f" | induced vel {h['vi']:.1f} m/s | flow Q {h['Q']:.2f} m3/s"
              f" | power {h['power_W']:.0f} W | endurance {h['endurance_min']:.0f} min")

    hr("2. FILTERS IN THE ROTOR STREAM")
    results = {}
    for d in drones:
        print(f"\n--- {d.name}")
        base = d.hover()
        for key, f in FILTERS.items():
            for ar in (1.0, 5.0):
                h = d.hover(f, area_ratio=ar)
                e = h["eff"]
                cadr = h["Q"] * e["PM2.5_mass"]
                results[(d.name, key, ar)] = (h, cadr)
                print(f"  {f.name:<40s} media/disk={ar:>3.0f}x  U={e['face_velocity']:.2f} m/s"
                      f"  dp={h['filter_dp_Pa']:7.0f} Pa  power={h['power_W']:7.0f} W"
                      f" ({h['power_W'] / base['power_W']:.1f}x)  endur={h['endurance_min']:5.1f} min")
                print(f"      eff 0.1um {e['0.1um']:.2f} | 0.3um {e['0.3um']:.2f} | 1um {e['1.0um']:.2f}"
                      f" | 2.5um {e['2.5um']:.2f} | 10um {e['10.0um']:.2f}"
                      f" | PM2.5-mass {e['PM2.5_mass']:.2f} | PM10-mass {e['PM10_mass']:.2f}"
                      f" | CADR(PM2.5) {cadr:.2f} m3/s")

    # pick the realistic best case: heavy octo, electret coarse, pleated 5x
    h_s, cadr_s = results[(drones[0].name, "electret_coarse", 5.0)]
    h_l, cadr_l = results[(drones[1].name, "electret_coarse", 5.0)]

    hr("3. ENERGY COST PER m3 OF CLEAN AIR (PM2.5)")
    print(f"  small quad + electret     : {h_s['power_W'] / cadr_s:7.0f} J/m3 (CADR {cadr_s:.2f} m3/s)")
    print(f"  heavy octo + electret     : {h_l['power_W'] / cadr_l:7.0f} J/m3 (CADR {cadr_l:.2f} m3/s)")
    print(f"  ground HEPA room purifier : {50 / 0.14:7.0f} J/m3 (50 W, CADR ~300 cfm = 0.14 m3/s)")
    print(f"  large ground fan + electret wall (2 kW, 30 m3/s, 60% eff): {2000 / 18:7.0f} J/m3")

    hr("4. CITY SCALE: 100 drones over Delhi NCT")
    area = 1.483e9           # m2
    L_city = math.sqrt(area)
    duty = 0.6               # fraction of time airborne (battery swaps)
    for N, cadr, label in ((100, cadr_l, "100 heavy octos"), (100, cadr_s, "100 small quads"),
                           (10000, cadr_l, "10,000 heavy octos")):
        for H in (150, 400):
            for wind in (0.5, 2.0):
                k_vent = wind / L_city
                red, tau = box_reduction(area * H, N * cadr * duty, k_vent)
                print(f"  {label:<20s} mixing ht {H:>3d} m, wind {wind:.1f} m/s:"
                      f"  PM2.5 reduction {red * 100:8.4f} %   (drone-only clean-out time {tau / 3.15e7:,.1f} yr)")

    hr("5. LOCAL CLEAN-AIR ZONE (park / school ground 200 m x 200 m)")
    for H in (20, 50):
        for wind in (0.3, 1.0):
            # exchange with surroundings: horizontal flushing + vertical entrainment
            k_vent = wind / 200 + 0.002
            for N in (10, 100):
                red, _ = box_reduction(200 * 200 * H, N * cadr_l * duty, k_vent)
                print(f"  {N:>3d} heavy octos, layer {H} m, wind {wind} m/s: PM2.5 reduction {red * 100:5.1f} %")

    hr("6. MIST / SPRAY SCAVENGING")
    print("  Collection efficiency E of one falling drop (Slinn 1983):")
    print("    drop D      0.1um    0.3um    1um      2.5um    10um   | PM2.5-mass PM10-mass  fall m/s  life@RH70%  life@RH95%")
    for D in (20e-6, 50e-6, 100e-6, 200e-6, 500e-6, 1e-3, 2e-3):
        es = [slinn_efficiency(x * 1e-6, D) for x in (0.1, 0.3, 1.0, 2.5, 10.0)]
        m25 = mass_scav_eff(D, 0.5e-6, 2.0, 2.5e-6)
        m10 = mass_scav_eff(D, 4.0e-6, 2.2, 10e-6)
        print(f"    {D * 1e6:6.0f} um  " + "  ".join(f"{x:7.4f}" for x in es)
              + f"  |  {m25:7.4f}   {m10:7.4f}   {drop_terminal_velocity(D):6.2f}"
              f"   {evap_time(D, 0.70):8.1f} s   {evap_time(D, 0.95):8.1f} s")

    print("\n  Water needed to remove 50% of the particles in the column the drops fall through")
    print("  (best case: no evaporation, drops traverse whole mixing layer).")
    print("  Removed fraction = 1 - exp(-1.5 * E * W / D),  W = water depth sprayed (m)")
    for D in (50e-6, 100e-6, 200e-6, 1e-3):
        for lab, E in (("PM2.5", mass_scav_eff(D, 0.5e-6, 2.0, 2.5e-6)),
                       ("PM10", mass_scav_eff(D, 4.0e-6, 2.2, 10e-6)),
                       ("PM2.5 charged mist (assume E=0.1)", 0.1)):
            if lab.startswith("PM2.5 charged") and D != 100e-6:
                continue
            W = math.log(2) * D / (1.5 * max(E, 1e-9))
            litres_delhi = W * area * 1e3
            sorties = litres_delhi / 40.0
            print(f"    D={D * 1e6:5.0f} um {lab:<34s} E={E:.4f}  water {W * 1e3:9.3f} mm"
                  f"  -> Delhi {litres_delhi / 1e9:9.2f} billion L  = {sorties:,.0f} sorties of a 40 L drone")
    print("\n  For scale: one 40 L agri-drone sortie over 1 ha = 0.004 mm of water.")


if __name__ == "__main__":
    main()
