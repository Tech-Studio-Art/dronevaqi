#!/usr/bin/env python3
"""
0-D / box models for the parts of Solutions 3 and 4 where CFD adds nothing:

 3b. Vehicle ram-air filters  -- Bernoulli ram pressure vs Davies filter
     pressure drop (Hinds 1999) -> self-driven flow, CADR per bus, road box.
 4b. Bipolar-ion agglomeration -- sectional Smoluchowski coagulation
     (Seinfeld & Pandis 2016, ch.13) with the Fuchs transition-regime kernel
     and the Fuchs (1964) electrostatic correction W^-1 = y/(exp(y)-1);
     particles diffusion-charged per White (1951).
 4c. Charged mist             -- Slinn (1983) collection efficiency plus the
     Coulomb term of Kraemer & Johnstone (1955), E_c = 4|K_c|, drops charged
     to 15 % of the Rayleigh limit (~2 mC/kg for 100 um, typical of
     induction-charged agricultural sprayers); evaporation per Maxwell; junction box.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from drone_aqi_sim import (FILTERS, slinn_efficiency, drop_terminal_velocity, evap_time,  # noqa: E402
                           cunningham, diffusivity, KB, T_AIR, MU_AIR, RHO_P, E_CHG, EPS0)

RHO = 1.25
C25, C10 = 300.0, 500.0
OUT = {}


# ================================================================ 3b vehicles
def vehicles():
    f = FILTERS["electret_coarse"]
    res = []
    for kmh in (15, 25, 40):
        V = kmh / 3.6
        p_ram = 0.8 * 0.5 * RHO * V ** 2          # Cp ~0.8 at a front-face intake
        A_in, AR = 0.4, 6.0                        # intake 0.4 m2, pleat ratio 6
        K_duct = 1.5                               # duct + louvre loss coeff on intake vel
        # solve p_ram = dp_filter(U/AR) + K_duct 0.5 rho U^2 for intake velocity U
        lo, hi = 0.0, V
        for _ in range(60):
            U = 0.5 * (lo + hi)
            r = f.pressure_drop(U / AR) + K_duct * 0.5 * RHO * U ** 2 - p_ram
            lo, hi = (U, hi) if r < 0 else (lo, U)
        Q = U * A_in
        e25 = f.mass_eff(U / AR, 0.5e-6, 2.0, 2.5e-6)
        e10 = f.mass_eff(U / AR, 4.0e-6, 2.2, 10e-6)
        res.append(dict(speed_kmh=kmh, ram_Pa=p_ram, intake_U=U, Q=Q, eff_PM25=e25, eff_PM10coarse=e10,
                        CADR_PM25=Q * e25, CADR_PM10=Q * (C25 * e25 + (C10 - C25) * e10) / C10))
    # road-corridor box: 1 km of arterial road, 30 m wide, 15 m mixing, cross-wind 1 m/s
    V_box = 1000 * 30 * 15
    k_vent = 1.0 / 30 + 0.005
    corridor = []
    for buses_per_km in (5, 20, 60):
        cadr = buses_per_km * res[1]["CADR_PM25"]
        k_d = cadr / V_box
        corridor.append(dict(buses_on_1km=buses_per_km, PM25_reduction_pct=100 * k_d / (k_d + k_vent)))
    return dict(per_bus=res, corridor=corridor,
                fleet_note="Delhi city-bus fleet ~7,000-8,000 buses (order of magnitude).")


# ================================================================ 4b coagulation
def fuchs_kernel(d1, d2):
    """Fuchs transition-regime Brownian coagulation kernel [m3/s] (Seinfeld & Pandis Table 13.1)."""
    def props(d):
        D = diffusivity(d)
        m = RHO_P * math.pi / 6 * d ** 3
        c = math.sqrt(8 * KB * T_AIR / (math.pi * m))
        lam = 8 * D / (math.pi * c)
        g = (math.sqrt(2) / (3 * d * lam)) * ((d + lam) ** 3 - (d ** 2 + lam ** 2) ** 1.5) - d
        return D, c, g
    D1, c1, g1 = props(d1)
    D2, c2, g2 = props(d2)
    ds = d1 + d2
    Ds = D1 + D2
    cs = math.sqrt(c1 ** 2 + c2 ** 2)
    gs = math.sqrt(g1 ** 2 + g2 ** 2)
    return 2 * math.pi * ds * Ds / (ds / (ds + 2 * gs) + 8 * Ds / (cs * ds))


def white_charge(d, Nit):
    """Diffusion charging (White 1951): mean charges after ion exposure Ni*t [s/m3]."""
    ci = 240.0  # m/s mean ion thermal speed
    KE = 1 / (4 * math.pi * EPS0)
    return d * KB * T_AIR / (2 * KE * E_CHG ** 2) * math.log(
        1 + math.pi * KE * d * ci * E_CHG ** 2 * Nit / (2 * KB * T_AIR))


def coagulation(hours=(1, 6), Nit_values=(0, 1e12, 1e13)):
    # Delhi winter fine mode: MMD 0.5 um, GSD 2, 300 ug/m3  -> sections
    edges = np.logspace(np.log10(0.02e-6), np.log10(20e-6), 61)
    dm = np.sqrt(edges[:-1] * edges[1:])
    vol = math.pi / 6 * dm ** 3
    w = np.exp(-(np.log(dm / 0.5e-6)) ** 2 / (2 * math.log(2.0) ** 2))
    mass = w / w.sum() * C25 * 1e-9          # kg/m3 per section
    N0 = mass / (RHO_P * vol)                # #/m3
    out = []
    for Nit in Nit_values:
        q = np.array([white_charge(d, Nit) if Nit > 0 else 0.0 for d in dm])
        K = np.zeros((len(dm), len(dm)))
        for i in range(len(dm)):
            for j in range(len(dm)):
                b = fuchs_kernel(dm[i], dm[j])
                if Nit > 0:
                    # half the population +, half -: average of attracting and repelling pairs
                    y = q[i] * q[j] * E_CHG ** 2 / (4 * math.pi * EPS0 * 0.5 * (dm[i] + dm[j]) * KB * T_AIR)
                    y = min(y, 50)
                    f_rep = y / math.expm1(y) if y > 1e-9 else 1.0
                    f_att = (-y) / math.expm1(-y) if y > 1e-9 else 1.0
                    b *= 0.5 * (f_rep + f_att)
                K[i, j] = b
        N = N0.copy()
        dt, t, rec = 10.0, 0.0, {}
        targets = sorted(h * 3600 for h in hours)
        # volume-conserving sectional scheme (Jacobson 2005 semi-implicit)
        f_split = np.zeros((len(dm), len(dm), len(dm)))
        for i in range(len(dm)):
            for j in range(len(dm)):
                v = vol[i] + vol[j]
                k = np.searchsorted(vol, v) - 1
                k = min(max(k, 0), len(dm) - 1)
                if k < len(dm) - 1:
                    fk = (vol[k + 1] - v) / (vol[k + 1] - vol[k]) * vol[k] / v
                    f_split[i, j, k] = fk
                    f_split[i, j, k + 1] = 1 - fk
                else:
                    f_split[i, j, k] = 1.0
        vv = vol[:, None] + vol[None, :]
        while targets:
            P = 0.5 * np.einsum("ijk,ij->k", f_split, K * vv * np.outer(N, N))
            L = K @ N
            N = (vol * N + dt * P) / (vol * (1 + dt * L))
            t += dt
            if t >= targets[0]:
                m = N * vol * RHO_P
                cum = np.cumsum(m) / m.sum()
                rec[f"{int(targets[0] / 3600)}h"] = dict(
                    mass_conservation=float(m.sum() / (N0 * vol * RHO_P).sum()),
                    number_change_pct=100 * (N.sum() / N0.sum() - 1),
                    MMD_um=float(np.interp(0.5, cum, dm) * 1e6),
                    mass_above_2p5um_pct=float(100 * m[dm > 2.5e-6].sum() / m.sum()))
                targets.pop(0)
            if t > 600:
                dt = 60.0
        out.append(dict(Ni_t=Nit, mean_charge_0p5um=white_charge(0.5e-6, Nit) if Nit else 0.0, **rec))
    m0 = N0 * vol * RHO_P
    return dict(initial_MMD_um=float(np.interp(0.5, np.cumsum(m0) / m0.sum(), dm) * 1e6), initial_number_per_cm3=float(N0.sum() * 1e-6), runs=out)


# ================================================================ 4c charged mist
def rayleigh_q(D, frac=0.15):
    return frac * 8 * math.pi * math.sqrt(EPS0 * 0.072 * (D / 2) ** 3)


def coulomb_eff(dp, D, n_p):
    U = drop_terminal_velocity(D)
    Q = rayleigh_q(D)
    Kc = cunningham(dp) * n_p * E_CHG * Q / (3 * math.pi ** 2 * EPS0 * MU_AIR * dp * D ** 2 * U)
    return min(4 * Kc, 1.0)


def mass_avg(fn, mmd, gsd, dlo, dhi):
    d = np.logspace(np.log10(dlo), np.log10(dhi), 120)
    w = np.exp(-(np.log(d / mmd)) ** 2 / (2 * math.log(gsd) ** 2))
    return float((w * np.array([fn(x) for x in d])).sum() / w.sum())


def charged_mist():
    rows = []
    for D in (50e-6, 100e-6, 200e-6):
        for mode, npf in (("uncharged", None), ("charged drop, ambient Boltzmann particles", "boltz"),
                          ("charged drop + pre-charged particles (Ni t=1e13)", "pre")):
            def E(dp):
                e = slinn_efficiency(dp, D)
                if npf == "boltz":
                    n = 2.37 * math.sqrt(dp * 1e6)
                    e = e + 0.5 * coulomb_eff(dp, D, n)      # only opposite-sign half attracted
                elif npf == "pre":
                    e = e + coulomb_eff(dp, D, max(white_charge(dp, 1e13), 1))
                return min(e, 1.0)
            e25 = mass_avg(E, 0.5e-6, 2.0, 0.02e-6, 2.5e-6)
            ec = mass_avg(E, 4e-6, 2.2, 2.5e-6, 10e-6)
            rows.append(dict(D_um=D * 1e6, mode=mode,
                             q_over_m_mC_per_kg=rayleigh_q(D) / (1000 * math.pi / 6 * D ** 3) * 1e3 if npf else 0, E_PM25=e25, E_coarse=ec,
                             fall_m_s=drop_terminal_velocity(D), life_s_RH85=evap_time(D, 0.85),
                             fall_before_evap_m=drop_terminal_velocity(D) * evap_time(D, 0.85) * 0.5))
    # junction box: 50 m x 50 m, 10 m layer, 1 m/s wind -> air replaced every 50 s
    box = []
    A, H, u = 2500.0, 10.0, 1.0
    k_vent = u / 50.0
    for r in rows:
        if r["D_um"] != 100:
            continue
        # scavenging coefficient Lambda = 1.5 E Wdot / D, Wdot = water depth rate (m/s);
        # only the part of the layer drops traverse before evaporating counts
        f_layer = min(r["fall_before_evap_m"] / H, 1.0)
        for lph in (1000, 5000, 20000):
            Wdot = lph / 1000 / 3600 / A
            lam25 = 1.5 * r["E_PM25"] * Wdot / 100e-6 * f_layer
            lamc = 1.5 * r["E_coarse"] * Wdot / 100e-6 * f_layer
            red25 = lam25 / (lam25 + k_vent)
            redc = lamc / (lamc + k_vent)
            red10 = (C25 * red25 + (C10 - C25) * redc) / C10
            box.append(dict(mode=r["mode"], water_L_per_h=lph, PM25_red_pct=100 * red25,
                            PM10_red_pct=100 * red10))
    return dict(drops=rows, junction_box=box)


if __name__ == "__main__":
    OUT["vehicles"] = vehicles()
    OUT["charged_mist"] = charged_mist()
    OUT["coagulation"] = coagulation()
    os.makedirs(os.path.join(HERE, "runs"), exist_ok=True)
    json.dump(OUT, open(os.path.join(HERE, "runs", "models0d.json"), "w"), indent=1, default=float)
    print(json.dumps(OUT, indent=1, default=float))
