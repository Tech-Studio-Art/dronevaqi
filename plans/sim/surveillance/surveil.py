#!/usr/bin/env python3
"""
24x7 pollution-source surveillance: detection, localisation and avoided emissions.

Methods (standard, from the literature):
  * Dispersion: steady Gaussian plume, ground-level source, ground reflection,
    Briggs urban dispersion coefficients (Hanna, Briggs & Hosker 1982):
      stability D (day/overcast):  sy = 0.16x(1+0.0004x)^-1/2,  sz = 0.14x(1+0.0003x)^-1/2
      stability E-F (night):       sy = 0.11x(1+0.0004x)^-1/2,  sz = 0.08x(1+0.0015x)^-1/2
  * Source: open burning of domestic waste, PM2.5 emission factor ~10 g/kg
    (Wiedinmyer et al. 2014), 100 kg/h burn rate -> Q = 0.28 g/s.
  * Network: square grid of low-cost PM sensors (Snyder et al. 2013) at spacing s.
    A sensor "detects" when the 15-min plume increment exceeds 25 ug/m3
    (well above low-cost sensor noise and short-term background swings).
  * Monte Carlo: random source position in a grid cell, random wind direction.
  * Localisation: the detecting sensors and the wind direction give an upwind
    search sector; a patrol drone (10 m/s, thermal camera) flies the sector.
  * Avoided emission = Q x (burn duration - response time).
"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(42)

EF = 10.0            # g PM2.5 / kg waste
BURN = 100.0         # kg/h
Q = EF * BURN / 3600  # g/s
THRESH = 25e-6       # g/m3
Z_SENSOR = 3.0       # m (pole / rooftop mount)
CLASSES = {"night (stable, E-F)": dict(u=1.0, a=(0.11, 0.0004, 0.08, 0.0015)),
           "day (neutral, D)": dict(u=2.5, a=(0.16, 0.0004, 0.14, 0.0003))}


def sig(x, a):
    sy = a[0] * x / np.sqrt(1 + a[1] * x)
    sz = a[2] * x / np.sqrt(1 + a[3] * x)
    return sy, sz


def conc(xd, yc, z, u, a, Hs=1.0):
    """Ground-reflected Gaussian plume, g/m3. xd downwind, yc crosswind (m)."""
    out = np.zeros_like(xd, dtype=float)
    m = xd > 1.0
    sy, sz = sig(xd[m], a)
    out[m] = (Q / (2 * np.pi * u * sy * sz) * np.exp(-yc[m] ** 2 / (2 * sy ** 2))
              * (np.exp(-(z - Hs) ** 2 / (2 * sz ** 2)) + np.exp(-(z + Hs) ** 2 / (2 * sz ** 2))))
    return out


def detect_prob(s, cls, n=4000):
    a, u = CLASSES[cls]["a"], CLASSES[cls]["u"]
    k = np.arange(-6, 7) * s
    gx, gy = np.meshgrid(k, k)
    gx, gy = gx.ravel(), gy.ravel()
    hits, nsens, dist = 0, [], []
    for _ in range(n):
        sx, sy_ = rng.uniform(0, s, 2)
        th = rng.uniform(0, 2 * np.pi)
        dx, dy = gx - sx, gy - sy_
        xd = dx * np.cos(th) + dy * np.sin(th)
        yc = -dx * np.sin(th) + dy * np.cos(th)
        c = conc(xd, yc, Z_SENSOR, u, a)
        det = c > THRESH
        if det.any():
            hits += 1
            nsens.append(int(det.sum()))
            dist.append(float(np.hypot(dx[det], dy[det]).min()))
    return dict(spacing_m=s, p_detect=hits / n, mean_sensors=float(np.mean(nsens)) if nsens else 0.0,
                mean_dist_m=float(np.mean(dist)) if dist else float("nan"))


def response_time_h(dist_m, base_detect_min=15, drone_speed=10.0, search_min=5, ground_team_min=30):
    """detection averaging + drone transit + sector search + ground team to suppress."""
    return (base_detect_min + dist_m / drone_speed / 60 + search_min + ground_team_min) / 60


def patrol(n_drones, area_km2=25.0, swath_km=0.3, speed_kmh=36.0, duty=0.7, T_burn=4.0,
           p_sensor=0.32, t_sensor_h=0.35, ground_h=0.5, n=20000):
    """Thermal-camera patrol: lawnmower route, revisit interval R = area / (drones x swath x speed x duty).
    A fire starts at a random time in the cycle -> detect delay ~ U(0, R). If a fixed sensor detects
    first (probability p_sensor, delay t_sensor_h incl. drone dispatch), the earlier time wins."""
    R = area_km2 / (n_drones * swath_km * speed_kmh * duty) if n_drones else float("inf")
    t_drone = rng.uniform(0, R, n) if n_drones else np.full(n, np.inf)
    t_sens = np.where(rng.uniform(0, 1, n) < p_sensor, t_sensor_h, np.inf)
    t_det = np.minimum(t_drone, t_sens)
    t_stop = np.minimum(t_det + ground_h, T_burn)
    burned = t_stop / T_burn
    return dict(drones=n_drones, revisit_h=R, p_detect_in_burn=float((t_det < T_burn).mean()),
                mean_detect_h=float(np.mean(np.minimum(t_det, T_burn))),
                emission_avoided_pct=float(100 * (1 - burned.mean())))


def main():
    spacings = [250, 500, 1000, 2000, 6000]
    det = {cls: [detect_prob(s, cls) for s in spacings] for cls in CLASSES}
    # avoided emissions per fire vs response time (burn duration 4 h if unattended)
    T_BURN = 4.0
    resp = np.linspace(0.25, 4, 16)
    avoided = [dict(response_h=float(r), avoided_kg=Q * max(T_BURN - r, 0) * 3.6) for r in resp]
    # plume map for figure (night)
    x = np.linspace(-200, 2000, 441)
    y = np.linspace(-500, 500, 201)
    X, Y = np.meshgrid(x, y)
    cn = conc(X, Y, Z_SENSOR, 1.0, CLASSES["night (stable, E-F)"]["a"]) * 1e6
    np.save(os.path.join(HERE, "plume_night.npy"), np.stack([X, Y, cn]))
    cd = conc(X, Y, Z_SENSOR, 2.5, CLASSES["day (neutral, D)"]["a"]) * 1e6
    np.save(os.path.join(HERE, "plume_day.npy"), np.stack([X, Y, cd]))
    # increment on the plume axis
    axis = {cls: {d: float(conc(np.array([float(d)]), np.array([0.0]), Z_SENSOR, CLASSES[cls]["u"],
                                CLASSES[cls]["a"])[0] * 1e6) for d in (100, 250, 500, 1000, 2000)}
            for cls in CLASSES}
    # network scenario for a 25 km2 district: 500 m grid with patrol drone
    s_pick = 500
    dn = next(d for d in det["night (stable, E-F)"] if d["spacing_m"] == s_pick)
    rt = response_time_h(dn["mean_dist_m"] + s_pick)
    out = dict(Q_g_s=Q, EF_g_kg=EF, burn_kg_h=BURN, threshold_ug_m3=THRESH * 1e6, sensor_height_m=Z_SENSOR,
               detection=det, axis_increment_ug_m3=axis, avoided=avoided,
               scenario=dict(spacing_m=s_pick, sensors_per_25km2=int((5000 / s_pick + 1) ** 2),
                             p_detect_night=dn["p_detect"], response_h=rt,
                             avoided_kg_per_fire=Q * max(T_BURN - rt, 0) * 3.6,
                             unattended_kg_per_fire=Q * T_BURN * 3.6))
    out["patrol"] = [patrol(k) for k in (0, 1, 2, 4)]
    json.dump(out, open(os.path.join(HERE, "results.json"), "w"), indent=1)
    for r in out["patrol"]:
        print("patrol", {k: round(v, 2) for k, v in r.items()})
    for cls in det:
        print(cls, [(d["spacing_m"], round(d["p_detect"], 2)) for d in det[cls]])
    print("axis", {k: {d: round(v) for d, v in a.items()} for k, a in axis.items()})
    print("scenario", out["scenario"])


if __name__ == "__main__":
    main()
