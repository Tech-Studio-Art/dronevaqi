#!/usr/bin/env python3
"""
24x7 pollution-source surveillance of a waste-burning hotspot.

Scenario (from the Drone V AQI team): small open fires of domestic waste, each
burning 100 g of waste per hour for 4 hours; up to 100 such fires in one
500 m x 500 m hotspot.

Methods (standard, from the literature):
  * Dispersion: steady Gaussian plume, ground-level sources, ground reflection,
    Briggs urban dispersion coefficients (Hanna, Briggs & Hosker 1982);
    the hotspot field is the superposition of all fires.
      stability D (day):   sy = 0.16x(1+0.0004x)^-1/2,  sz = 0.14x(1+0.0003x)^-1/2, u = 2.5 m/s
      stability E-F (night): sy = 0.11x(1+0.0004x)^-1/2, sz = 0.08x(1+0.0015x)^-1/2, u = 1.0 m/s
  * Emission factor ~10 g PM2.5 per kg waste (Wiedinmyer et al. 2014)
    -> 1 g/h PM2.5 per fire, 100 g/h per hotspot.
  * Sensors: low-cost PM sensors (Snyder et al. 2013) on a square grid, 3 m high.
    One small fire is far below sensor noise; the network detects the HOTSPOT:
    an alert when the 1-hour mean increment at a sensor exceeds a threshold
    (10 or 25 ug/m3 above neighbouring sensors).
  * Response: patrol drones with thermal cameras map the fires in the alerted
    hotspot; ground teams stop the fires one by one.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(42)

EF = 10.0                # g PM2.5 / kg waste
BURN_G_H = 100.0         # g waste / h per fire
Q1 = EF * BURN_G_H / 1000 / 3600   # g/s PM2.5 per fire
N_FIRES = 100
ZONE = 500.0
T_BURN = 4.0             # h without action
Z_SENSOR = 3.0
CLASSES = {"night (stable, E-F)": dict(u=1.0, a=(0.11, 0.0004, 0.08, 0.0015)),
           "day (neutral, D)": dict(u=2.5, a=(0.16, 0.0004, 0.14, 0.0003))}


def conc_one(xd, yc, z, u, a, Q, Hs=1.0):
    out = np.zeros_like(xd, dtype=float)
    m = xd > 1.0
    x = xd[m]
    sy = a[0] * x / np.sqrt(1 + a[1] * x)
    sz = a[2] * x / np.sqrt(1 + a[3] * x)
    out[m] = (Q / (2 * np.pi * u * sy * sz) * np.exp(-yc[m] ** 2 / (2 * sy ** 2))
              * (np.exp(-(z - Hs) ** 2 / (2 * sz ** 2)) + np.exp(-(z + Hs) ** 2 / (2 * sz ** 2))))
    return out


def hotspot_field(px, py, fx, fy, theta, cls):
    """ug/m3 increment at points (px,py) from fires (fx,fy), wind towards angle theta."""
    a, u = CLASSES[cls]["a"], CLASSES[cls]["u"]
    c = np.zeros_like(px, dtype=float)
    ct, st = np.cos(theta), np.sin(theta)
    for x0, y0 in zip(fx, fy):
        dx, dy = px - x0, py - y0
        c += conc_one(dx * ct + dy * st, -dx * st + dy * ct, Z_SENSOR, u, a, Q1)
    return c * 1e6


def fires(n=N_FIRES):
    return rng.uniform(-ZONE / 2, ZONE / 2, n), rng.uniform(-ZONE / 2, ZONE / 2, n)


def field_stats():
    out = {}
    for cls in CLASSES:
        fx, fy = fires()
        x = np.linspace(-600, 2000, 261)
        y = np.linspace(-700, 700, 141)
        X, Y = np.meshgrid(x, y)
        C = hotspot_field(X, Y, fx, fy, 0.0, cls)
        inz = (np.abs(X) <= ZONE / 2) & (np.abs(Y) <= ZONE / 2)
        dw = lambda d: float(C[(np.abs(X - (ZONE / 2 + d)) < 6) & (np.abs(Y) < 50)].mean())
        out[cls] = dict(zone_mean=float(C[inz].mean()), zone_p90=float(np.percentile(C[inz], 90)),
                        downwind_250m=dw(250), downwind_1km=dw(1000))
        if cls.startswith("night"):
            np.save(os.path.join(HERE, "hotspot_night.npy"), np.stack([X, Y, C]))
            np.save(os.path.join(HERE, "hotspot_fires.npy"), np.stack([fx, fy]))
    return out


def detect_prob(spacing, cls, thresh, n=400):
    k = np.arange(-8, 9) * spacing
    gx, gy = np.meshgrid(k, k)
    gx, gy = gx.ravel(), gy.ravel()
    hit = 0
    for _ in range(n):
        ox, oy = rng.uniform(0, spacing, 2)
        th = rng.uniform(0, 2 * np.pi)
        fx, fy = fires()
        c = hotspot_field(gx + ox, gy + oy, fx, fy, th, cls)
        hit += c.max() > thresh
    return hit / n


def response(n_drones, n_teams, p_sensor, area_km2=25.0, swath_km=0.3, speed_kmh=36.0, duty=0.7,
             t_sensor_h=1.0, map_h=0.1, min_per_fire_drone=6.0, min_per_fire_foot=20.0, n=4000):
    """Emission prevented in one hotspot. Detection: first of (sensor alert at t_sensor_h, if the
    network sees the hotspot) and (patrol drone pass, uniform in its revisit interval). With a drone,
    fires are mapped (map_h) and teams go straight to each fire; without a drone, teams search on foot."""
    R = area_km2 / (n_drones * swath_km * speed_kmh * duty) if n_drones else np.inf
    t_dr = rng.uniform(0, R, n) if n_drones else np.full(n, np.inf)
    t_se = np.where(rng.uniform(0, 1, n) < p_sensor, t_sensor_h, np.inf)
    t_det = np.minimum(t_dr, t_se)
    per_fire_h = (min_per_fire_drone if n_drones else min_per_fire_foot) / 60
    start = t_det + (map_h if n_drones else 0)
    k = np.arange(1, N_FIRES + 1)
    t_stop = start[:, None] + np.ceil(k / n_teams)[None, :] * per_fire_h
    burned = np.minimum(t_stop, T_BURN) / T_BURN
    return dict(drones=n_drones, teams=n_teams, revisit_h=float(R), p_found=float((t_det < T_BURN).mean()),
                mean_detect_h=float(np.mean(np.minimum(t_det, T_BURN))),
                emission_prevented_pct=float(100 * (1 - burned.mean())))


def main():
    out = dict(Q_per_fire_g_h=Q1 * 3600, n_fires=N_FIRES, zone_m=ZONE, burn_h=T_BURN,
               hotspot_g_h=Q1 * 3600 * N_FIRES, hotspot_g_per_night=Q1 * 3600 * N_FIRES * T_BURN)
    out["single_fire_100m_night_ug_m3"] = float(conc_one(np.array([100.0]), np.array([0.0]), Z_SENSOR, 1.0,
                                                          CLASSES["night (stable, E-F)"]["a"], Q1)[0] * 1e6)
    out["field"] = field_stats()
    out["detection"] = {}
    for cls in CLASSES:
        for th in (10.0, 25.0):
            out["detection"][f"{cls}|{th:g}"] = [dict(spacing_m=s, p=detect_prob(s, cls, th)) for s in (250, 500, 1000, 2000)]
    p500 = next(d["p"] for d in out["detection"]["night (stable, E-F)|10"] if d["spacing_m"] == 500)
    out["response"] = [response(d, t, p500) for t in (2, 4) for d in (0, 1, 2, 4)]
    json.dump(out, open(os.path.join(HERE, "results.json"), "w"), indent=1)
    print("per fire g/h", out["Q_per_fire_g_h"], "single fire @100 m night", round(out["single_fire_100m_night_ug_m3"], 2))
    print("field", json.dumps(out["field"], indent=0))
    for k, v in out["detection"].items():
        print(k, [(d["spacing_m"], round(d["p"], 2)) for d in v])
    for r in out["response"]:
        print({k: round(v, 2) for k, v in r.items()})


if __name__ == "__main__":
    main()
