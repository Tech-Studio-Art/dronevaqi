#!/usr/bin/env python3
"""
Validation / verification of the fence model (step 1b-iii).

(a) Air flow through a porous screen. Actuator (momentum) theory for a screen with
    pressure-loss coefficient K gives the through-flow fraction
        u_s / U = 1 / (1 + K/4)
    (Taylor's screen theory; reviewed by Laws & Livesey 1978). We run the 2-D fence with
    K = 0.5, 1.5, 3, 6 and no charge, and compare the mean through-flow velocity with the
    velocity in the same zone without a fence (K = 0).
(b) Particle capture inside the fence. The Deutsch-Anderson law gives the single-pass
    penetration exp(-lambda L / u_s). We compare it with the CFD concentration drop across
    the fence (medium mesh of the grid study, PM2.5, wind 1 m/s).
"""
import json
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cases2d  # noqa: E402

OUT = os.path.join(HERE, "runs", "validation", "fence")


def through(K):
    cases2d.RUNS = os.path.join(OUT, f"K{K:g}")
    os.makedirs(cases2d.RUNS, exist_ok=True)
    cases2d.ESP["K"] = K
    cases2d.ESP["w25"] = 0.0
    return K, cases2d.esp_fence(1.0, "PM25")["u_through"]


def main():
    os.makedirs(OUT, exist_ok=True)
    Ks = [0.0, 0.5, 1.5, 3.0, 6.0]
    with ThreadPoolExecutor(max_workers=2) as ex:
        res = dict(ex.map(through, Ks))
    u0 = res[0.0]
    a = [dict(K=K, cfd=res[K] / u0, theory=1 / (1 + K / 4)) for K in Ks if K > 0]
    # (b) capture check on the medium grid-study case
    case = os.path.join(HERE, "runs", "gci", "fence_medium", "espfence_U1.0_PM25")
    d = cases2d.load(case)
    H, L = 6.0, 0.5
    lam = 20.0 * 0.04
    zsel = (d["z"] > 0.5) & (d["z"] < H - 0.5)
    up = zsel & (d["x"] > -0.15) & (d["x"] < 0)
    dn = zsel & (d["x"] > L) & (d["x"] < L + 0.15)
    inside = zsel & (d["x"] > 0) & (d["x"] < L)
    us = float((d["U"][inside, 0] * d["V"][inside]).sum() / d["V"][inside].sum())
    c_up = float((d["T"][up] * d["V"][up]).sum() / d["V"][up].sum())
    c_dn = float((d["T"][dn] * d["V"][dn]).sum() / d["V"][dn].sum())
    b = dict(u_through=us, cfd_penetration=c_dn / c_up, deutsch_anderson=math.exp(-lam * L / us))
    out = dict(screen=a, capture=b)
    json.dump(out, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    for r in a:
        print(f"K={r['K']:4.1f}  through-flow CFD {r['cfd']:.3f}  theory {r['theory']:.3f}  diff {100 * (r['cfd'] / r['theory'] - 1):+.1f} %")
    print(f"capture: CFD penetration {b['cfd_penetration']:.3f}  Deutsch-Anderson {b['deutsch_anderson']:.3f}")


if __name__ == "__main__":
    main()
