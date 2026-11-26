#!/usr/bin/env python3
"""
Validation (step 1b-i): penetration of the drone's downward jet into stable air.

One drone hovers at 220 m in a deep, uniformly stratified layer (dtheta/dz = 0.05 K/m,
the same N as the 5 K / 100 m inversion), far from the ground, so the jet stops freely.
Same solver, turbulence model and actuator disk as mixing_axi.py; large isolated cell
(R = 150 m). The penetration depth is where the time-mean centreline velocity changes from
downward to zero. It is compared with the Morton-Taylor-Turner (1956) entrainment model
(mtt.py) over the usual range of jet entrainment coefficients.
"""
import glob
import json
import math
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mixing_axi as mx  # noqa: E402
from foam import run, write  # noqa: E402
import mtt  # noqa: E402

RUNS = os.path.join(HERE, "runs", "validation", "jetpen")
ZH, H, R = 220.0, 300.0, 150.0
CASES = {"heavy": dict(T=50 * 9.81, Rd=math.sqrt(4 * 1.474 / math.pi), A=4 * 1.474, zh=ZH),
         "quad": dict(T=2.0 * 9.81, Rd=0.75, A=4 * 0.0507, zh=ZH)}


def mesh(case, c, model):
    a = math.radians(mx.WEDGE / 2)
    cs, sn = math.cos(a), math.sin(a)
    r0, dr = 0.02, c["Rd"] / 4
    r1 = 10 * c["Rd"]
    n1 = int(math.ceil((r1 - r0) / dr))
    dz = 0.5 if c["T"] > 100 else 0.35
    zb0 = ZH - (90 if c["T"] > 100 else 55)
    zb1 = ZH + 10
    nzb = int(math.ceil((zb1 - zb0) / dz))
    v = []
    for z in (0, H):
        v += [(r0 * cs, -r0 * sn, z), (R * cs, -R * sn, z), (R * cs, R * sn, z), (r0 * cs, r0 * sn, z)]
    vs = " ".join(f"({x:.6f} {y:.6f} {z:.6f})" for x, y, z in v)
    nr = n1 + 45
    nz = 40 + nzb + 30
    gr = f"(({(r1 - r0) / (R - r0):.6f} {n1 / nr:.6f} 1) ({(R - r1) / (R - r0):.6f} {45 / nr:.6f} 10))"
    gz = f"(({zb0 / H:.6f} {40 / nz:.6f} 0.12) ({(zb1 - zb0) / H:.6f} {nzb / nz:.6f} 1) ({(H - zb1) / H:.6f} {30 / nz:.6f} 8))"
    write(case, "system/blockMeshDict", f"""
convertToMeters 1;
vertices ( {vs} );
blocks ( hex (0 1 2 3 4 5 6 7) ({nr} 1 {nz}) simpleGrading ({gr} 1 {gz}) );
boundary (
  ground {{ type wall; faces ( (0 3 2 1) ); }}
  top    {{ type wall; faces ( (4 5 6 7) ); }}
  outer  {{ type symmetry; faces ( (1 2 6 5) ); }}
  axis   {{ type symmetry; faces ( (0 4 7 3) ); }}
  front  {{ type wedge; faces ( (0 1 5 4) ); }}
  back   {{ type wedge; faces ( (3 7 6 2) ); }}
);
""")
    write(case, "system/topoSetDict", f"""
actions (
 {{ name disk; type cellSet; action new; source cylinderToCell;
   sourceInfo {{ p1 (0 0 {ZH - dz}); p2 (0 0 {ZH + dz}); radius {c['Rd']}; }} }}
 {{ name diskZone; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set disk; }} }}
);
""")


def solve(name, t_end=900.0):
    c = CASES[name]
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    mx.H_DOM, mx.INV_H, mx.D_THETA, mx.R_CELL, mx.T_END = H, 1000.0, 50.0, R, t_end
    mx.build(case, c)
    mesh(case, c, None)
    cd = os.path.join(case, "system/controlDict")
    s = open(cd).read().replace("writeInterval 300;", "writeInterval 60;")
    open(cd, "w").write(s)
    run(case, "blockMesh")
    run(case, "topoSet")
    mx.set_initial_theta(case)
    run(case, "buoyantBoussinesqPimpleFoam", "solve")
    return analyse(name)


def analyse(name, t_avg=(600, 900)):
    import pyvista as pv
    c = CASES[name]
    case = os.path.join(RUNS, name)
    run(case, "foamToVTK -no-point-data -no-boundary -overwrite -name VTK_all -fields '(U)'", "vtk.all")
    times = sorted(float(d) for d in os.listdir(case) if d.replace(".", "").isdigit() and float(d) > 0)
    files = sorted(glob.glob(os.path.join(case, "VTK_all", "*", "internal.vtu")),
                   key=lambda f: int(f.split("_")[-1].split("/")[0]))[1:]
    prof, depth_t = [], []
    for f, t in zip(files, times):
        g = pv.read(f)
        cc = np.asarray(g.cell_centers().points)
        r = np.hypot(cc[:, 0], cc[:, 1])
        m = r < 0.5 * c["Rd"]
        z, w = cc[m, 2], np.asarray(g.cell_data["U"])[m, 2]
        o = np.argsort(z)
        z, w = z[o], w[o]
        below = z < ZH - 3
        wref = -w[below][np.argmin(abs(z[below] - (ZH - 5)))]
        moving = below & (w < -0.05 * wref)
        zstop = z[moving].min() if moving.any() else ZH
        depth_t.append((t, ZH - zstop))
        if t_avg[0] <= t <= t_avg[1]:
            prof.append(w)
    wm = np.mean(prof, axis=0)
    below = z < ZH - 3
    wref = -wm[below][np.argmin(abs(z[below] - (ZH - 5)))]
    moving = below & (wm < -0.05 * wref)
    depth_mean = ZH - z[moving].min()
    rho = 1.25
    N = math.sqrt(9.81 / mx.THETA0 * 0.05)
    M0 = c["T"] / rho
    vi = math.sqrt(c["T"] / (2 * rho * c["A"]))
    Q0 = c["A"] * vi
    L = (M0 / N ** 2) ** 0.25
    mttr = {a: mtt.penetration(M0, Q0, N, a) for a in (0.05, 0.065, 0.076, 0.09, 0.1)}
    # entrainment coefficient that makes MTT match the CFD depth
    lo, hi = 0.01, 0.3
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if mtt.penetration(M0, Q0, N, mid) > depth_mean:
            lo = mid
        else:
            hi = mid
    res = dict(case=name, depth_cfd_mean_m=float(depth_mean), depth_max_over_time_m=float(max(d for _, d in depth_t)),
               depth_vs_time=depth_t, C_cfd=float(depth_mean / L), length_scale_m=L,
               mtt_depth_by_alpha={str(k): v for k, v in mttr.items()}, alpha_matching_cfd=0.5 * (lo + hi))
    json.dump(res, open(os.path.join(case, "result.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor
    names = sys.argv[1:] or list(CASES)
    os.makedirs(RUNS, exist_ok=True)
    with ThreadPoolExecutor(max_workers=2) as ex:
        for r in ex.map(solve, names):
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != "depth_vs_time"})
