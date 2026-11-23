#!/usr/bin/env python3
"""
Solution 2 in 3-D: bus shelter (8 m long, 2.5 m deep, roof at 2.5 m, back wall,
OPEN ends) with a filtered-air "shower" from the roof.

Same methods as cases2d.py (steady RANS k-epsilon, Richards & Hoxey ABL inlet,
passive scalar Sc_t = 0.7). Wind directions:
  parallel : wind blows along the shelter, straight through the open ends
             (worst case: mean flow carries ambient air through the shelter)
  into     : wind blows into the open front
Supply air: filtered to 5 % of ambient (HEPA/F9-class unit).
"""
import json
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from foam import write, run  # noqa: E402
from cases2d import common, load, spacing  # noqa: E402

RUNS = os.path.join(HERE, "runs", "shelter3d")


def lattice3d(case, xb, yb, zb, hx, hy, hz, solid, special):
    nx, ny, nz = len(xb) - 1, len(yb) - 1, len(zb) - 1
    vid, verts = {}, []
    for k, z in enumerate(zb):
        for j, y in enumerate(yb):
            for i, x in enumerate(xb):
                vid[(i, j, k)] = len(verts)
                verts.append((x, y, z))
    blocks, faces = [], {}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (i, j, k) in solid:
                    continue
                v = [vid[(i, j, k)], vid[(i + 1, j, k)], vid[(i + 1, j + 1, k)], vid[(i, j + 1, k)],
                     vid[(i, j, k + 1)], vid[(i + 1, j, k + 1)], vid[(i + 1, j + 1, k + 1)], vid[(i, j + 1, k + 1)]]
                blocks.append(f"hex ({' '.join(map(str, v))}) ({spacing(xb[i], xb[i + 1], hx[i])} "
                              f"{spacing(yb[j], yb[j + 1], hy[j])} {spacing(zb[k], zb[k + 1], hz[k])}) simpleGrading (1 1 1)")
                nb = {"xm": (i - 1, j, k), "xp": (i + 1, j, k), "ym": (i, j - 1, k), "yp": (i, j + 1, k),
                      "zm": (i, j, k - 1), "zp": (i, j, k + 1)}
                fd = {"xm": (v[0], v[4], v[7], v[3]), "xp": (v[1], v[2], v[6], v[5]),
                      "ym": (v[0], v[1], v[5], v[4]), "yp": (v[3], v[7], v[6], v[2]),
                      "zm": (v[0], v[3], v[2], v[1]), "zp": (v[4], v[5], v[6], v[7])}
                for side, (a, b, c) in nb.items():
                    out = not (0 <= a < nx and 0 <= b < ny and 0 <= c < nz)
                    if out or (a, b, c) in solid:
                        p = special(i, j, k, side)
                        if p is None:
                            p = ({"xm": "inlet", "xp": "outlet", "ym": "sides", "yp": "sides",
                                  "zm": "ground", "zp": "top"}[side] if out else "walls")
                        faces.setdefault(p, []).append(fd[side])
    types = {"ground": "wall", "walls": "wall"}
    bnd = "\n".join(f"  {p} {{ type {types.get(p, 'patch')}; faces ( "
                    + " ".join("(" + " ".join(map(str, f)) + ")" for f in fl) + " ); }" for p, fl in faces.items())
    write(case, "system/blockMeshDict", "convertToMeters 1;\nvertices (\n"
          + "\n".join(f"({x} {y} {z})" for x, y, z in verts) + "\n);\nblocks (\n" + "\n".join(blocks)
          + "\n);\nedges ();\nboundary (\n" + bnd + "\n);\nmergePatchPairs ();\n")
    return list(faces)


def shelter3d(Uref, w, wind):
    name = f"s3d_U{Uref}_w{w}_{wind}"
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    f = 0.1   # fine cell
    if wind == "parallel":
        # shelter length along x (flow direction), depth along y, back wall at y 2.5-2.6
        xb = [-25, -6, -1, 0, 0.5, 7.5, 8, 9, 20, 50]
        yb = [-15, -4, -1, 0, 0.3, 2.3, 2.5, 2.6, 4, 15]
        hx = [3, 1, 0.25, f, f, f, 0.25, 1, 3]
        hy = [3, 1, 0.25, f, f, f, f, 0.5, 3]
        xin = range(3, 6)          # i indices of the shelter footprint (x 0-8)
        yin = range(3, 7)          # j indices incl. back wall (y 0-2.6)
        wall = {(i, 6) for i in xin}
        sup = {(4, 4)}
        inside = dict(x=(0, 8), y=(0, 2.5))
    else:
        # wind into opening: depth along x (opening at x=0 faces upwind), length along y
        xb = [-25, -6, -1, 0, 0.3, 2.3, 2.5, 2.6, 4, 20, 50]
        yb = [-15, -4, -1, 0, 0.5, 7.5, 8, 9, 20]
        hx = [3, 1, 0.25, f, f, f, f, 0.5, 1, 3]
        hy = [3, 1, 0.25, f, f, f, 0.25, 1]
        xin = range(3, 7)
        yin = range(3, 6)
        wall = {(6, j) for j in yin}
        sup = {(4, 4)}
        inside = dict(x=(0, 2.5), y=(0, 8))
    zb = [0, 0.5, 1.8, 2.5, 2.6, 4, 8, 20]
    hz = [f, f, f, f, 0.25, 0.75, 2.5]
    assert len(hx) == len(xb) - 1 and len(hy) == len(yb) - 1
    solid = {(i, j, 3) for i in xin for j in yin}                   # roof slab
    solid |= {(i, j, k) for (i, j) in wall for k in (0, 1, 2)}       # back wall
    sup_ij = sup

    def special(i, j, k, side):
        if side == "zp" and k == 2 and (i, j) in sup_ij and w > 0:
            return "supply"
        return None
    patches = lattice3d(case, xb, yb, zb, hx, hy, hz, solid, special)
    common(case, patches, Uref, supply={"supply": ((0, -w), 0.05)} if w > 0 else None, n_iter=1500)
    write(case, "system/decomposeParDict", "numberOfSubdomains 1; method scotch;\n")
    run(case, "blockMesh")
    run(case, "buoyantBoussinesqSimpleFoam", "solve")
    import pyvista as pv
    from foam import vtk_zone
    g = vtk_zone(case, None, "VTK_full", fields=("U", "T"))
    cc = np.asarray(g.cell_centers().points)
    V = np.asarray(g.cell_data["Volume"])
    T = np.asarray(g.cell_data["T"])
    m = ((cc[:, 0] > inside["x"][0]) & (cc[:, 0] < inside["x"][1]) & (cc[:, 1] > inside["y"][0])
         & (cc[:, 1] < inside["y"][1]) & (cc[:, 2] > 0.5) & (cc[:, 2] < 1.8))
    bz = float((T[m] * V[m]).sum() / V[m].sum())
    Q = w * 7.0 * 2.0                    # supply area 7 m x 2 m
    return dict(case=name, Uref=Uref, w_supply=w, wind=wind, cells=int(len(V)), breathing_zone_C=bz,
                reduction_pct=100 * (1 - bz), Q_m3s=Q, fan_kW=Q * 350 / 0.5 / 1000,
                fan_kW_lowdp=Q * 150 / 0.5 / 1000)


JOBS = [(1.0, 0.25, "parallel"), (1.0, 0.5, "parallel"), (1.0, 1.0, "parallel"),
        (2.0, 0.5, "parallel"), (1.0, 0.5, "into")]

if __name__ == "__main__":
    os.makedirs(RUNS, exist_ok=True)
    nproc = int(os.environ.get("NPROC", "4"))
    jobs = JOBS if len(sys.argv) == 1 else [(float(sys.argv[1]), float(sys.argv[2]), sys.argv[3])]

    def go(a):
        try:
            return shelter3d(*a)
        except Exception as e:
            return dict(case=str(a), error=str(e))
    with ThreadPoolExecutor(max_workers=nproc) as ex:
        res = list(ex.map(go, jobs))
    json.dump(res, open(os.path.join(RUNS, "summary.json" if len(sys.argv) == 1 else "test.json"), "w"), indent=1)
    for r in res:
        print(r)
