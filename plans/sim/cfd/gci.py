#!/usr/bin/env python3
"""
Grid-independence study (step 1a).

Three meshes with a constant refinement ratio r = sqrt(2) (coarse, medium = base mesh,
fine). The discretisation uncertainty is reported as the Grid Convergence Index (GCI)
following the procedure of Celik et al. (2008), J. Fluids Eng. 130:078001:
  eps21 = f2 - f1, eps32 = f3 - f2 (1 = fine, 2 = medium, 3 = coarse)
  apparent order   p = |ln|eps32/eps21|| / ln r        (constant r)
  extrapolated     f_ext = (r^p f1 - f2) / (r^p - 1)
  GCI_fine         = 1.25 |(f1 - f2)/f1| / (r^p - 1)
Oscillatory convergence (eps32/eps21 < 0) is flagged; then p is not defined and the
GCI is given with p = 1 (conservative) and a safety factor of 3.

  python3 gci.py fence    # run the fence meshes (2-D, minutes)
  python3 gci.py report   # collect fence + inversion-mixing results, write runs/gci/summary.json
"""
import json
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
R = math.sqrt(2)
OUT = os.path.join(HERE, "runs", "gci")
LEVELS = (("coarse", 1 / R), ("medium", 1.0), ("fine", R))


def gci(f1, f2, f3, r=R):
    e21, e32 = f2 - f1, f3 - f2
    res = dict(fine=f1, medium=f2, coarse=f3)
    if abs(e21) < 1e-12:
        return dict(res, p=None, f_ext=f1, gci_fine_pct=0.0, note="fine = medium")
    ratio = e32 / e21
    if ratio < 0:
        g = 3.0 * abs((f1 - f2) / f1) / (r - 1)
        return dict(res, p=None, f_ext=None, gci_fine_pct=100 * g, note="oscillatory convergence")
    p = abs(math.log(abs(ratio))) / math.log(r)
    note = "monotonic convergence" if abs(e21) < abs(e32) else "monotonic divergence"
    if abs(e21) >= abs(e32):
        g = 3.0 * abs((f1 - f2) / f1) / (r - 1)
        return dict(res, p=p, f_ext=None, gci_fine_pct=100 * g, note=note)
    fext = (r ** p * f1 - f2) / (r ** p - 1)
    g = 1.25 * abs((f1 - f2) / f1) / (r ** p - 1)
    return dict(res, p=p, f_ext=fext, gci_fine_pct=100 * g, note=note)


def fence():
    import cases2d
    import esp_road
    out = {}
    for tag, f in LEVELS:
        cases2d.MESH_F = f
        cases2d.RUNS = os.path.join(OUT, f"fence_{tag}")
        esp_road.RUNS = os.path.join(OUT, f"road_{tag}")
        os.makedirs(cases2d.RUNS, exist_ok=True)
        os.makedirs(esp_road.RUNS, exist_ok=True)
        jobs = [("amb", None), ("none", "road"), ("esp", "bg"), ("esp", "road")]

        def go(j):
            if j[0] == "amb":
                return "amb", cases2d.esp_fence(1.0, "PM25")
            return f"{j[0]}_{j[1]}", esp_road.case(*j)[1]
        with ThreadPoolExecutor(max_workers=int(os.environ.get("NPROC", "2"))) as ex:
            res = dict(ex.map(go, jobs))
        amb = res["amb"]["pedestrian_C"]
        s = 0.4 / res["none_road"][2.5]
        cn = {x: 1 + s * res["none_road"][x] for x in (2.5, 40)}
        ce = {x: res["esp_bg"][x] + s * res["esp_road"][x] for x in (2.5, 40)}
        out[tag] = dict(amb_red_2m=100 * (1 - amb["x=2m"]), amb_red_40m=100 * (1 - amb["x=40m"]),
                        road_red_2p5m=100 * (1 - ce[2.5] / cn[2.5]), road_red_40m=100 * (1 - ce[40] / cn[40]),
                        single_pass_eff=100 * res["amb"]["single_pass_eff"])
        print(tag, {k: round(v, 2) for k, v in out[tag].items()}, flush=True)
    cases2d.MESH_F = 1.0
    json.dump(out, open(os.path.join(OUT, "fence_levels.json"), "w"), indent=1)


def mixing_levels():
    import mixing_axi as mx
    out = {}
    for tag, name in (("coarse", "heavy_z60_gc"), ("medium", "heavy_z60"), ("fine", "heavy_z60_gf")):
        p = os.path.join(mx.RUNS, name, "result.json")
        if not os.path.exists(p):
            return None
        r = json.load(open(p))
        s = r["series"]
        out[tag] = dict(cells=r["cells"],
                        red_10min=100 * (1 - s[2]["z<2m"] / s[0]["z<2m"]),
                        red_30min=100 * (1 - s[-1]["z<2m"] / s[0]["z<2m"]),
                        red_30min_0_10m=100 * (1 - s[-1]["z<10m"] / s[0]["z<10m"]))
    return out


def report():
    summary = {}
    fl = os.path.join(OUT, "fence_levels.json")
    if os.path.exists(fl):
        L = json.load(open(fl))
        summary["fence"] = {k: gci(L["fine"][k], L["medium"][k], L["coarse"][k]) for k in L["medium"]
                            if all(L[t][k] == L[t][k] for t in ("fine", "medium", "coarse"))}   # skip NaN (no cells in window)
    ml = mixing_levels()
    if ml:
        summary["mixing"] = {k: gci(ml["fine"][k], ml["medium"][k], ml["coarse"][k]) for k in ml["medium"] if k != "cells"}
        summary["mixing_cells"] = {k: v["cells"] for k, v in ml.items()}
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    for grp, d in summary.items():
        if isinstance(d, dict) and grp != "mixing_cells":
            for k, v in d.items():
                print(f"{grp:7s} {k:18s} coarse {v['coarse']:7.2f}  medium {v['medium']:7.2f}  fine {v['fine']:7.2f}  "
                      f"p {v['p'] if v['p'] is None else round(v['p'], 2)}  GCI {v['gci_fine_pct']:.2f} %  {v['note']}")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    {"fence": fence, "report": report}[sys.argv[1]]()
