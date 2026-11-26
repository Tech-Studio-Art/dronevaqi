#!/usr/bin/env python3
"""
ESP wind fence beside a road (2-D RANS, same methods as cases2d.py).

Geometry (wind blows from the road towards the footpath):
  road (2 lanes, x = -12 ... -4 m) | verge | fence (x = 0 ... 0.5 m, 6 m high) | footpath (x = 1 ... 4 m)
Traffic emission = uniform volume source in the lowest 1.5 m over the road
(line-source treatment used in near-road CFD, e.g. Baldauf et al. 2008;
Hagler et al. 2011). Pollutant is linear, so background and road parts are run
separately and superposed:
    C = C_background + s * C_road(unit source)
with s chosen so that WITHOUT a fence the footpath is 40 % above background
(an assumed, typical near-road increment; stated in the report).
Cases: no fence, passive porous fence (same K, no charge), ESP fence.
"""
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from foam import run, write  # noqa: E402
from cases2d import lattice_mesh, common, porous_fv, load, region_mean, ESP  # noqa: E402

RUNS = os.path.join(HERE, "runs", "esp_road")
U = 1.0


def case(mode, part, w25=ESP["w25"]):
    name = f"{mode}_{part}" + ("" if w25 == ESP["w25"] else f"_w{w25}")
    c = os.path.join(RUNS, name)
    shutil.rmtree(c, ignore_errors=True)
    os.makedirs(c)
    H, L = ESP["H"], ESP["L"]
    xb = [-60, -16, -12, -4, -2, 0, L, 1, 4, 15, 60, 150]
    zb = [0, 1.5, H, H + 2, 15, 60]
    hx = [3, 1, 0.25, 0.25, 0.25, 0.1, 0.1, 0.25, 0.5, 1.5, 3]
    hz = [0.1, 0.2, 0.25, 0.75, 3]
    patches = lattice_mesh(c, xb, zb, hx, hz)
    acts = [f"""{{ name road; type cellSet; action new; source boxToCell; sourceInfo {{ box (-12 -1 0) (-4 2 1.5); }} }}
 {{ name roadZone; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set road; }} }}"""]
    if mode != "none":
        acts.append(f"""{{ name fence; type cellSet; action new; source boxToCell; sourceInfo {{ box (0 -1 0) ({L} 2 {H}); }} }}
 {{ name fenceZone; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set fence; }} }}""")
    write(c, "system/topoSetDict", "actions (\n" + "\n".join(acts) + "\n);\n")
    fv = ""
    if mode != "none":
        lam = ESP["a"] * w25 if mode == "esp" else 0.0
        fv = porous_fv("fenceZone", ESP["K"], L, lam)
    if part == "road":
        fv += """road
{
    type scalarSemiImplicitSource; active yes;
    scalarSemiImplicitSourceCoeffs
    {
        timeStart 0; duration 1e9; selectionMode cellZone; cellZone roadZone; volumeMode specific;
        injectionRateSuSp { T (0.01 0); }
    }
}
"""
    common(c, patches, U, fvoptions=fv)
    if part == "road":   # background zero for the road-only part
        p = os.path.join(c, "0", "T")
        s = open(p).read().replace("internalField   uniform 1;", "internalField   uniform 0;")
        s = s.replace("type fixedValue; value uniform 1;", "type fixedValue; value uniform 0;")
        s = s.replace("inletValue uniform 1; value uniform 1;", "inletValue uniform 0; value uniform 0;")
        open(p, "w").write(s)
    run(c, "blockMesh")
    run(c, "topoSet")
    run(c, "buoyantBoussinesqSimpleFoam", "solve")
    d = load(c)
    prof = {x: region_mean(d, x - 0.5, x + 0.5, 1.0, 2.0) for x in (2.5, 5, 10, 20, 40)}
    return name, prof


if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor
    os.makedirs(RUNS, exist_ok=True)
    jobs = [("none", "road"), ("passive", "road"), ("esp", "road"), ("passive", "bg"), ("esp", "bg"),
            ("esp", "road", ESP["w25"] / 2), ("esp", "bg", ESP["w25"] / 2)]
    with ThreadPoolExecutor(max_workers=int(os.environ.get("NPROC", "3"))) as ex:
        res = dict(ex.map(lambda j: case(*j), jobs))
    res["none_bg"] = {x: 1.0 for x in (2.5, 5, 10, 20, 40)}
    s = 0.4 / res["none_road"][2.5]           # footpath +40 % without a fence
    out = {"scale": s, "raw": {k: {str(x): v for x, v in p.items()} for k, p in res.items()}, "total": {}}
    for mode, bg, rd in (("no fence", "none_bg", "none_road"), ("passive fence", "passive_bg", "passive_road"),
                         ("ESP fence", "esp_bg", "esp_road"), ("ESP fence, half migration velocity", "esp_bg_w0.02", "esp_road_w0.02")):
        out["total"][mode] = {str(x): res[bg][x] + s * res[rd][x] for x in (2.5, 5, 10, 20, 40)}
    json.dump(out, open(os.path.join(RUNS, "summary.json"), "w"), indent=1)
    for k, v in out["total"].items():
        print(k, {x: round(c * 300) for x, c in v.items()}, "ug/m3 PM2.5")
