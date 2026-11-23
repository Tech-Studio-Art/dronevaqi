#!/usr/bin/env python3
"""Full-field PM2.5 tracer for figure-making: restart each filter case from its final
flow, reset the tracer to ambient, apply the (stable) PM2.5 filter sink, run 400 its,
export VTK_pm25. Flow is unchanged (beta = 0)."""
import os, re, sys, json
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import drone3d as d
from foam import run, last_time, vtk_zone

def go(name):
    case = os.path.join(d.RUNS, name)
    r = json.load(open(os.path.join(case, "result.json")))
    lam = next(h for h in r["history"] if h["pass_"] == "PM25")["lam"]
    fv = os.path.join(case, "constant/fvOptions")
    s = open(fv).read()
    s = re.sub(r"injectionRateSuSp \{ T \(0 -?[0-9.e+-]+\); \}", f"injectionRateSuSp {{ T (0 {-lam:.6g}); }}", s)
    open(fv, "w").write(s)
    fs = os.path.join(case, "system/fvSchemes")
    txt = open(fs).read().replace("div(phi,T) bounded Gauss limitedLinear 1;", "div(phi,T) bounded Gauss upwind;")
    open(fs, "w").write(txt)
    t = last_time(case)
    tp = os.path.join(case, t, "T")
    s = re.sub(r"internalField\s+nonuniform List<scalar>\s*\n\d+\s*\n\(.*?\n\)\s*;", "internalField   uniform 1;",
               open(tp).read(), flags=re.S)
    s = re.sub(r"value\s+nonuniform List<scalar>\s*\n\d+\s*\n\(.*?\n\)\s*;", "value uniform 1;", s, flags=re.S)
    open(tp, "w").write(s)
    cd = os.path.join(case, "system/controlDict")
    s = re.sub(r"endTime \d+;", f"endTime {int(float(t)) + 400};", open(cd).read())
    s = re.sub(r"writeInterval \d+;", "writeInterval 400;", s)
    open(cd, "w").write(s)
    run(case, "buoyantBoussinesqSimpleFoam", "solve.pm25field")
    vtk_zone(case, None, "VTK_pm25")
    st = d.stats(case, r["filter"], "pm25field")
    print(name, "PM2.5 C_avg in filter", st["T_avg"], "(was", next(h for h in r["history"] if h["pass_"] == "PM25")["T_avg"], ")")

if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=int(os.environ.get("NPROC", "2"))) as ex:
        list(ex.map(go, sys.argv[1:] or ["quad_above", "quad_below", "heavy_above", "heavy_below"]))
