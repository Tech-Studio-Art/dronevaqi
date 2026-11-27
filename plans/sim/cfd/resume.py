#!/usr/bin/env python3
"""Resume an interrupted inversion-mixing case from its last saved time, then analyse it.
Usage: python3 resume.py <case dir relative to runs/> [<case dir> ...]"""
import json
import math
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mixing_axi as mx  # noqa: E402
from foam import run  # noqa: E402


def resume(rel):
    case = os.path.join(HERE, "runs", rel)
    name = os.path.basename(case)
    app = re.search(r"application\s+(\w+);", open(os.path.join(case, "system/controlDict")).read()).group(1)
    n = sum(1 for f in os.listdir(case) if f.startswith("log.solve"))
    run(case, app, f"solve.resume{n}")
    cells = int(re.search(r"nCells:\s*(\d+)", open(os.path.join(case, "log.blockMesh")).read()).group(1))
    key = name.replace("_cg1", "")
    c = mx.CASES.get(key, mx.CASES.get(name, {}))
    res = dict(case=name, cells=cells, thrust_N=c.get("T"), hover_z=c.get("zh"), series=mx.analyse(case))
    if "mixing_buoy" in rel:
        res["Cg"] = float(re.search(r"Cg\s+([0-9.]+);", open(os.path.join(case, "constant/turbulenceProperties")).read()).group(1))
    elif c.get("T"):
        res.update(disk_R=c["Rd"], jet_penetration_theory_m=mx.jet_penetration(c["T"]))
    json.dump(res, open(os.path.join(case, "result.json"), "w"), indent=1)
    s = res["series"]
    return name, [round(100 * x["z<2m"] / s[0]["z<2m"], 1) for x in s]


if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=len(sys.argv) - 1) as ex:
        for r in ex.map(resume, sys.argv[1:]):
            print(*r)
