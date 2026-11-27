#!/usr/bin/env python3
"""Rewrite the numbers and the performance chart of plans/docs/inversion-drones.html from the
CFD results (best estimate = calibrated buoyant k-epsilon; upper limit = standard k-epsilon)."""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PLANS = os.path.dirname(os.path.dirname(HERE))
CFD = os.path.join(PLANS, "sim", "cfd")
PAGE = os.path.join(PLANS, "docs", "inversion-drones.html")

mb = json.load(open(os.path.join(CFD, "runs", "mixing_buoy", "summary.json")))
ms = json.load(open(os.path.join(CFD, "runs", "mixing", "summary.json")))


def rel(series):
    return [100 * x["z<2m"] / series[0]["z<2m"] for x in series]


def pts(vals):
    # chart: x 0..30 min -> 60..540 px; y 76..100 % -> 260..20 px (10 px per %)
    return " ".join(f"{60 + 80 * i},{20 + (100 - v) * 10:.0f}" for i, v in enumerate(vals))


best = {k: rel(mb[k]["series"]) for k in ("no_drone", "quad_z30", "heavy_z60", "heavy_z90") if k in mb}
if "heavy_z90" not in best:
    best["heavy_z90"] = rel(ms["heavy_z90"]["series"])
upper = rel(ms["heavy_z60"]["series"])
d60, d60u, dq = 100 - best["heavy_z60"][-1], 100 - upper[-1], 100 - best["quad_z30"][-1]

s = open(PAGE).read()
# polylines (order in the page: no drones s4, quads s2, heavy 90 s3, heavy 60 s1)
lines = {"--s4": best["no_drone"], "--s2": best["quad_z30"], "--s3": best["heavy_z90"], "--s1": best["heavy_z60"]}
for var, vals in lines.items():
    s, n = re.subn(rf'(<polyline fill="none" stroke="var\({var}\)" stroke-width="[0-9.]+" points=")[^"]*(")',
                   lambda m: m.group(1) + pts(vals) + m.group(2), s, count=1)
    assert n == 1, var
dashed = f'<polyline class="upper" fill="none" stroke="var(--s1)" stroke-width="1.6" stroke-dasharray="5 4" points="{pts(upper)}"/>'
s = re.sub(r'\s*<polyline class="upper"[^>]*/>', "", s)
s = s.replace('<polyline fill="none" stroke="var(--s1)"', dashed + '\n        <polyline fill="none" stroke="var(--s1)"', 1)
# end markers and labels
s = re.sub(r'<circle cx="540" cy="\d+" r="4" fill="var\(--s1\)"', f'<circle cx="540" cy="{20 + (100 - best["heavy_z60"][-1]) * 10:.0f}" r="4" fill="var(--s1)"', s)
s = re.sub(r'<circle cx="540" cy="\d+" r="4" fill="var\(--s3\)"', f'<circle cx="540" cy="{20 + (100 - best["heavy_z90"][-1]) * 10:.0f}" r="4" fill="var(--s3)"', s)
s = re.sub(r'<circle cx="540" cy="\d+" r="4" fill="var\(--s2\)"', f'<circle cx="540" cy="{20 + (100 - best["quad_z30"][-1]) * 10:.0f}" r="4" fill="var(--s2)"', s)
s = re.sub(r'<text x="552" y="\d+" class="tk">quads</text>', f'<text x="552" y="{24 + (100 - best["quad_z30"][-1]) * 10:.0f}" class="tk">quads</text>', s)
s = re.sub(r'<text x="552" y="\d+" class="tk">heavy 60 m</text>', f'<text x="552" y="{20 + (100 - best["heavy_z60"][-1]) * 10 - 6:.0f}" class="tk">heavy 60 m</text>', s)
s = re.sub(r'<text x="552" y="\d+" class="tk">heavy 90 m</text>', f'<text x="552" y="{20 + (100 - best["heavy_z90"][-1]) * 10 + 14:.0f}" class="tk">heavy 90 m</text>', s)
# legend
s = re.sub(r'(<span><i style="background:var\(--s1\)"></i>)100 heavy @ 60 m: [0-9.]+ %(</span>)',
           rf'\g<1>100 heavy @ 60 m: {best["heavy_z60"][-1]:.1f} % (dashed: standard k-ε {upper[-1]:.1f} %)\g<2>', s)
s = re.sub(r'(<span><i style="background:var\(--s3\)"></i>)100 heavy @ 90 m: [0-9.]+ %(</span>)',
           rf'\g<1>100 heavy @ 90 m: {best["heavy_z90"][-1]:.1f} %\g<2>', s)
s = re.sub(r'(<span><i style="background:var\(--s2\)"></i>)100 quads @ 30 m: [0-9.]+ %(</span>)',
           rf'\g<1>100 quads @ 30 m: {best["quad_z30"][-1]:.1f} %\g<2>', s)
s = re.sub(r'(<span><i style="background:var\(--s4\)"></i>)no drones: [0-9.]+ %(</span>)',
           rf'\g<1>no drones: {best["no_drone"][-1]:.1f} %\g<2>', s)
# headline numbers and text
s = re.sub(r'Breathing-zone PM after 30 min, 100 heavy drones: <b>[^<]*</b>',
           f'Breathing-zone PM after 30 min, 100 heavy drones: <b>−{d60:.0f} % (at 60 m) to −{100 - best["heavy_z90"][-1]:.0f} % (at 90 m)</b>', s)
s = re.sub(r'In the CFD, 100 quads cut breathing-zone PM by [0-9.]+ % in 30 minutes against [0-9.]+ % for 100 heavy drones\.',
           f'In the CFD, 100 quads cut breathing-zone PM by {dq:.0f} % in 30 minutes against {d60:.0f} % for 100 heavy drones.', s)
s = re.sub(r'<td>−[0-9.]+ % in 30 min <small>\(100 quads; too weak to be the main mixer\)</small></td>',
           f'<td>−{dq:.0f} % in 30 min <small>(100 quads; too weak to be the main mixer)</small></td>', s)
s = s.replace('Standard k-ε does not damp turbulence in stable air, so these figures are an upper bound.',
              'Solid lines use a buoyant k-ε model calibrated against measured stable-air profiles (best estimate); '
              'the dashed line is the standard k-ε model, which does not damp turbulence in stable air (comparison). '
              'Checks of mesh, jet penetration and turbulence model: sim/cfd/REPORT.md, step 1.')
s = s.replace("their jets need about 15 minutes to work down to street level",
              "their jets need about 20 minutes to work down to street level")
s = s.replace("URANS with Boussinesq buoyancy and standard k-ε in OpenFOAM",
              "URANS with Boussinesq buoyancy and two k-ε turbulence models in OpenFOAM")
open(PAGE, "w").write(s)
print(f"best 60 m {d60:.1f} %, upper {d60u:.1f} %, quad {dq:.1f} %, 90 m {100 - best['heavy_z90'][-1]:.1f} %")
