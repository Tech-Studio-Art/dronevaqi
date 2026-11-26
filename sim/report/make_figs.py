#!/usr/bin/env python3
"""Figures for the pilot-proposal report (Research Dex assets)."""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle, FancyArrowPatch  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
CFD = os.path.join(SIM, "cfd")
sys.path.insert(0, CFD)
OUT = os.path.join(HERE, "assets")
os.makedirs(OUT, exist_ok=True)

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
BLUES = ["#ffffff", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("seq", BLUES)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "axes.facecolor": "white", "lines.linewidth": 2,
                     "legend.frameon": False, "savefig.dpi": 160})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- inversion mixing
def mixing():
    m = json.load(open(os.path.join(CFD, "runs", "mixing", "summary.json")))
    p2 = os.path.join(CFD, "runs", "mixing", "heavy_z60_prt2", "result.json")
    if os.path.exists(p2):
        m["heavy_z60_prt2"] = json.load(open(p2))
    lab = [("no_drone", "No drones"), ("quad_z30", "100 quad drones at 30 m"),
           ("heavy_z60", "100 heavy drones at 60 m"), ("heavy_z90", "100 heavy drones at 90 m"),
           ("heavy_z60_prt2", "100 heavy drones at 60 m, conservative (Pr_t = 2)")]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    for i, (k, l) in enumerate(lab):
        if k not in m:
            continue
        s = m[k]["series"]
        y = [100 * x["z<2m"] / s[0]["z<2m"] for x in s]
        t = np.arange(len(y)) * 5
        ax.plot(t, y, color=["#8a8f9c", CAT[1], CAT[0], CAT[2], CAT[5]][i], marker="o", ms=4,
                ls="--" if "prt2" in k else "-", label=f"{l}: {y[-1]:.0f} %")
    ax.set_xlabel("Time after drones start (min)")
    ax.set_ylabel("PM at 0–2 m, % of start value")
    ax.set_ylim(75, 101)
    ax.legend(loc="lower left", fontsize=8.5)
    save(fig, "aqipilot_fig2_mixing_time.png")

    import pyvista as pv
    import mixing_axi as mx
    vt = sorted(glob.glob(os.path.join(CFD, "runs", "mixing", "heavy_z60", "VTK_all", "*", "internal.vtu")),
                key=lambda f: int(f.split("_")[-1].split("/")[0]))
    if vt:
        fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4), sharey=True, constrained_layout=True)
        for ax, f, ttl in zip(axs, (vt[0], vt[-1]), ("Before (t = 0)", "After 30 min of mixing")):
            g = pv.read(f)
            cc = np.asarray(g.cell_centers().points)
            r = np.hypot(cc[:, 0], cc[:, 1])
            C = mx.pm_from_theta(np.asarray(g.cell_data["T"]))
            k = cc[:, 2] < 140
            a = ax.tricontourf(r[k], cc[k, 2], 100 * C[k], levels=np.linspace(45, 100, 12), cmap=SEQ)
            ax.set_title(ttl, fontsize=10)
            ax.set_xlabel("Distance from drone axis (m)")
            ax.grid(False)
            ax.plot([0], [60], marker="v", color=CAT[1], ms=9)
        axs[0].set_ylabel("Height (m)")
        fig.colorbar(a, ax=axs, label="PM, % of start street-level value")
        save(fig, "aqipilot_fig3_mixing_field.png")


# ---------------------------------------------------------------- ESP fence
def esp_schematic():
    fig, ax = plt.subplots(figsize=(7.4, 3.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.6, 4.2)
    ax.axis("off")
    ax.add_patch(Rectangle((2.6, 0.4), 0.15, 3.2, color="#9aa3b2"))
    ax.text(2.67, 0.05, "Grounded\nouter mesh", ha="center", va="top", fontsize=8, color=INK2)
    for y in np.linspace(0.6, 3.4, 8):
        ax.plot([3.25, 3.25], [y - 0.08, y + 0.08], color=CAT[1], lw=3)
    ax.text(3.25, 3.75, "Ionizer\n8–12 kV", ha="center", fontsize=8, color=CAT[1])
    for i, y in enumerate(np.linspace(0.5, 3.5, 13)):
        ax.plot([3.8, 5.8], [y, y], color=CAT[0] if i % 2 else "#4a4f5c", lw=2)
    ax.text(4.8, 3.75, "Collector plates\n(alternate HV / ground)", ha="center", fontsize=8, color=CAT[0])
    ax.add_patch(Rectangle((6.05, 0.4), 0.15, 3.2, color="#9aa3b2"))
    ax.text(6.12, 0.05, "Rear\nmesh", ha="center", va="top", fontsize=8, color=INK2)
    for y in (1.0, 2.0, 3.0):
        ax.add_patch(FancyArrowPatch((0.3, y), (2.4, y), arrowstyle="-|>", mutation_scale=14, color="#8a8f9c", lw=1.6))
        ax.add_patch(FancyArrowPatch((6.4, y), (8.8, y), arrowstyle="-|>", mutation_scale=14, color=CAT[2], lw=1.6))
    ax.text(1.35, -0.25, "Wind with PM\n(no fan)", ha="center", fontsize=8.5, color=INK2)
    ax.text(7.6, -0.25, "Cleaner air", ha="center", fontsize=8.5, color=CAT[2])
    rng = np.random.default_rng(1)
    ax.scatter(rng.uniform(0.4, 2.3, 40), rng.uniform(0.6, 3.4, 40), s=6, color="#5b6070")
    ax.scatter(rng.uniform(6.5, 8.6, 12), rng.uniform(0.6, 3.4, 12), s=6, color="#5b6070")
    ax.text(9.2, 2.0, "Solar panel +\nHV supply\n140 W per m", fontsize=8, color=INK2, va="center")
    save(fig, "aqipilot_fig4_esp_principle.png")


def esp_results():
    road = json.load(open(os.path.join(CFD, "runs", "esp_road", "summary.json")))["total"]
    amb = {}
    for r in json.load(open(os.path.join(CFD, "runs", "2d", "summary_all.json"))):
        if r["case"].startswith("espfence"):
            amb[r["case"]] = r
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.4), constrained_layout=True)
    cols = {"no fence": "#8a8f9c", "passive fence": CAT[3], "ESP fence": CAT[0], "ESP fence, half migration velocity": CAT[5]}
    for k, v in road.items():
        x = [float(a) for a in v]
        axs[0].plot(x, [300 * c for c in v.values()], marker="o", ms=4, color=cols[k], label=k,
                    ls="--" if "half" in k else "-")
    axs[0].axhline(300, color=INK2, lw=0.8, ls=":")
    axs[0].text(2.6, 307, "background 300", fontsize=7.5, color=INK2)
    axs[0].set_xscale("log")
    axs[0].set_xticks([2.5, 5, 10, 20, 40])
    axs[0].set_xticklabels(["2.5", "5", "10", "20", "40"])
    axs[0].minorticks_off()
    axs[0].set_xlabel("Distance behind fence (m)")
    axs[0].set_ylabel("PM2.5 at 1–2 m (µg/m³)")
    axs[0].set_title("(a) Fence between road and footpath, wind 1 m/s", fontsize=9)
    axs[0].legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    for i, (U, size, ls) in enumerate(((1.0, "PM25", "-"), (1.0, "coarse", "--"), (2.0, "PM25", "-"), (2.0, "coarse", "--"))):
        r = amb.get(f"espfence_U{U}_{size}")
        if not r:
            continue
        xs = [float(k[2:-1]) for k in r["pedestrian_C"]]
        axs[1].plot(xs, [100 * (1 - v) for v in r["pedestrian_C"].values()], color=CAT[0] if U == 1 else CAT[1], ls=ls,
                    marker="o", ms=4, label=f"{'PM2.5' if size == 'PM25' else 'PM10 coarse part'}, {U:g} m/s")
    axs[1].set_xscale("log")
    axs[1].set_xticks([2, 5, 10, 20, 40, 80])
    axs[1].set_xticklabels(["2", "5", "10", "20", "40", "80"])
    axs[1].minorticks_off()
    axs[1].set_xlabel("Distance behind fence (m)")
    axs[1].set_ylabel("Decrease at 1–2 m height (%)")
    axs[1].set_title("(b) Background air only", fontsize=9)
    axs[1].legend(fontsize=7.5)
    save(fig, "aqipilot_fig5_esp_results.png")


# ---------------------------------------------------------------- surveillance
def surveillance():
    S = os.path.join(SIM, "surveillance")
    r = json.load(open(os.path.join(S, "results.json")))
    X, Y, C = np.load(os.path.join(S, "plume_night.npy"))
    fig, ax = plt.subplots(figsize=(7.4, 3.3), constrained_layout=True)
    from matplotlib.colors import BoundaryNorm
    lv = [5, 10, 25, 50, 100, 250, 500, 1000, 5000]
    cm = LinearSegmentedColormap.from_list("p", ["#e8f1fc", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
    C = np.where(C < 5, np.nan, C)
    a = ax.contourf(X, Y, C, levels=lv, cmap=cm, norm=BoundaryNorm(lv, cm.N))
    ax.contour(X, Y, C, levels=[25], colors=[CAT[1]], linewidths=1.5)
    g = np.arange(-500, 2001, 500)
    gx, gy = np.meshgrid(g, np.arange(-500, 501, 500))
    ax.scatter(gx, gy, marker="s", s=22, color=INK, zorder=3, label="Sensor, 500 m grid")
    ax.plot([0], [0], marker="*", ms=14, color=CAT[1], zorder=4, label="Waste fire, 100 kg/h")
    ax.set_xlim(-150, 1600)
    ax.set_ylim(-550, 550)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("Downwind distance (m), night, wind 1 m/s")
    ax.set_ylabel("Crosswind (m)")
    ax.legend(loc="upper right", fontsize=8)
    cb = fig.colorbar(a, ax=ax, shrink=0.9)
    cb.set_label("PM2.5 added (µg/m³)")
    cb.set_ticks(lv[:-1])
    save(fig, "aqipilot_fig6_fire_plume.png")

    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.2), constrained_layout=True)
    for i, (cls, lab) in enumerate((("night (stable, E-F)", "Night, stable air"), ("day (neutral, D)", "Day, neutral air"))):
        d = r["detection"][cls]
        axs[0].plot([x["spacing_m"] for x in d], [100 * x["p_detect"] for x in d], marker="o", ms=4, color=CAT[i], label=lab)
    axs[0].set_xscale("log")
    axs[0].set_xticks([250, 500, 1000, 2000, 6000])
    axs[0].set_xticklabels(["250", "500", "1000", "2000", "6000"])
    axs[0].minorticks_off()
    axs[0].set_xlabel("Sensor spacing (m)")
    axs[0].set_ylabel("Fires found by sensors (%)")
    axs[0].set_title("(a) Fixed sensors only", fontsize=9)
    axs[0].legend(fontsize=8)
    p = r["patrol"]
    lbl = ["Sensors\nonly", "+1 drone", "+2 drones", "+4 drones"]
    vals = [x["emission_avoided_pct"] for x in p]
    axs[1].bar(lbl, vals, color=[CAT[3]] + [CAT[0]] * 3, width=0.6)
    for i, v in enumerate(vals):
        axs[1].text(i, v + 1.5, f"{v:.0f} %", ha="center", fontsize=8.5)
    axs[1].set_ylim(0, 100)
    axs[1].set_ylabel("PM2.5 emission prevented per fire (%)")
    axs[1].set_title("(b) 25 km² district, 500 m sensor grid", fontsize=9)
    save(fig, "aqipilot_fig7_surveillance.png")


def design_sheet():
    """Screenshot the drone design sheets from docs/inversion-drones.html."""
    import subprocess
    page = os.path.join(os.path.dirname(SIM), "docs", "inversion-drones.html")
    js = f"""
const {{ chromium }} = require('playwright');
(async () => {{
  const b = await chromium.launch();
  const p = await b.newPage({{ viewport: {{ width: 1180, height: 900 }}, deviceScaleFactor: 1.6, colorScheme: 'light' }});
  await p.goto('file://{page}');
  await p.waitForTimeout(800);
  const secs = await p.$$('section.sheet');
  const sv = await secs[0].$('svg'); await sv.screenshot({{ path: '{OUT}/aqipilot_fig1_mission.png' }});
  const d1 = await secs[1].$('.drawings'); await d1.screenshot({{ path: '{OUT}/aqipilot_fig1b_im50.png' }});
  await b.close();
}})();
"""
    f = os.path.join(HERE, "_shot.js")
    open(f, "w").write(js)
    subprocess.run(["node", f], check=True)
    os.remove(f)


if __name__ == "__main__":
    design_sheet()
    mixing()
    esp_schematic()
    esp_results()
    surveillance()
    print(sorted(os.listdir(OUT)))
