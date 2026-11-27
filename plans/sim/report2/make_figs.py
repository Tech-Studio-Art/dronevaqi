#!/usr/bin/env python3
"""Figures for report 2 (drone fleet in wind). Data: plans/sim/cfd/runs/fleet3d*/ (summary.json, cells.json,
calm_rhythm.json)."""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
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
DIV = LinearSegmentedColormap.from_list("div", ["#b3471c", "#eb6834", "#f6c8b0", "#f7f7f5", "#cde2fb", "#3987e5", "#104281"])
INK, INK2, GRID, MUTED = "#0b0b0b", "#52514e", "#e4e3df", "#8a8f9c"
plt.rcParams.update({"font.size": 9.5, "font.family": "DejaVu Sans", "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "axes.facecolor": "white", "lines.linewidth": 1.8,
                     "legend.frameon": False, "savefig.dpi": 160})
X0, X1, Y1 = -296.875, 1003.125, 403.125
RHO, G, TH0 = 1.25, 9.81, 285.0
N_BV = math.sqrt(G / TH0 * 5.0 / 100.0)             # buoyancy frequency of the inversion
FM1 = 50 * 9.81 / (math.pi * RHO)                    # Briggs momentum flux of one drone (m4/s2)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight")
    plt.close(fig)


def load(d):
    p = os.path.join(CFD, "runs", d, "summary.json")
    return json.load(open(p)) if os.path.exists(p) else None


def J(*p):
    return json.load(open(os.path.join(CFD, "runs", *p)))


def full(m):
    m = np.asarray(m)
    return np.vstack([m[::-1], m])


# ---------------------------------------------------------------- 1 layout
def fig1_layout():
    import fleet3d as f
    fig = plt.figure(figsize=(7.6, 4.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.95, 1], hspace=0.42, wspace=0.08)
    titles = ["(a) Dispersed fleet: 100 sites, 50 m grid", "(b) K4: 4 clusters of 25 drones", "(c) K1: 1 cluster of 100 drones"]
    for n in range(3):
        a = fig.add_subplot(gs[0, n])
        a.add_patch(Rectangle((X0, -Y1), X1 - X0, 2 * Y1, fc="#f4f7fb", ec=INK2, lw=0.7))
        a.add_patch(Rectangle((0, -250), 500, 500, fc="none", ec=INK, lw=1.0))
        a.axhline(0, color=INK2, lw=0.6, ls="-.")
        if n == 0:
            for x, y, g in f.fleet_sites():
                for yy in (y, -y):
                    a.plot(x, yy, "o", ms=2.6, mfc=CAT[0] if g == "A" else "white", mec=CAT[0], mew=0.8)
        else:
            boxes = f.CLUSTERS["K4" if n == 1 else "K1"]["boxes"]
            for (bx0, by0), (bx1, by1) in boxes:
                lo = -by1 if by0 > 0 else -by1
                a.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fc=CAT[1], ec=CAT[1]))
                a.add_patch(Rectangle((bx0, -by1), bx1 - bx0, by1 - by0 if by0 > 0 else by1, fc=CAT[1], ec=CAT[1]))
        a.annotate("", xy=(-100, 330), xytext=(-280, 330), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
        a.text(-280, 230, "Wind", fontsize=7)
        a.set_xlim(X0 - 10, X1 + 10)
        a.set_ylim(-Y1 - 10, Y1 + 10)
        a.set_aspect("equal")
        a.set_title(titles[n], fontsize=7.8)
        a.set_xticks([0, 500, 1000])
        a.tick_params(labelsize=7)
        if n:
            a.set_yticklabels([])
        else:
            a.set_ylabel("y (m)", fontsize=8)
        a.set_xlabel("x (m)", fontsize=8)
        a.grid(False)
    b = fig.add_subplot(gs[1, :])
    z = np.linspace(0, 200, 201)
    th = 285 + 5 * np.minimum(z, 100) / 100
    b.imshow(np.tile(th[:, None], (1, 2)), extent=[X0, X1, 0, 200], origin="lower", aspect="auto",
             cmap=LinearSegmentedColormap.from_list("t", ["#dbe8f8", "#fbe3d4"]), alpha=0.9)
    b.axhline(100, color=INK2, lw=0.8, ls="--")
    b.text(X1 - 10, 104, "Top of the inversion (100 m)", ha="right", fontsize=7.5, color=INK2)
    for x in range(25, 500, 50):
        b.plot(x, 60, marker="s", ms=3.5, color=CAT[0])
    zz = np.linspace(0.5, 190, 60)
    u = 0.0888 / 0.41 * np.log((zz + 0.1) / 0.1)
    b.plot(X0 + 30 + 60 * u, zz, color=INK, lw=1.1)
    b.text(X0 + 30, 160, "Inlet wind,\nU10 = 1 m/s", fontsize=7)
    b.text(250, 70, "Drones hover at 60 m", ha="center", fontsize=7.5, color=CAT[0])
    b.text(760, 30, "PM: 100 % at street level,\n50 % at 100 m and above", ha="center", fontsize=7.2, color=INK2)
    b.set_xlim(X0 - 10, X1 + 10)
    b.set_ylim(0, 200)
    b.set_xlabel("x, along the wind (m)", fontsize=8)
    b.set_ylabel("Height (m)", fontsize=8)
    b.set_title("(d) Side view", fontsize=7.8)
    b.grid(False)
    save(fig, "aqifleet_fig1_layout.png")


# ---------------------------------------------------------------- 2 disk check
def fig2_disk():
    c = J("fleet3d", "cells.json")
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    rows = [("resolved_axi", INK, "-", "o"), ("cell", CAT[1], "-", "s"), ("cell_fine", CAT[1], "--", "s"),
            ("cell_r3p5", CAT[2], "--", "^"), ("cell_1col", CAT[0], "-", "o")]
    for k, col, ls, mk in rows:
        s = c[k]["series"]
        y = [100 * x["z<2m"] / s[0]["z<2m"] for x in s]
        lab = c[k]["label"][0].upper() + c[k]["label"][1:]
        ax.plot(np.arange(len(y)) * 5, y, color=col, ls=ls, marker=mk, ms=3.5,
                label=f"{lab}, {c[k]['disk_area_m2']:.0f} m²: {y[-1]:.0f} %")
    ax.set_xlabel("Time after drones start (min)")
    ax.set_ylabel("PM at street level,\n% of start value")
    ax.set_ylim(78, 101)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.2)
    save(fig, "aqifleet_fig2_disk_check.png")


# ---------------------------------------------------------------- 3 calm rhythm
def fig3_rhythm():
    r = J("fleet3d", "calm_rhythm.json")["cases"]
    fig, ax = plt.subplots(figsize=(7.0, 3.3))
    sty = {"pF100": (CAT[0], "-", "o"), "pR600": (CAT[2], "-", "s"), "pC50": (CAT[1], "-", "^"), "pR60": (CAT[3], "--", "v")}
    for k, (col, ls, mk) in sty.items():
        d = r[k]["decrease_pct"]
        lab = r[k]["label"][0].upper() + r[k]["label"][1:]
        ax.plot(np.arange(len(d)) * 5, d, color=col, ls=ls, marker=mk, ms=3.5,
                label=f"{lab}: {d[-1]:.0f} % at 30 min, {r[k]['mean_30min_pct']:.0f} % mean")
    ax.set_xlabel("Time after drones start (min)")
    ax.set_ylabel("Decrease of PM at\nstreet level (%)")
    ax.set_ylim(-0.5, 21)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, fontsize=7.2)
    save(fig, "aqifleet_fig3_rhythm.png")


# ---------------------------------------------------------------- 4 dispersed fleet in wind
def fig4_dispersed(S, S2):
    r = J("fleet3d", "calm_rhythm.json")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 3.6), gridspec_kw=dict(width_ratios=[1, 1.6]))
    zc = np.array(r["z"])
    a.plot(r["cases"]["pF100"]["profile_z_30min"], zc, color=MUTED, lw=1.6, label="Calm air, 100 drones, 30 min")
    z = np.array(S["z"])
    rows = [(S, "F100", CAT[0], "-", "Wind 1 m/s, 100 drones"), (S, "C50", CAT[1], "--", "Wind 1 m/s, 50 drones"),
            (S2, "F100", CAT[5], ":", "Wind 2 m/s, 100 drones")]
    for src, k, col, ls, lab in rows:
        if src and k in src["cases"]:
            c = src["cases"][k]
            a.plot(c["profile_z"], z, color=col, ls=ls, lw=1.6, label=f"{lab}, {c['minutes']:.0f} min")
    a.axvline(0, color=INK2, lw=0.6)
    a.axhline(60, color=CAT[0], lw=0.6, ls="--")
    a.set_ylim(0, 100)
    a.set_xlim(-25, 22)
    a.set_xlabel("Decrease of PM (%)\n(negative: increase)")
    a.set_ylabel("Height (m)")
    a.legend(fontsize=6.6, loc="lower left")
    a.set_title("(a) Zone-mean change of PM", fontsize=8.5)
    c = S["cases"]["F100"]
    x = np.array(S["x"])
    dec = 100 * (1 - np.array(c["section_pm"]) / np.array(c["section_ref"]))
    cs = b.pcolormesh(x, z, dec, cmap=DIV, norm=TwoSlopeNorm(0, -12, 12), shading="nearest", rasterized=True)
    u, w = np.array(c["section_u"]), np.array(c["section_w"])
    sx = slice(1, None, 5)
    zi = [k for k in range(len(z)) if z[k] < 95][::2]
    XX, ZZ = np.meshgrid(x[sx], z[zi])
    b.quiver(XX, ZZ, u[zi][:, sx], w[zi][:, sx], color=INK, scale=45, width=0.0016, alpha=0.7)
    for xx in range(25, 500, 50):
        b.plot(xx, 60, "s", ms=2.5, color=INK)
    b.set_xlim(-150, 1003)
    b.set_ylim(0, 100)
    b.set_xlabel("x, along the wind (m)")
    b.set_title(f"(b) 100 drones, wind 1 m/s, {c['minutes']:.0f} min, section at y = {c['section_y']:.0f} m", fontsize=8.5)
    b.grid(False)
    cb = fig.colorbar(cs, ax=b, shrink=0.9, pad=0.02)
    cb.set_label("Decrease of PM (%)", fontsize=8)
    fig.tight_layout()
    save(fig, "aqifleet_fig4_dispersed.png")


# ---------------------------------------------------------------- 5 penetration against Briggs
def penetration_points(S, S2):
    """(label, momentum flux Fm, wind at the drones U60, penetration below the drones, reached ground?)"""
    pts = []
    for src in (S, S2):
        if not src:
            continue
        u60 = src.get("U60") or 1.43
        for k, nd in (("F100", 1), ("C50", 1), ("K4", 25), ("K1", 100)):
            if k in src["cases"]:
                c = src["cases"][k]
                pen = 60 - c["mix_bottom_m"]
                pts.append((f"{k}, wind {src['U10']:g} m/s", nd * FM1, u60, pen, c["mix_bottom_m"] <= 1.0))
    return pts


def fig5_briggs(S, S2):
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    xx = np.linspace(0, 80, 50)
    ax.plot(xx, 1.5 * xx, color=INK, lw=1.2, label="Briggs: Δh = 1.5 (Fm / (U N))^(1/3)")
    ax.axhline(60, color=CAT[1], lw=0.8, ls="--")
    ax.text(2, 62, "Ground (60 m below the drones)", fontsize=7.5, color=CAT[1])
    mk = {"F100": "o", "C50": "s", "K4": "D", "K1": "^"}
    for lab, fm, u, pen, ground in penetration_points(S, S2):
        k = lab.split(",")[0]
        col = CAT[0] if "1 m/s" in lab else CAT[5]
        xv = (fm / (u * N_BV)) ** (1 / 3)
        ax.plot(xv, min(pen, 60), mk[k], ms=6, color=col, mfc=col if not ground else "white", mew=1.5, label=f"CFD {lab}")
    ax.set_xlabel("(Fm / (U N))^(1/3) (m)")
    ax.set_ylabel("Depth of mixing below\nthe drones (m)")
    ax.set_xlim(0, 75)
    ax.set_ylim(0, 110)
    ax.legend(fontsize=6.8, loc="upper left")
    save(fig, "aqifleet_fig5_briggs.png")


# ---------------------------------------------------------------- 6 cluster sections
def fig6_sections(S, S2=None):
    keys = [k for k in ("K1", "K4") if k in S["cases"]]
    x, z = np.array(S["x"]), np.array(S["z"])
    fig, axs = plt.subplots(len(keys), 1, figsize=(7.4, 2.5 * len(keys)), sharex=True, squeeze=False)
    for ax, k in zip(axs[:, 0], keys):
        c = S["cases"][k]
        dec = 100 * (1 - np.array(c["section_pm"]) / np.array(c["section_ref"]))
        cs = ax.pcolormesh(x, z, dec, cmap=DIV, norm=TwoSlopeNorm(0, -25, 25), shading="nearest", rasterized=True)
        u, w = np.array(c["section_u"]), np.array(c["section_w"])
        sx = slice(1, None, 5)
        zi = [n for n in range(len(z)) if z[n] < 95][::2]
        XX, ZZ = np.meshgrid(x[sx], z[zi])
        ax.quiver(XX, ZZ, u[zi][:, sx], w[zi][:, sx], color=INK, scale=45, width=0.0016, alpha=0.7)
        ax.set_ylim(0, 100)
        ax.set_xlim(-150, 1003)
        ax.set_ylabel("Height (m)")
        name = "K1: 1 cluster of 100 drones" if k == "K1" else "K4: 4 clusters of 25 drones"
        ax.set_title(f"{name}, wind 1 m/s, {c['minutes']:.0f} min, section at y = {c['section_y']:.1f} m", fontsize=8.5)
        ax.grid(False)
        cb = fig.colorbar(cs, ax=ax, shrink=0.9, pad=0.01)
        cb.set_label("Decrease of PM (%)", fontsize=7.5)
    axs[-1, 0].set_xlabel("x, along the wind (m)")
    save(fig, "aqifleet_fig6_sections.png")


# ---------------------------------------------------------------- 7 maps
def fig7_maps(S, S2=None):
    rows = [(S, "K1", "K1, wind 1 m/s"), (S, "K4", "K4, wind 1 m/s")]
    if S2 and "K1" in S2["cases"]:
        rows.append((S2, "K1", "K1, wind 2 m/s"))
    rows = [r for r in rows if r[1] in r[0]["cases"]]
    fig, axs = plt.subplots(len(rows), 1, figsize=(7.6, 2.75 * len(rows)), sharex=True, squeeze=False)
    lev = [-2, 2, 5, 10, 15, 20, 25, 30]
    norm = BoundaryNorm(lev, SEQ.N)
    for ax, (src, k, lab) in zip(axs[:, 0], rows):
        x, y = np.array(src["x"]), np.array(src["y"])
        yy = np.concatenate([-y[::-1], y])
        c = src["cases"][k]
        im = ax.pcolormesh(x, yy, full(c["map"]), cmap=SEQ, norm=norm, shading="nearest", rasterized=True)
        ax.add_patch(Rectangle((0, -250), 500, 500, fc="none", ec=INK, lw=0.8))
        ax.set_title(f"{lab}, {c['minutes']:.0f} min: zone mean decrease {c['zone_last']:.1f} %, largest {c['max_ground']:.0f} %", fontsize=8.3)
        ax.set_aspect("equal")
        ax.set_xlim(X0, X1)
        ax.set_ylim(-Y1, Y1)
        ax.set_ylabel("y (m)")
        ax.grid(False)
    axs[-1, 0].set_xlabel("x, along the wind (m)")
    cb = fig.colorbar(im, ax=axs[:, 0], shrink=0.85, pad=0.02)
    cb.set_label("Decrease of PM at street level (%)")
    save(fig, "aqifleet_fig7_maps.png")


# ---------------------------------------------------------------- 8 profiles and placement
def fig8_profiles(S, S2=None):
    rows = [(S, "K1", CAT[0], "-", "K1, wind 1 m/s"), (S, "K4", CAT[1], "-", "K4, wind 1 m/s"),
            (S, "F100", MUTED, "--", "Dispersed 100, wind 1 m/s")]
    if S2 and "K1" in S2["cases"]:
        rows.insert(2, (S2, "K1", CAT[5], "-.", "K1, wind 2 m/s"))
    rows = [r for r in rows if r[0] and r[1] in r[0]["cases"]]
    fig = plt.figure(figsize=(7.6, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.7, 1], hspace=0.45)
    a, b, c = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    for src, k, col, ls, lab in rows:
        d = src["cases"][k]
        a.plot(src["x"], d["profile_x"], color=col, ls=ls, lw=1.5, label=lab)
        w = d["window"]
        b.plot([float(v) for v in w], list(w.values()), color=col, ls=ls, lw=1.5)
        c.plot(src["y"], d["profile_y"], color=col, ls=ls, lw=1.5)
    a.axvspan(0, 500, color="#eef4fc", zorder=0)
    a.text(250, a.get_ylim()[1] * 0.92 if a.get_ylim()[1] > 0 else 1, "Fleet zone", ha="center", fontsize=7.5, color=CAT[0])
    a.axhline(0, color=INK2, lw=0.6)
    a.set_xlim(X0, X1)
    a.set_xlabel("x, along the wind (m)")
    a.set_ylabel("Decrease (%)")
    a.set_title("(a) Along the wind: mean over the zone width (y = −250 m to 250 m)", fontsize=8.5)
    a.legend(fontsize=7, loc="upper right")
    b.axvline(0, color=INK2, lw=0.6, ls="--")
    b.axhline(0, color=INK2, lw=0.6)
    b.set_xlabel("Upwind edge of the target area (m from\nthe upwind edge of the fleet zone)")
    b.set_ylabel("Mean decrease in a\n500 m × 500 m target (%)")
    b.set_title("(b) Position of the target area", fontsize=8.5)
    c.axvspan(0, 250, color="#eef4fc", zorder=0)
    c.axhline(0, color=INK2, lw=0.6)
    c.set_xlim(0, 403)
    c.set_xlabel("y, from the wind axis (m)")
    c.set_title("(c) Across the wind (x = 0–500 m)", fontsize=8.5)
    save(fig, "aqifleet_fig8_profiles.png")


# ---------------------------------------------------------------- 9 time series
def fig9_time(S, S2=None):
    cell = J("fleet3d", "cells.json")["cell_1col"]["series"]
    cdec = [100 * (1 - x["z<2m"] / cell[0]["z<2m"]) for x in cell]
    fig, ax = plt.subplots(figsize=(7.0, 3.3))
    ax.plot(np.arange(len(cdec)) * 5, cdec, color=MUTED, marker="o", ms=3, lw=1.3,
            label=f"Calm air, 100 drones (unit cell, no edges): {cdec[-1]:.0f} %")
    rows = [(S, "K1", CAT[0], "-", "K1, wind 1 m/s"), (S, "K4", CAT[1], "-", "K4, wind 1 m/s"),
            (S2, "K1", CAT[5], "-.", "K1, wind 2 m/s"), (S, "F100", INK2, "--", "Dispersed 100, wind 1 m/s")]
    for src, k, col, ls, lab in rows:
        if src and k in src["cases"]:
            d = src["cases"][k]
            ax.plot(d["t_min"], d["series"], color=col, ls=ls, marker="o", ms=3, lw=1.5, label=f"{lab}: {d['series'][-1]:.1f} %")
    ax.set_xlabel("Time after drones start (min)")
    ax.set_ylabel("Zone-mean decrease of\nPM at street level (%)")
    ax.set_ylim(bottom=-0.5)
    ax.legend(fontsize=7, loc="upper left")
    save(fig, "aqifleet_fig9_time.png")


if __name__ == "__main__":
    S, S2 = load(os.environ.get("FLEET_SUMMARY_DIR", "fleet3d")), load(os.environ.get("FLEET_SUMMARY_DIR2", "fleet3d_u2"))
    fig1_layout()
    fig2_disk()
    fig3_rhythm()
    if S:
        fig4_dispersed(S, S2)
        fig5_briggs(S, S2)
        if "K1" in S["cases"] or "K4" in S["cases"]:
            fig6_sections(S, S2)
            fig7_maps(S, S2)
            fig8_profiles(S, S2)
            fig9_time(S, S2)
    print(sorted(os.listdir(OUT)))
