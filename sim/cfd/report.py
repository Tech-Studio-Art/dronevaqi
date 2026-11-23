#!/usr/bin/env python3
"""Make figures + Markdown tables from all simulation outputs."""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# validated reference palette (dataviz skill): categorical order, blue sequential ramp
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
BLUES = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("seq", BLUES)
SEQ_R = LinearSegmentedColormap.from_list("seq_r", BLUES[::-1])
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "lines.linewidth": 2, "legend.frameon": False})


def load_json(p):
    return json.load(open(p)) if os.path.exists(p) else None


def slice_plane(vtu, normal, origin, field):
    import pyvista as pv
    g = pv.read(vtu)
    s = g.slice(normal=normal, origin=origin)
    s = s.cell_data_to_point_data()
    return s.points, np.asarray(s.point_data[field])


# ------------------------------------------------------------------ drones
def drones():
    rows = []
    for key in ("quad", "heavy"):
        s = 5.0 if key == "heavy" else 1.0
        arm = 0.8 if key == "heavy" else 0.225
        data = {}
        for fpos in ("none", "above", "below"):
            case = os.path.join(RUNS, "drone3d", f"{key}_{fpos}")
            r = load_json(os.path.join(case, "result.json"))
            if not r:
                continue
            rows.append(r)
            vtu = sorted(glob.glob(os.path.join(case, "VTK_full", "*", "internal.vtu")))[-1]
            pts, U = slice_plane(vtu, "y", (0, arm, 0), "U")
            pm = None
            v2 = sorted(glob.glob(os.path.join(case, "VTK_pm25", "*", "internal.vtu")))
            if v2:
                _, pm = slice_plane(v2[-1], "y", (0, arm, 0), "T")
            m = (np.abs(pts[:, 0]) < 1.5 * s) & (pts[:, 2] > -4 * s) & (pts[:, 2] < 1.2 * s)
            data[fpos] = (pts[m], np.linalg.norm(U, axis=1)[m], None if pm is None else pm[m], r)
        if not data:
            continue
        umax = max(v[1].max() for v in data.values())
        fig, axs = plt.subplots(2, 3, figsize=(10, 6.4), constrained_layout=True)
        for c, fpos in enumerate(("none", "above", "below")):
            if fpos not in data:
                continue
            P, um, pm, r = data[fpos]
            a = axs[0, c].tricontourf(P[:, 0], P[:, 2], um, levels=np.linspace(0, umax, 13), cmap=SEQ)
            if pm is not None:
                b = axs[1, c].tricontourf(P[:, 0], P[:, 2], 100 * (1 - pm), levels=np.linspace(0, 12, 13), cmap=SEQ)
            else:
                axs[1, c].text(0.5, 0.5, "no filter:\nPM unchanged", ha="center", va="center",
                               transform=axs[1, c].transAxes, color=INK2)
            title = {"none": "open rotors (no filter)", "above": "ducted filter ABOVE rotors",
                     "below": "ducted filter BELOW rotors"}[fpos]
            axs[0, c].set_title(title, color=INK)
            for ax in axs[:, c]:
                ax.set_aspect("equal")
                ax.set_xlabel("x (m)")
                ax.set_ylabel("z (m)")
                ax.grid(False)
        fig.colorbar(a, ax=axs[0, :], label="|U| (m/s)", shrink=0.8)
        fig.colorbar(b, ax=axs[1, :], label="PM2.5 reduction vs ambient (%)", shrink=0.8)
        fig.suptitle(data["none"][3]["drone"] + " - plane through one rotor row", color=INK)
        fig.savefig(os.path.join(FIG, f"drone3d_{key}.png"), dpi=150)
        plt.close(fig)
    return rows


# ------------------------------------------------------------------ mixing
def mixing():
    res = [load_json(p) for p in sorted(glob.glob(os.path.join(RUNS, "mixing", "*", "result.json")))]
    res = [r for r in res if r]
    if not res:
        return []
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
    order = ["no_drone", "quad_z30", "heavy_z60", "heavy_z90"]
    labels = {"no_drone": "no drones", "quad_z30": "100 quads @ 30 m", "heavy_z60": "100 heavy @ 60 m",
              "heavy_z90": "100 heavy @ 90 m"}
    for i, name in enumerate(order):
        r = next((x for x in res if x["case"] == name), None)
        if not r:
            continue
        t = np.arange(len(r["series"])) * 300 / 60
        for ax, key in zip(axs, ("z<2m", "z<10m")):
            ax.plot(t, [100 * s[key] for s in r["series"]], color=CAT[i], label=labels[name], marker="o", ms=4)
    for ax, key in zip(axs, ("0-2 m (breathing zone)", "0-10 m")):
        ax.set_title(f"PM in {key}, unit-cell average", color=INK)
        ax.set_xlabel("minutes after drones start")
        ax.set_ylabel("% of initial surface PM")
    axs[1].legend(loc="lower left")
    fig.savefig(os.path.join(FIG, "mixing_timeseries.png"), dpi=150)
    plt.close(fig)
    # final field heavy_z60
    case = os.path.join(RUNS, "mixing", "heavy_z60")
    vt = sorted(glob.glob(os.path.join(case, "VTK_all", "*", "internal.vtu")),
                key=lambda f: int(f.split("_")[-1].split("/")[0]) if f.split("_")[-1].split("/")[0].isdigit() else 0)
    if vt:
        import pyvista as pv
        fig, axs = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True, sharey=True)
        import mixing_axi as mx
        for ax, f, ttl in zip(axs, (vt[0], vt[-1]), ("t = 0", "t = 30 min")):
            g = pv.read(f)
            cc = np.asarray(g.cell_centers().points)
            r_ = np.hypot(cc[:, 0], cc[:, 1])
            C = mx.pm_from_theta(np.asarray(g.cell_data["T"]))
            m = cc[:, 2] < 140
            a = ax.tricontourf(r_[m], cc[m, 2], C[m], levels=np.linspace(0.45, 1.0, 12), cmap=SEQ)
            ax.set_title(f"PM / surface PM, {ttl}", color=INK)
            ax.set_xlabel("r (m) from drone axis")
            ax.grid(False)
        axs[0].set_ylabel("z (m)")
        fig.colorbar(a, ax=axs, label="PM / initial surface PM")
        fig.savefig(os.path.join(FIG, "mixing_field_heavy_z60.png"), dpi=150)
        plt.close(fig)
    return res


# ------------------------------------------------------------------ 2-D
def two_d():
    out = {}
    for p in glob.glob(os.path.join(RUNS, "2d", "summary_*.json")):
        for r in json.load(open(p)):
            out[r["case"]] = r
    sh = [r for r in out.values() if r["case"].startswith("shelter") and "error" not in r]
    if sh:
        fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True, sharey=True)
        for ax, wf in zip(axs, ("front", "back")):
            for i, U in enumerate((0.5, 1.0, 2.0)):
                pts = sorted((r["w_supply"], r["reduction_pct"]) for r in sh if r["wind_from"] == wf and r["Uref"] == U)
                if pts:
                    ax.plot(*zip(*pts), color=CAT[i], marker="o", ms=5, label=f"wind {U} m/s @10 m")
            ax.set_title(f"wind {'into opening' if wf == 'front' else 'over back wall'}", color=INK)
            ax.set_xlabel("filtered supply velocity from roof (m/s)")
        axs[0].set_ylabel("breathing-zone PM reduction (%)")
        axs[1].legend()
        fig.savefig(os.path.join(FIG, "shelter.png"), dpi=150)
        plt.close(fig)
    es = [r for r in out.values() if r["case"].startswith("espfence") and "error" not in r]
    if es:
        fig, ax = plt.subplots(figsize=(6, 3.4), constrained_layout=True)
        i = 0
        for U in (1.0, 2.0):
            for size, ls in (("PM25", "-"), ("coarse", "--")):
                r = next((x for x in es if x["Uref"] == U and x["size"] == size), None)
                if not r:
                    continue
                xs = [float(k[2:-1]) for k in r["pedestrian_C"]]
                ax.plot(xs, [100 * (1 - v) for v in r["pedestrian_C"].values()], color=CAT[0 if U == 1.0 else 1],
                        ls=ls, marker="o", ms=4, label=f"{'PM2.5' if size == 'PM25' else 'coarse PM'}, {U} m/s")
        ax.set_xscale("log")
        ax.set_xlabel("distance downwind of 6 m ESP fence (m)")
        ax.set_ylabel("PM reduction at 1-2 m height (%)")
        ax.legend()
        fig.savefig(os.path.join(FIG, "esp_fence.png"), dpi=150)
        plt.close(fig)
    return out


if __name__ == "__main__":
    import sys
    sys.path.insert(0, HERE)
    d = drones()
    m = mixing()
    t = two_d()
    json.dump(dict(drones=d, mixing=m, twod=t), open(os.path.join(RUNS, "all_results.json"), "w"),
              indent=1, default=float)
    print("figures:", sorted(os.listdir(FIG)))
