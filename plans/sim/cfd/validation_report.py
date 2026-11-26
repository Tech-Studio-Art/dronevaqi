#!/usr/bin/env python3
"""Collect step-1b validation results -> runs/validation/summary.json + figures/validation_*.png"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
V = os.path.join(HERE, "runs", "validation")
FIG = os.path.join(HERE, "figures")
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
INK2, GRID = "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9.5, "axes.grid": True, "grid.color": GRID, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False, "figure.facecolor": "white",
                     "lines.linewidth": 2, "savefig.dpi": 160})


def main():
    out = {}
    # (i) jet penetration
    jp = {n: json.load(open(os.path.join(V, "jetpen", n, "result.json"))) for n in ("heavy", "quad")}
    out["jet_penetration"] = {n: {k: v for k, v in r.items() if k != "depth_vs_time"} for n, r in jp.items()}
    # (ii) stable surface layer
    mf = json.load(open(os.path.join(V, "most_fit.json")))
    out["surface_layer_fits"] = {k: dict(intercept=v["intercept"], slope=v["slope"], n=v["n"]) for k, v in mf.items()}
    if "phih_slope" in mf.get("Cg=1.5", {}):
        out["surface_layer_fits"]["Cg=1.5"].update(phih_intercept=mf["Cg=1.5"]["phih_intercept"],
                                                   phih_slope=mf["Cg=1.5"]["phih_slope"])
    # (iii) fence
    out["fence"] = json.load(open(os.path.join(V, "fence", "summary.json")))
    json.dump(out, open(os.path.join(V, "summary.json"), "w"), indent=1)

    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), constrained_layout=True)
    # (a) surface layer
    ax = axs[0]
    for i, (lab, key) in enumerate((("standard k-ε", "std"), ("buoyant k-ε, Cg = 1", "Cg=1"), ("buoyant k-ε, Cg = 1.5", "Cg=1.5"))):
        zs, ps = [], []
        for r in mf[key]["recs"]:
            if r["hours_cooling"] < 3:
                continue
            for z, zt, p in zip(r["z"], r["zeta"], r["phim"]):
                if 20 < z < 80 and 0.05 < zt < 1.0:
                    zs.append(zt)
                    ps.append(p)
        ax.scatter(zs, ps, s=8, color=CAT[i], label=lab, alpha=0.8)
    zz = np.linspace(0, 1, 20)
    ax.plot(zz, 1 + 4.7 * zz, color="#333", lw=1.2, ls="--", label="Businger et al. 1971")
    ax.plot(zz, 1 + 5 * zz, color="#333", lw=1.2, ls=":", label="Dyer 1974")
    ax.plot(zz, 1 + 6 * zz, color="#888", lw=1.2, ls="-.", label="Högström 1988")
    ax.set_xlabel("stability ζ = z/Λ")
    ax.set_ylabel("φm = κz/u* dU/dz")
    ax.set_title("(a) Stable surface layer", fontsize=9.5)
    ax.legend(fontsize=7)
    # (b) jet penetration
    ax = axs[1]
    for i, n in enumerate(("heavy", "quad")):
        r = jp[n]
        t, d = zip(*r["depth_vs_time"])
        ax.plot(np.array(t) / 60, d, color=CAT[i], marker="o", ms=3, label=f"CFD, {n} drone")
        lo, hi = r["mtt_depth_by_alpha"]["0.1"], r["mtt_depth_by_alpha"]["0.05"]
        ax.axhspan(lo, hi, color=CAT[i], alpha=0.12)
        ax.axhline(r["mtt_depth_by_alpha"]["0.076"], color=CAT[i], lw=1, ls="--")
    ax.set_xlabel("time (min)")
    ax.set_ylabel("jet penetration depth (m)")
    ax.set_title("(b) Jet into stable air vs MTT model", fontsize=9.5)
    ax.legend(fontsize=7.5, loc="center right")
    # (c) fence through-flow
    ax = axs[2]
    sc = out["fence"]["screen"]
    K = np.linspace(0, 6.5, 50)
    ax.plot(K, 1 / (1 + K / 4), color="#333", lw=1.2, ls="--", label="screen theory 1/(1+K/4)")
    ax.plot([r["K"] for r in sc], [r["cfd"] for r in sc], "o", color=CAT[0], label="CFD, 6 m fence")
    ax.axvline(1.5, color=CAT[1], lw=1, ls=":")
    ax.text(1.6, 0.95, "design K = 1.5", fontsize=7.5, color=CAT[1])
    ax.set_xlabel("pressure-loss coefficient K")
    ax.set_ylabel("through-flow / free wind")
    ax.set_title("(c) Air flow through the fence", fontsize=9.5)
    ax.legend(fontsize=7.5)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "validation_step1b.png"), bbox_inches="tight")
    print(json.dumps({k: v for k, v in out.items()}, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
