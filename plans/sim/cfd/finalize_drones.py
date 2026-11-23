#!/usr/bin/env python3
"""Final CADR bookkeeping for drone3d filter cases.

PM2.5: CADR = lambda * integral(C dV) over the filter (resolved tracer sink).
Coarse (2.5-10 um): single-pass capture ~98 %, sink lambda ~ 60-320 1/s. OpenFOAM's
semiImplicitSource applies a negative (sink) coefficient explicitly, which is
unstable at such stiffness. The PM2.5 tracer shows that the air entering the filter
is ambient (inlet concentration recovered from the plug-flow average
C_avg = C_in (1-P)/(-ln P) is ~1.00, i.e. no re-ingestion of cleaned air), so the
coarse removal is closed by mass balance: CADR_c = Q_filter * eta_c * C_in,
with Q_filter from the CFD through-flow.
"""
import json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import drone3d as d

out = []
for key in ("quad", "heavy"):
    for fpos in ("none", "above", "below"):
        p = os.path.join(d.RUNS, f"{key}_{fpos}", "result.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        if fpos != "none":
            h25 = next(h for h in r["history"] if h["pass_"] == "PM25")
            hf = next(h for h in r["history"] if h["pass_"] == "coarse")   # final flow state
            P25 = 1 - h25["eff"]
            C_in = h25["T_avg"] * (-math.log(P25)) / (1 - P25)
            A_f = hf["filt_vol"] / d.geometry(d.DRONES[key])["Ls"]
            Q_f = hf["uz_f"] * A_f
            cc = Q_f * hf["eff"] * C_in
            c25 = r["CADR_PM25"]
            c10 = (d.C_PM25 * c25 + (d.C_PM10 - d.C_PM25) * cc) / d.C_PM10
            # finish the thrust iteration analytically (filter load at final velocity)
            F = d.MU_AIR * hf["uz_f"] * d.D_MEDIA * d.MEDIA.L / d.PLEAT * A_f
            r.update(C_in_filter=C_in, Q_filter=Q_f, CADR_coarse=cc, CADR_PM10=c10,
                     g_per_h_PM10=c10 * d.C_PM10 * 3.6e-3, F_filter_N=F,
                     thrust_converged_N=r["weight_N"] + F,
                     power_W_converged=r["power_W"] * ((r["weight_N"] + F) / r["thrust_N"]) ** 1.5,
                     )
            r["endurance_min_converged"] = d.DRONES[key]["batt_wh"] / r["power_W_converged"] * 60
            json.dump(r, open(p, "w"), indent=1, default=float)
        # consistent power: momentum theory at the CFD-determined total disk thrust
        # (weight + aerodynamic filter load), same FM / coaxial / drive factors for all
        cfg = d.DRONES[key]
        A = 4 * math.pi * (cfg["D"] / 2) ** 2
        Tt = r.get("thrust_converged_N", r["weight_N"])
        P = Tt ** 1.5 / math.sqrt(2 * d.RHO_AIR * A) / cfg["fm"] * cfg["coax"] / cfg["eta"]
        r.update(power_mt_W=P, endurance_mt_min=cfg["batt_wh"] / P * 60, disk_area=A)
        json.dump(r, open(p, "w"), indent=1, default=float)
        out.append({k: v for k, v in r.items() if k != "history"})
json.dump(out, open(os.path.join(d.RUNS, "summary.json"), "w"), indent=1, default=float)
for r in out:
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
