#!/usr/bin/env python3
"""
3-D RANS CFD of hovering multirotors with/without rotor-stream filters.

Literature methods used (no new pollution model is invented):
  * Rotors: uniformly loaded actuator disk with prescribed thrust (Froude-Rankine
    momentum theory as a volumetric momentum source; standard for multirotor
    downwash CFD). Swirl neglected.
  * Turbulence: steady RANS, k-omega SST (Menter 1994).
  * Filter: porous medium, Darcy-Forchheimer source (OpenFOAM
    explicitPorositySource). Darcy coefficient from the Davies (1973) fibrous
    filter correlation, pleat area ratio 5.
  * Pollutant: Eulerian passive-scalar advection-diffusion, eddy diffusivity
    nu_t/Sc_t with Sc_t = 0.7 (Tominaga & Stathopoulos 2007). Solved as the
    scalar T of buoyantBoussinesqSimpleFoam with beta = 0 (no buoyancy), so T
    is a pure passive tracer with Pr_t = Sc_t = 0.7. Two size classes are run
    separately (fine PM2.5, coarse 2.5-10 um); PM10 = PM2.5 + coarse.
    Filter capture: first-order sink lambda = -ln(P) U/L in the filter zone,
    with single-pass penetration P from single-fibre filtration theory (Hinds
    1999) at the CFD face velocity. Settling (<= 5 mm/s at 10 um) is neglected
    against ~5 m/s rotor velocities.

Cases: {quad 2 kg, heavy 50 kg MTOW (25 kg payload)} x {no filter, filter
above rotors, filter below rotors}. Disk thrust is iterated so that net lift
(rotor thrust - aerodynamic load on the filter) equals weight.
"""
import json
import math
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from foam import write, run, field, write_field, vtk_zone, zone_stats  # noqa: E402
from drone_aqi_sim import FILTERS, RHO_AIR, MU_AIR, G  # noqa: E402

NU = MU_AIR / RHO_AIR
RUNS = os.path.join(HERE, "runs", "drone3d")
C_PM25, C_PM10 = 300.0, 500.0   # ug/m3, Delhi Nov episode: PM10 > PM2.5

DRONES = {
    "quad": dict(label="Quad 2 kg, 4 x 10 in rotors", mass=2.0, filt_mass=0.15,
                 D=0.254, arm=0.225, s=1.0, coax=1.0, fm=0.6, eta=0.85, batt_wh=95),
    "heavy": dict(label="Heavy-lift, 25 kg payload (50 kg MTOW), 4 coaxial 54 in rotor pairs",
                  mass=50.0, filt_mass=2.0, D=1.37, arm=0.80, s=5.0, coax=1.22, fm=0.6,
                  eta=0.85, batt_wh=2200),
}
PLEAT = 5.0
MEDIA = FILTERS["electret_coarse"]
D_MEDIA = 64 * MEDIA.alpha ** 1.5 * (1 + 56 * MEDIA.alpha ** 3) / MEDIA.d_f ** 2


def coarse_eff(U):
    d = np.logspace(math.log10(2.5e-6), math.log10(10e-6), 60)
    w = np.exp(-(np.log(d / 4.0e-6)) ** 2 / (2 * math.log(2.2) ** 2))
    e = np.array([MEDIA.efficiency(x, U) for x in d])
    return float((w * e).sum() / w.sum())


def fine_eff(U):
    return MEDIA.mass_eff(U, 0.5e-6, 2.0, 2.5e-6)


def geometry(cfg):
    s = cfg["s"]
    a = cfg["arm"]
    return dict(h=0.15 * s, L=3.0 * s, zlo=-8.1 * s, zhi=3.0 * s,
                rotors=[(sx * a, sy * a) for sx in (-1, 1) for sy in (-1, 1)], R=cfg["D"] / 2,
                disk=(-0.04 * s, 0.04 * s), above=(0.0375 * s, 0.1125 * s), below=(-0.1125 * s, -0.0375 * s),
                Ls=0.075 * s)


def write_mesh_dicts(case, cfg, fpos):
    g = geometry(cfg)
    s, h, L = cfg["s"], g["h"], g["L"]
    nx, nz = int(round(2 * L / h)), int(round((g["zhi"] - g["zlo"]) / h))
    write(case, "system/blockMeshDict", f"""
convertToMeters 1;
vertices ( ({-L} {-L} {g['zlo']}) ({L} {-L} {g['zlo']}) ({L} {L} {g['zlo']}) ({-L} {L} {g['zlo']})
           ({-L} {-L} {g['zhi']}) ({L} {-L} {g['zhi']}) ({L} {L} {g['zhi']}) ({-L} {L} {g['zhi']}) );
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {nx} {nz}) simpleGrading (1 1 1) );
boundary ( open {{ type patch; faces ( (0 3 2 1) (4 5 6 7) (0 1 5 4) (1 2 6 5) (2 3 7 6) (3 0 4 7) ); }} );
""")
    for r, (b0, b1) in {"r1": ((-0.85, -0.85, -2.6), (0.85, 0.85, 0.7)),
                        "r2": ((-0.5, -0.5, -0.5), (0.5, 0.5, 0.3))}.items():
        write(case, f"system/topoSetDict.{r}", f"""
actions ( {{ name ref; type cellSet; action new; source boxToCell;
  sourceInfo {{ box ({b0[0] * s} {b0[1] * s} {b0[2] * s}) ({b1[0] * s} {b1[1] * s} {b1[2] * s}); }} }} );
""")
    write(case, "system/refineMeshDict", """
set ref; coordinateSystem global;
globalCoeffs { tan1 (1 0 0); tan2 (0 1 0); }
directions ( tan1 tan2 normal ); useHexTopology yes; geometricCut no; writeMesh no;
""")
    acts = []
    zones = [("disk", g["disk"])] + ([(fpos, g[fpos])] if fpos != "none" else [])
    for i, (x, y) in enumerate(g["rotors"]):
        for zn, (z0, z1) in zones:
            acts.append(f"""{{ name {zn}{i}; type cellSet; action new; source cylinderToCell;
  sourceInfo {{ p1 ({x} {y} {z0}); p2 ({x} {y} {z1}); radius {g['R']}; }} }}""")
    acts.append("{ name disk; type cellSet; action new; source cellToCell; sourceInfo { sets (disk0 disk1 disk2 disk3); } }")
    acts.append("{ name diskZone; type cellZoneSet; action new; source setToCellZone; sourceInfo { set disk; } }")
    if fpos != "none":
        # ducted housing: solid shroud around rotor + filter (modelled as an
        # impermeable porous annulus one cell thick)
        z0 = min(g["disk"][0], g[fpos][0]); z1 = max(g["disk"][1], g[fpos][1])
        for i, (x, y) in enumerate(g["rotors"]):
            acts.append(f"""{{ name shroud{i}; type cellSet; action new; source cylinderAnnulusToCell;
  sourceInfo {{ p1 ({x} {y} {z0}); p2 ({x} {y} {z1}); outerRadius {g['R'] + 0.05 * s}; innerRadius {g['R']}; }} }}""")
        acts.append("{ name shroud; type cellSet; action new; source cellToCell; sourceInfo { sets (shroud0 shroud1 shroud2 shroud3); } }")
        acts.append("{ name shroudZone; type cellZoneSet; action new; source setToCellZone; sourceInfo { set shroud; } }")
        acts.append("{ name filt; type cellSet; action new; source cellToCell; sourceInfo { sets ("
                    + " ".join(f"{fpos}{i}" for i in range(4)) + "); } }")
        acts.append("{ name filtZone; type cellZoneSet; action new; source setToCellZone; sourceInfo { set filt; } }")
    write(case, "system/topoSetDict", "actions (\n" + "\n".join(acts) + "\n);\n")


def write_setup(case, cfg, fpos, thrust_total, lam, end_time, write_interval):
    s = cfg["s"]
    fo = [f"""disks
{{
    type vectorSemiImplicitSource; active yes;
    vectorSemiImplicitSourceCoeffs
    {{
        timeStart 0; duration 1e9; selectionMode cellZone; cellZone diskZone; volumeMode absolute;
        injectionRateSuSp {{ U ((0 0 {-thrust_total / RHO_AIR:.6g}) 0); }}
    }}
}}"""]
    if fpos != "none":
        d_slab = D_MEDIA * MEDIA.L / (PLEAT * geometry(cfg)["Ls"])
        fo.append(f"""filter
{{
    type explicitPorositySource; active yes;
    explicitPorositySourceCoeffs
    {{
        selectionMode cellZone; cellZone filtZone; type DarcyForchheimer;
        DarcyForchheimerCoeffs
        {{
            d d [0 -2 0 0 0 0 0] ({1000 * d_slab:.6g} {1000 * d_slab:.6g} {d_slab:.6g});
            f f [0 -1 0 0 0 0 0] (0 0 0);
            coordinateSystem {{ origin (0 0 0); e1 (1 0 0); e2 (0 1 0); }}
        }}
    }}
}}""")
        fo.append("""shroud
{
    type explicitPorositySource; active yes;
    explicitPorositySourceCoeffs
    {
        selectionMode cellZone; cellZone shroudZone; type DarcyForchheimer;
        DarcyForchheimerCoeffs
        {
            d d [0 -2 0 0 0 0 0] (1e10 1e10 1e10);
            f f [0 -1 0 0 0 0 0] (0 0 0);
            coordinateSystem { origin (0 0 0); e1 (1 0 0); e2 (0 1 0); }
        }
    }
}""")
        if lam > 0:
            fo.append(f"""capture
{{
    type scalarSemiImplicitSource; active yes;
    scalarSemiImplicitSourceCoeffs
    {{
        timeStart 0; duration 1e9; selectionMode cellZone; cellZone filtZone; volumeMode specific;
        injectionRateSuSp {{ T (0 {-lam:.6g}); }}
    }}
}}""")
    write(case, "constant/fvOptions", "\n".join(fo) + "\n")
    write(case, "constant/transportProperties",
          f"transportModel Newtonian;\nnu {NU:.6g};\nbeta 0;\nTRef 0;\nPr 1.0;\nPrt 0.7;\n")
    write(case, "constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n",
          cls="uniformDimensionedVectorField")
    write(case, "constant/turbulenceProperties",
          "simulationType RAS;\nRAS { RASModel kOmegaSST; turbulence on; printCoeffs off; }\n")
    write(case, "system/controlDict", f"""
application buoyantBoussinesqSimpleFoam; startFrom latestTime; stopAt endTime; endTime {end_time};
deltaT 1; writeControl timeStep; writeInterval {write_interval}; purgeWrite 1; writeFormat ascii;
writePrecision 8; runTimeModifiable false;
""")
    write(case, "system/fvSchemes", """
ddtSchemes { default steadyState; }
gradSchemes { default cellLimited Gauss linear 1; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U);
  div(phi,k) bounded Gauss upwind; div(phi,omega) bounded Gauss upwind;
  div(phi,T) bounded Gauss limitedLinear 1;
  div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.5; }
wallDist { method meshWave; }
""")
    write(case, "system/fvSolution", """
solvers {
  p_rgh { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.05; }
  "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
  T { solver PBiCGStab; preconditioner DILU; tolerance 1e-9; relTol 0.05; maxIter 200; }
}
SIMPLE { nNonOrthogonalCorrectors 1; consistent yes; pRefCell 0; pRefValue 0;
  residualControl { p_rgh 1e-5; U 1e-6; T 1e-7; } }
relaxationFactors { fields { p_rgh 0.3; } equations { U 0.7; k 0.5; omega 0.5; T 0.7; } }
""")


def write_initial(case, cfg):
    s = cfg["s"]
    write_field(case, "0", field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)",
                                 {"open": "        type pressureInletOutletVelocity; value uniform (0 0 0);"}))
    write_field(case, "0", field("p_rgh", "volScalarField", "[0 2 -2 0 0 0 0]", "0",
                                 {"open": "        type totalPressure; p0 uniform 0; value uniform 0;"}))
    write_field(case, "0", field("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0",
                                 {"open": "        type calculated; value uniform 0;"}))
    for n, v, d in (("k", 1e-4, "[0 2 -2 0 0 0 0]"), ("omega", 0.1 / s, "[0 0 -1 0 0 0 0]"),
                    ("T", 1.0, "[0 0 0 1 0 0 0]")):
        write_field(case, "0", field(n, "volScalarField", d, v,
                                     {"open": f"        type inletOutlet; inletValue uniform {v}; value uniform {v};"}))
    for n in ("nut", "alphat"):
        write_field(case, "0", field(n, "volScalarField", "[0 2 -1 0 0 0 0]", "0",
                                     {"open": "        type calculated; value uniform 0;"}))


def stats(case, fpos, tag):
    d = zone_stats(vtk_zone(case, "diskZone", f"VTK_disk_{tag}"))
    out = dict(vdisk=abs(d["U_avg"][2]), disk_vol=d["V"])
    if fpos != "none":
        f = zone_stats(vtk_zone(case, "filtZone", f"VTK_filt_{tag}"))
        out.update(uz_f=abs(f["U_avg"][2]), filt_vol=f["V"], T_int=f["T_int"], T_avg=f["T_avg"])
    return out


def solve_case(key, fpos, n_main=1500, n_more=700):
    cfg = DRONES[key]
    g = geometry(cfg)
    case = os.path.join(RUNS, f"{key}_{fpos}")
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    W = (cfg["mass"] + (cfg["filt_mass"] if fpos != "none" else 0.0)) * G
    T = W
    Ls = g["Ls"]
    d_slab = D_MEDIA * MEDIA.L / (PLEAT * Ls)
    write_mesh_dicts(case, cfg, fpos)
    write_setup(case, cfg, fpos, T, 0.0, n_main, n_main)
    write_initial(case, cfg)
    run(case, "blockMesh")
    for r in ("r1", "r2"):
        run(case, f"topoSet -dict system/topoSetDict.{r}", f"topoSet.{r}")
        run(case, "refineMesh -overwrite -dict system/refineMeshDict", f"refineMesh.{r}")
    run(case, "topoSet")
    run(case, "buoyantBoussinesqSimpleFoam", "solve.0")
    hist = [dict(pass_="flow", thrust=T, **stats(case, fpos, "0"))]
    out = dict(case=f"{key}_{fpos}", drone=cfg["label"], filter=fpos, weight_N=W)
    t_end = n_main
    if fpos != "none":
        res = {}
        for size, eff_fn in (("PM25", fine_eff), ("coarse", coarse_eff)):
            h = hist[-1]
            F = MU_AIR * h["uz_f"] * d_slab * h["filt_vol"]    # N, Darcy load on filter
            T = W + F
            e = eff_fn(h["uz_f"] / PLEAT)
            lam = -math.log(max(1 - e, 1e-9)) * h["uz_f"] / Ls
            t_end += n_more
            write_setup(case, cfg, fpos, T, lam, t_end, n_more)
            run(case, "buoyantBoussinesqSimpleFoam", f"solve.{size}")
            st = stats(case, fpos, size)
            hist.append(dict(pass_=size, thrust=T, lam=lam, eff=e, **st))
            # CADR [m3/s]: removal rate per unit ambient concentration (zone mass balance)
            res[size] = dict(eff=e, CADR=lam * st["T_int"], C_in_filter=st["T_avg"])
        h = hist[-1]
        dp = MU_AIR * h["uz_f"] * Ls * d_slab
        c25, cc = res["PM25"]["CADR"], res["coarse"]["CADR"]
        c10 = (C_PM25 * c25 + (C_PM10 - C_PM25) * cc) / C_PM10
        out.update(face_vel=h["uz_f"], dp_Pa=dp, F_filter_N=MU_AIR * h["uz_f"] * d_slab * h["filt_vol"],
                   eff_PM25=res["PM25"]["eff"], eff_coarse=res["coarse"]["eff"],
                   C_filter_PM25=res["PM25"]["C_in_filter"],
                   CADR_PM25=c25, CADR_PM10=c10,
                   g_per_h_PM25=c25 * C_PM25 * 3.6e-3, g_per_h_PM10=c10 * C_PM10 * 3.6e-3)
    h = hist[-1]
    A = 4 * math.pi * g["R"] ** 2
    Ttot = h["thrust"]
    # power: actuator-disk power T*v_disk (includes filter loss, since the disk must
    # supply the filter pressure drop), corrected by figure of merit and coaxial loss
    p_aero = Ttot * h["vdisk"] / cfg["fm"] * cfg["coax"]
    out.update(thrust_N=Ttot, v_disk=h["vdisk"], v_momentum=math.sqrt(Ttot / (2 * RHO_AIR * A)),
               Q_rotor=A * h["vdisk"], power_W=p_aero / cfg["eta"],
               endurance_min=cfg["batt_wh"] / (p_aero / cfg["eta"]) * 60, history=hist)
    json.dump(out, open(os.path.join(case, "result.json"), "w"), indent=1, default=float)
    # full-field export of y = rotor-row plane for figures
    vtk_zone(case, None, "VTK_full")
    return out


def main(keys=None):
    os.makedirs(RUNS, exist_ok=True)
    jobs = [(k, f) for k in (keys or DRONES) for f in ("none", "above", "below")]
    with ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(lambda a: solve_case(*a), jobs))
    json.dump(results, open(os.path.join(RUNS, "summary.json"), "w"), indent=1, default=float)
    for r in results:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != "history"})


if __name__ == "__main__":
    if len(sys.argv) == 4:
        print(solve_case(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[3]) // 2))
    else:
        main()
