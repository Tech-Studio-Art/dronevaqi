#!/usr/bin/env python3
"""
Report 2 - three-dimensional fleet study (network effects).

Physics (same as the report-1 best estimate): buoyantPimpleFoam, buoyant k-epsilon with Cg = 1.5
(calibrated against the Monin-Obukhov flux-profile relations, step 1b), Boussinesq equation of
state and a large specific heat (T = potential temperature), Pr_t = Sc_t = 0.7. PM is read from
theta (same transport equation). Night inversion: +5 K over 0-100 m, neutral above; PM falls
linearly from the street value to 50 % at 100 m.

Drones: each drone is an actuator disk of one cell column (6.25 m x 6.25 m, 39 m2, two cells high)
with thrust = weight (50 kg, 490.5 N), hovering at 60 m. The far field of a jet depends on its
momentum flux; a larger disk gives a slower ("lazy") jet with less penetration. The disk is checked
against the resolved axisymmetric result in a periodic single-drone cell ('cell', 'cellR' cases):
a 7 m radius disk (156 m2) delays the onset by about 15 min, the one-column disk (39 m2) does not.

Cases
  cell, cell_fine   periodic 50 m x 50 m cell, no wind (verification of mesh + coarse disk)
  base              fleet domain, wind, no drones (spin-up 0-1800 s, then reference 1800-3600 s)
  F100              100 drones on a 50 m grid over 500 m x 500 m, continuous
  C50               50 drones (quasi-checkerboard group A), continuous
  R60, R600         50 drones that alternate between the A and B sites, period 60 s or 600 s
Fleet domain: half of the zone (mirror symmetry about the wind axis y = 0), x -297..1003 m,
y 0..403 m, z 0..200 m; the grid is shifted so that each drone site is a cell centre. Inlet: Richards & Hoxey log wind profile (U10 = 1 m/s, z0 = 0.1 m) and the
inversion profile of theta; outlet: hydrostatic p_rgh profile of the undisturbed inversion.
"""
import glob
import json
import math
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from foam import write, run, field, write_field  # noqa: E402

RUNS = os.path.join(HERE, "runs", os.environ.get("FLEET_DIR", "fleet3d"))
RHO0, NU = 1.25, 1.42e-5
TH0, DTH, HINV, PMR = 285.0, 5.0, 100.0, 0.5
CP_BIG, CG, PRT, Z0 = 1.0e5, 1.5, 0.7, 0.1
U10 = float(os.environ.get("FLEET_U10", "1.0"))
THRUST, ZH, RD, SP = 50 * 9.81, 60.0, 3.5, 50.0
T_SPIN, T_RUN = 1800.0, 1800.0
# Small ambient source in the epsilon equation (m2/s3). It stops a known defect of the buoyant k-epsilon
# model: in weakly unstable air with horizontal flow, the buoyancy term acts in the k equation but not in
# the epsilon equation, so k/epsilon can grow without limit (seen in case F100 at 2.5 min). The value is
# 1e4 times smaller than epsilon in the jets and in the surface layer.
EPS_AMB = float(os.environ.get("FLEET_EPSAMB", "0"))
# Turbulence model. buoyantKEpsilonR = buoyantKEpsilon with C3 = 1 in unstable air (Rodi 1987); identical in
# stable air. Built from turbulence/ (build.sh) into turbulence/lib/libbuoyantKEpsilonR.so.
RAS = os.environ.get("FLEET_RAS", "buoyantKEpsilonR")
RASLIB = os.path.join(HERE, "turbulence", "lib", "libbuoyantKEpsilonR.so")
# Lower bounds of k and epsilon. They stop the k-epsilon stiffness when k -> 0 with finite epsilon (seen above a
# drone disk that starts). The largest nut at the bounds is Cmu*KMIN^2/EPSMIN = 9e-6 m2/s (below the molecular
# viscosity of air), so the bounds do not mix the inversion.
KMIN = float(os.environ.get("FLEET_KMIN", "1e-6"))
EPSMIN = float(os.environ.get("FLEET_EPSMIN", "1e-8"))


def turbulence_dicts(case):
    """turbulenceProperties for RAS and the libs line of controlDict."""
    write(case, "constant/turbulenceProperties",
          f"simulationType RAS;\nRAS {{ RASModel {RAS}; turbulence on; printCoeffs on; kMin {KMIN:g}; epsilonMin {EPSMIN:g}; "
          f"{RAS}Coeffs {{ Cg {CG}; Prt {PRT}; C3u 1.0; }} }}\n")
    libs = '"libatmosphericModels.so"' + (f' "{RASLIB}"' if RAS == "buoyantKEpsilonR" else "")
    cd = os.path.join(case, "system/controlDict")
    if os.path.exists(cd):
        s = re.sub(r"libs \([^)]*\);", f"libs ({libs});", open(cd).read())
        open(cd, "w").write(s)
KAPPA, CMU = 0.41, 0.09


def theta(z):
    return TH0 + DTH * np.minimum(z, HINV) / HINV


def pm(th):
    return 1.0 - (1.0 - PMR) * (th - TH0) / DTH


def ustar(u10=U10):
    return KAPPA * u10 / math.log((10 + Z0) / Z0)


def prgh(z):
    """Hydrostatic p_rgh of the undisturbed inversion (Boussinesq EoS, hRef = 0)."""
    c = RHO0 / TH0 * 9.81 * DTH / HINV
    zz = np.minimum(z, HINV)
    return 1e5 - c * zz ** 2 / 2


# ------------------------------------------------------------------ mesh
ZSEG = [(50.0, 20, 3.0), (22.0, 8, 1.0), (128.0, 16, 4.0)]          # 0-50, 50-72, 72-200 m


def grading(segs):
    L = sum(s[0] for s in segs)
    n = sum(s[1] for s in segs)
    g = " ".join(f"({s[0] / L:.10f} {s[1] / n:.10f} {s[2]})" for s in segs)
    return L, n, f"({g})"


def block_mesh(case, x0, xsegs, y0, ysegs, patches):
    Lx, nx, gx = grading(xsegs)
    Ly, ny, gy = grading(ysegs)
    Lz, nz, gz = grading(ZSEG)
    x1, y1 = x0 + Lx, y0 + Ly
    v = [(x0, y0, 0), (x1, y0, 0), (x1, y1, 0), (x0, y1, 0), (x0, y0, Lz), (x1, y0, Lz), (x1, y1, Lz), (x0, y1, Lz)]
    vs = " ".join(f"({a} {b} {c})" for a, b, c in v)
    faces = {"xmin": "(0 4 7 3)", "xmax": "(1 2 6 5)", "ymin": "(0 1 5 4)", "ymax": "(3 7 6 2)",
             "ground": "(0 3 2 1)", "top": "(4 5 6 7)"}
    bnd = []
    for name, (typ, fkeys, extra) in patches.items():
        bnd.append(f"  {name} {{ type {typ}; {extra} faces ( {' '.join(faces[k] for k in fkeys)} ); }}")
    write(case, "system/blockMeshDict", f"""
convertToMeters 1;
vertices ( {vs} );
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading ({gx} {gy} {gz}) );
boundary (
{chr(10).join(bnd)}
);
""")
    return nx * ny * nz


# ------------------------------------------------------------------ drones
def fleet_sites():
    """Half zone (y > 0): columns i = 0..9 (x = 25..475), rows r = 0..4 (y = 25..225).
    Group A = (i + r) even, B = odd: a checkerboard that is mirror-symmetric about y = 0."""
    return [(25 + SP * i, 25 + SP * r, "A" if (i + r) % 2 == 0 else "B") for i in range(10) for r in range(5)]


def zones(case, sites, rd=None):
    rd = RD if rd is None else rd
    acts = []
    groups = {}
    for k, (x, y, g) in enumerate(sites):
        acts.append(f"""{{ name d{k}; type cellSet; action new; source cylinderToCell;
  sourceInfo {{ p1 ({x} {y} {ZH - 3}); p2 ({x} {y} {ZH + 3}); radius {rd}; }} }}""")
        groups.setdefault(g, []).append(f"d{k}")
    for g, names in groups.items():
        acts.append(f"{{ name grp{g}; type cellSet; action new; source cellToCell; sourceInfo {{ sets ({' '.join(names)}); }} }}")
        acts.append(f"{{ name grp{g}Zone; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set grp{g}; }} }}")
    write(case, "system/topoSetDict", "actions (\n" + "\n".join(acts) + "\n);\n")
    return {g: len(n) for g, n in groups.items()}


def drone_options(case, counts, schedule):
    """schedule: list of (group, t_start, duration). Force in N (compressible momentum equation)."""
    fo = []
    for j, (g, t0, dur) in enumerate(schedule):
        fo.append(f"""on{j}_{g}
{{
    type vectorSemiImplicitSource; active yes;
    vectorSemiImplicitSourceCoeffs
    {{
        timeStart {t0:g}; duration {dur:g}; selectionMode cellZone; cellZone grp{g}Zone; volumeMode absolute;
        injectionRateSuSp {{ U ((0 0 {-THRUST * counts[g]:.6g}) 0); }}
    }}
}}""")
    if EPS_AMB > 0:
        fo.append(f"""ambientEpsilon
{{
    type scalarSemiImplicitSource; active yes;
    scalarSemiImplicitSourceCoeffs
    {{
        selectionMode all; volumeMode specific;
        injectionRateSuSp {{ epsilon ({RHO0 * EPS_AMB:.6g} 0); }}
    }}
}}""")
    write(case, "constant/fvOptions", "\n".join(fo) + "\n")


CLUSTERS = {   # half domain (y > 0); boxes select cell centres; drones = drones inside the half domain
    "K4": dict(boxes=[((-3.0, 47.0), (28.0, 78.0)), ((-3.0, 172.0), (28.0, 203.0))], drones=50,
               text="4 clusters of 25 drones (5 x 5 cells, 31 m x 31 m) at the upwind edge, 125 m apart across the wind"),
    "K1": dict(boxes=[((-3.0, 0.0), (59.0, 34.0))], drones=50,
               text="1 cluster of 100 drones (63 m x 69 m) at the upwind edge, on the wind axis"),
}


def cluster_zones(case, kind):
    acts = []
    for n, ((x0, y0), (x1, y1)) in enumerate(CLUSTERS[kind]["boxes"]):
        acts.append(f"{{ name k{n}; type cellSet; action new; source boxToCell; "
                    f"sourceInfo {{ box ({x0} {y0} {ZH - 3}) ({x1} {y1} {ZH + 3}); }} }}")
    names = " ".join(f"k{n}" for n in range(len(CLUSTERS[kind]["boxes"])))
    acts.append(f"{{ name grpK; type cellSet; action new; source cellToCell; sourceInfo {{ sets ({names}); }} }}")
    acts.append("{ name grpKZone; type cellZoneSet; action new; source setToCellZone; sourceInfo { set grpK; } }")
    write(case, "system/topoSetDict", "actions (\n" + "\n".join(acts) + "\n);\n")


def schedule_for(kind, t0, t1):
    if kind == "none":
        return []
    if kind in CLUSTERS:
        return [("K", t0, 1e9)]
    if kind == "F100":
        return [("A", t0, 1e9), ("B", t0, 1e9)]
    if kind == "C50":
        return [("A", t0, 1e9)]
    period = {"R60": 60.0, "R600": 600.0}[kind]
    out, t, g = [], t0, "A"
    while t < t1 + 1:
        out.append((g, t, period / 2))
        g = "B" if g == "A" else "A"
        t += period / 2
    return out


# ------------------------------------------------------------------ case set-up
def dicts(case, end_time, write_interval, nproc):
    write(case, "constant/thermophysicalProperties", f"""
thermoType {{ type heRhoThermo; mixture pureMixture; transport const; thermo hConst;
  equationOfState Boussinesq; specie specie; energy sensibleEnthalpy; }}
dpdt no;
mixture {{
  specie {{ molWeight 28.96; }}
  equationOfState {{ rho0 {RHO0}; T0 {TH0}; beta {1 / TH0:.8g}; }}
  thermodynamics {{ Cp {CP_BIG}; Hf 0; }}
  transport {{ mu {NU * RHO0:.6g}; Pr 0.7; }}
}}
""")
    write(case, "constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n", cls="uniformDimensionedVectorField")
    write(case, "system/controlDict", f"""
application buoyantPimpleFoam; startFrom latestTime; stopAt endTime; endTime {end_time:g};
deltaT 0.5; writeControl adjustableRunTime; writeInterval {write_interval:g}; purgeWrite 0; writeFormat binary;
writePrecision 8; runTimeModifiable false; adjustTimeStep yes; maxCo {float(os.environ.get("FLEET_MAXCO", "1.5")):g}; maxDeltaT 2;
libs ("libatmosphericModels.so");
""")
    turbulence_dicts(case)
    write(case, "system/fvSchemes", """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss limitedLinearV 1; div(phi,h) Gauss limitedLinear 1;
  div(phi,K) Gauss linear; div(phi,k) Gauss limitedLinear 1; div(phi,epsilon) Gauss limitedLinear 1;
  div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write(case, "system/fvSolution", """
solvers {
  "rho.*" { solver diagonal; }
  p_rgh { solver GAMG; smoother DIC; tolerance 1e-7; relTol 0.01; }
  p_rghFinal { $p_rgh; relTol 0; }
  "(U|k|epsilon)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
  "(U|k|epsilon)Final" { $U; relTol 0; }
  h { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0.001; maxIter 200; }
  hFinal { $h; relTol 0; }
}
PIMPLE { momentumPredictor yes; nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0;
  pRefCell 0; pRefValue 1e5; }
""")
    write(case, "system/decomposeParDict", decomp(nproc))


def decomp(nproc):
    return f"numberOfSubdomains {nproc};\nmethod simple;\nsimpleCoeffs {{ n ({nproc} 1 1); delta 0.001; }}\n"


def fields(case, kind):
    us = ustar()
    k0, e0 = us ** 2 / math.sqrt(CMU), us ** 3 / (KAPPA * 10)
    abl = (f"flowDir (1 0 0); zDir (0 0 1); Uref {U10}; Zref 10; z0 uniform {Z0}; zGround uniform 0; "
           f"kappa {KAPPA}; Cmu {CMU};")
    th_tab = "table ( (0 %g) (%g %g) (1000 %g) )" % (TH0, HINV, TH0 + DTH, TH0 + DTH)
    zz = np.linspace(0, HINV, 51)
    p_tab = "table ( " + " ".join(f"({z:g} {prgh(z):.6f})" for z in zz) + f" (1000 {prgh(1000):.6f}) )"
    if kind == "cell":
        side = {"xmin": "cyclic", "xmax": "cyclic", "ymin": "cyclic", "ymax": "cyclic"}
        S = {p: f"        type {t};" for p, t in side.items()}
        U = dict(S, ground="        type noSlip;", top="        type slip;")
        P = dict(S, ground="        type fixedFluxPressure; value uniform 1e5;",
                 top="        type fixedFluxPressure; value uniform 1e5;")
        T = dict(S, ground="        type zeroGradient;", top="        type zeroGradient;")
        K = dict(S, ground=f"        type kqRWallFunction; value uniform 1e-5;", top="        type zeroGradient;")
        E = dict(S, ground=f"        type epsilonWallFunction; value uniform 1e-7;", top="        type zeroGradient;")
        NT = dict(S, ground=f"        type nutkWallFunction; value uniform 0;", top="        type calculated; value uniform 0;")
        AT = dict(S, ground=f"        type compressible::alphatJayatillekeWallFunction; Prt {PRT}; value uniform 0;",
                  top="        type calculated; value uniform 0;")
        u0, kk, ee = "(0 0 0)", "1e-5", "1e-7"
    else:
        S = {"ymin": "        type symmetry;", "ymax": "        type symmetry;"}
        U = dict(S, inlet=f"        type atmBoundaryLayerInletVelocity; {abl} value uniform ({U10} 0 0);",
                 outlet="        type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0);",
                 ground="        type noSlip;", top="        type slip;")
        P = dict(S, inlet="        type fixedFluxPressure; value uniform 1e5;",
                 outlet=f"        type fixedProfile; profile {p_tab}; direction (0 0 1); origin 0;",
                 ground="        type fixedFluxPressure; value uniform 1e5;",
                 top="        type fixedFluxPressure; value uniform 1e5;")
        T = dict(S, inlet=f"        type fixedProfile; profile {th_tab}; direction (0 0 1); origin 0;",
                 outlet="        type zeroGradient;", ground="        type zeroGradient;", top="        type zeroGradient;")
        K = dict(S, inlet=f"        type atmBoundaryLayerInletK; {abl} value uniform {k0};",
                 outlet=f"        type inletOutlet; inletValue uniform {k0}; value uniform {k0};",
                 ground=f"        type kqRWallFunction; value uniform {k0};", top="        type zeroGradient;")
        E = dict(S, inlet=f"        type atmBoundaryLayerInletEpsilon; {abl} value uniform {e0};",
                 outlet=f"        type inletOutlet; inletValue uniform {e0}; value uniform {e0};",
                 ground=f"        type epsilonWallFunction; value uniform {e0};", top="        type zeroGradient;")
        NT = dict(S, inlet="        type calculated; value uniform 0;", outlet="        type calculated; value uniform 0;",
                  ground=f"        type nutkAtmRoughWallFunction; z0 uniform {Z0}; value uniform 0;",
                  top="        type calculated; value uniform 0;")
        AT = dict(S, inlet="        type calculated; value uniform 0;", outlet="        type calculated; value uniform 0;",
                  ground=f"        type compressible::alphatJayatillekeWallFunction; Prt {PRT}; value uniform 0;",
                  top="        type calculated; value uniform 0;")
        u0, kk, ee = f"({U10} 0 0)", f"{k0}", f"{e0}"
    write_field(case, "0", field("U", "volVectorField", "[0 1 -1 0 0 0 0]", u0, U))
    write_field(case, "0", field("p_rgh", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5", P))
    write_field(case, "0", field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5",
                                 {p: ("        type calculated; value uniform 1e5;" if "cyclic" not in b and "symmetry" not in b else b)
                                  for p, b in P.items()}))
    write_field(case, "0", field("T", "volScalarField", "[0 0 0 1 0 0 0]", f"{TH0}", T))
    write_field(case, "0", field("k", "volScalarField", "[0 2 -2 0 0 0 0]", kk, K))
    write_field(case, "0", field("epsilon", "volScalarField", "[0 2 -3 0 0 0 0]", ee, E))
    write_field(case, "0", field("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "0", NT))
    write_field(case, "0", field("alphat", "volScalarField", "[1 -1 -1 0 0 0 0]", "0", AT))


def set_initial(case, kind):
    """Nonuniform initial theta (and, with wind, the log wind and k-epsilon profiles) from cell centres."""
    import pyvista as pv
    run(case, "foamToVTK -time 0 -no-point-data -no-boundary -overwrite -name VTK_mesh -fields '(T)'", "vtk.mesh")
    g = pv.read(glob.glob(os.path.join(case, "VTK_mesh", "*", "internal.vtu"))[0])
    z = np.asarray(g.cell_centers().points)[:, 2]
    vals = {"T": ("List<scalar>", [f"{v:.6f}" for v in theta(z)]),
            "p_rgh": ("List<scalar>", [f"{v:.6f}" for v in prgh(z)])}
    if kind != "cell":
        us = ustar()
        u = us / KAPPA * np.log((z + Z0) / Z0)
        vals["U"] = ("List<vector>", [f"({v:.6f} 0 0)" for v in u])
        vals["epsilon"] = ("List<scalar>", [f"{us ** 3 / (KAPPA * (zz + Z0)):.8g}" for zz in z])
    for name, (kind_, lst) in vals.items():
        p = os.path.join(case, "0", name)
        s = open(p).read()
        body = f"nonuniform {kind_}\n{len(lst)}\n(\n" + "\n".join(lst) + "\n)"
        s = re.sub(r"internalField\s+uniform [^;]+;", f"internalField   {body};", s, count=1)
        open(p, "w").write(s)
    return len(z)


def solve_parallel(case, nproc, log):
    if nproc > 1:
        run(case, "decomposePar -force -latestTime", f"decompose.{log}")
        run(case, f"mpirun --oversubscribe --bind-to none -np {nproc} buoyantPimpleFoam -parallel", log)
        run(case, "reconstructPar -newTimes", f"reconstruct.{log}")
        for d in glob.glob(os.path.join(case, "processor*")):
            shutil.rmtree(d)
    else:
        run(case, "buoyantPimpleFoam", log)


# ------------------------------------------------------------------ cases
def build_cell(name, fine=False, n=None, rd=None, pos=(25.0, 25.0)):
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    n = n or (16 if fine else 8)
    cyc = lambda a, b: ("cyclic", [a], f"neighbourPatch {b};")
    ncell = block_mesh(case, 0.0, [(50.0, n, 1.0)], 0.0, [(50.0, n, 1.0)],
                       {"xmin": cyc("xmin", "xmax"), "xmax": cyc("xmax", "xmin"),
                        "ymin": cyc("ymin", "ymax"), "ymax": cyc("ymax", "ymin"),
                        "ground": ("wall", ["ground"], ""), "top": ("wall", ["top"], "")})
    counts = zones(case, [(pos[0], pos[1], "A")], rd)
    dicts(case, T_RUN, 300, 1)
    drone_options(case, counts, [("A", 0, 1e9)])
    fields(case, "cell")
    run(case, "blockMesh")
    run(case, "topoSet")
    set_initial(case, "cell")
    return case, ncell


def build_pcell(name, kind):
    """Calm-air rhythm test: periodic 100 m x 100 m cell with the A/B checkerboard of the fleet (4 sites on
    cell centres). kind: F100 (all sites), C50 (A sites), R60 / R600 (A and B alternate)."""
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    cyc = lambda a, b: ("cyclic", [a], f"neighbourPatch {b};")
    ncell = block_mesh(case, 0.0, [(100.0, 16, 1.0)], 0.0, [(100.0, 16, 1.0)],
                       {"xmin": cyc("xmin", "xmax"), "xmax": cyc("xmax", "xmin"),
                        "ymin": cyc("ymin", "ymax"), "ymax": cyc("ymax", "ymin"),
                        "ground": ("wall", ["ground"], ""), "top": ("wall", ["top"], "")})
    c0, c1 = 28.125, 78.125
    counts = zones(case, [(c0, c0, "A"), (c1, c1, "A"), (c1, c0, "B"), (c0, c1, "B")], RD)
    dicts(case, T_RUN, 300, 1)
    drone_options(case, counts, schedule_for(kind, 0.0, T_RUN))
    fields(case, "cell")
    run(case, "blockMesh")
    run(case, "topoSet")
    set_initial(case, "cell")
    return case, ncell


def build_base(nproc):
    case = os.path.join(RUNS, "base")
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    # x0 and the half-width first row put every site (x = 25 + 50 i, y = 25 + 50 r) on a cell centre
    ncell = block_mesh(case, -296.875, [(200.0, 16, 1.0), (900.0, 144, 1.0), (200.0, 16, 1.0)],
                       0.0, [(3.125, 1, 1.0), (300.0, 48, 1.0), (100.0, 8, 1.0)],
                       {"inlet": ("patch", ["xmin"], ""), "outlet": ("patch", ["xmax"], ""),
                        "ymin": ("symmetry", ["ymin"], ""), "ymax": ("symmetry", ["ymax"], ""),
                        "ground": ("wall", ["ground"], ""), "top": ("patch", ["top"], "")})
    counts = zones(case, fleet_sites())
    dicts(case, T_SPIN, 300, nproc)
    drone_options(case, counts, [])
    fields(case, "fleet")
    run(case, "blockMesh")
    run(case, "topoSet")
    set_initial(case, "fleet")
    return case, ncell, counts


def branch(kind, nproc):
    """Copy the spun-up base case and run one fleet configuration from T_SPIN to T_SPIN + T_RUN."""
    base = os.path.join(RUNS, "base")
    case = os.path.join(RUNS, kind)
    if os.environ.get("FLEET_RESUME") and os.path.isdir(os.path.join(case, "constant")):
        n = len(glob.glob(os.path.join(case, "log.solve.resume*")))
        solve_parallel(case, nproc, f"solve.resume{n + 1}")  # startFrom latestTime
        return case
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    for d in ("constant", "system", f"{T_SPIN:g}"):
        shutil.copytree(os.path.join(base, d), os.path.join(case, d))
    counts = {"A": sum(1 for s in fleet_sites() if s[2] == "A"), "B": sum(1 for s in fleet_sites() if s[2] == "B")}
    turbulence_dicts(case)
    if kind in CLUSTERS:
        cluster_zones(case, kind)
        run(case, "topoSet", "topoSet")
        counts = {"K": CLUSTERS[kind]["drones"]}
    drone_options(case, counts, schedule_for(kind, T_SPIN, T_SPIN + T_RUN))
    cd = os.path.join(case, "system/controlDict")
    s = re.sub(r"endTime [0-9.e+]+;", f"endTime {T_SPIN + T_RUN:g};", open(cd).read())
    if os.environ.get("FLEET_MAXCO"):                     # smaller time step (stability), same model
        s = re.sub(r"maxCo [0-9.e+]+;", f"maxCo {float(os.environ['FLEET_MAXCO']):g};", s)
    open(cd, "w").write(s)
    write(case, "system/decomposeParDict", decomp(nproc))
    solve_parallel(case, nproc, "solve")
    return case


# ------------------------------------------------------------------ analysis
def _load_times(case, name, fields=("T",), times=None):
    import pyvista as pv
    sel = f"-time '{times}'" if times else ""
    run(case, f"foamToVTK {sel} -no-point-data -no-boundary -overwrite -name {name} -fields '({' '.join(fields)})'",
        f"vtk.{name}")
    fs = sorted(glob.glob(os.path.join(case, name, "*", "internal.vtu")),
                key=lambda f: int(f.split("_")[-1].split("/")[0]) if f.split("_")[-1].split("/")[0].isdigit() else 0)
    tdirs = sorted(float(d) for d in os.listdir(case) if re.fullmatch(r"[0-9.e+]+", d))
    if times:
        lo, hi = [float(x) if x else None for x in times.split(":")]
        tdirs = [t for t in tdirs if (lo is None or t >= lo - 1e-6) and (hi is None or t <= hi + 1e-6)]
    out = []
    for f, t in zip(fs, tdirs):
        g = pv.read(f).compute_cell_sizes(length=False, area=False, volume=True)
        out.append((t, g))
    return out


def analyse_cell(name):
    case = os.path.join(RUNS, name)
    series = []
    for t, g in _load_times(case, "VTK_all"):
        z = np.asarray(g.cell_centers().points)[:, 2]
        V = np.asarray(g.cell_data["Volume"])
        C = pm(np.asarray(g.cell_data["T"]))
        rec = {"t": t}
        for lab, zm in (("z<2m", 2.0), ("z<10m", 10.0)):
            m = z < zm
            rec[lab] = float((C[m] * V[m]).sum() / V[m].sum())
        series.append(rec)
    res = dict(case=name, series=series, decrease_30min_pct=100 * (1 - series[-1]["z<2m"] / series[0]["z<2m"]))
    json.dump(res, open(os.path.join(case, "result.json"), "w"), indent=1)
    return res


# ------------------------------------------------------------------ fleet analysis
ZONE_X, ZONE_Y, DOWN_X = (0.0, 500.0), 250.0, (500.0, 1000.0)


def _axes(g):
    """Axis centres and (i, j, k) index of each cell of the single-block structured mesh."""
    c = np.asarray(g.cell_centers().points)
    ax, ix = [], []
    for d in range(3):
        u, inv = np.unique(np.round(c[:, d], 3), return_inverse=True)
        ax.append(u)
        ix.append(inv)
    return ax, ix


def extract_fleet(kind):
    """PM (from theta) on the structured grid at each output time -> <case>/pm.npz."""
    case = os.path.join(RUNS, kind)
    data = _load_times(case, "VTK_all", fields=("T", "U"), times=f"{T_SPIN:g}:")
    (x, y, z), (i, j, k) = _axes(data[0][1])
    V = np.zeros((len(z), len(y), len(x)))
    V[k, j, i] = np.asarray(data[0][1].cell_data["Volume"])
    P = np.zeros((len(data), len(z), len(y), len(x)), dtype=np.float32)
    U = np.zeros((len(data), 3, len(z), len(y), len(x)), dtype=np.float32)
    for n, (t, g) in enumerate(data):
        P[n][k, j, i] = pm(np.asarray(g.cell_data["T"]))
        u = np.asarray(g.cell_data["U"])
        for d in range(3):
            U[n, d][k, j, i] = u[:, d]
    times = np.array([t for t, _ in data])
    np.savez_compressed(os.path.join(case, "pm.npz"), t=times, x=x, y=y, z=z, V=V, pm=P, U=U[-1])
    shutil.rmtree(os.path.join(case, "VTK_all"), ignore_errors=True)
    return times


def _red(P, P0, A, m):
    """Reduction (%) of the area-weighted mean PM in mask m, relative to the reference P0."""
    return float(100 * (1 - (P[m] * A[m]).sum() / (P0[m] * A[m]).sum()))


SPECS = {   # case -> reference case, number of drones (full zone), turbulence model of the run
    "fleet3d": {"F100": dict(ref="none", drones=100, model="final"),
                "K4": dict(ref="none", drones=100, model="final"),
                "K1": dict(ref="none", drones=100, model="final"),
                "C50": dict(ref="none_noamb", drones=50, model="original"),
                "R600": dict(ref="none_noamb", drones=50, model="original"),
                "R60": dict(ref="none_noamb", drones=50, model="original")},
    "fleet3d_u2": {"F100": dict(ref="none", drones=100, model="original"),
                   "K1": dict(ref="none", drones=100, model="final")},
}
SECTION_Y = {"K1": 12.5, "K4": 62.5}      # vertical section through the cluster; others: first drone row


def summarise(spec):
    """Metrics of each case relative to its no-drone reference at the same time (ground layer, 0-1.37 m)."""
    out = dict(U10=U10, cases={})
    refs = {}
    for kind, sp in spec.items():
        if not os.path.exists(os.path.join(RUNS, kind, "pm.npz")):
            continue
        if sp["ref"] not in refs:
            refs[sp["ref"]] = np.load(os.path.join(RUNS, sp["ref"], "pm.npz"))
        ref, d = refs[sp["ref"]], np.load(os.path.join(RUNS, kind, "pm.npz"))
        n = min(len(ref["t"]), len(d["t"]))
        assert np.allclose(ref["t"][:n], d["t"][:n]), kind
        t, x, y, z = ref["t"][:n], ref["x"], ref["y"], ref["z"]
        A = ref["V"][0] / (2 * z[0])
        P0, P = ref["pm"][:n, 0].astype(float), d["pm"][:n, 0].astype(float)
        X, Y = np.meshgrid(x, y)
        zone = (X > ZONE_X[0]) & (X < ZONE_X[1]) & (Y < ZONE_Y)
        down = (X > DOWN_X[0]) & (X < DOWN_X[1]) & (Y < ZONE_Y)
        out.update(x=x.tolist(), y=y.tolist(), z=z.tolist(), ground_layer_m=float(2 * z[0]))
        k60 = int(np.argmin(abs(z - ZH)))                 # wind at the drones (reference, zone mean)
        out["U60"] = float(ref["U"][0][k60][zone].mean())
        series = [_red(P[m], P0[m], A, zone) for m in range(n)]
        rmap = 100 * (1 - P[-1] / P0[-1])
        my = Y < ZONE_Y
        rx = 100 * (1 - (P[-1] * A * my).sum(0) / (P0[-1] * A * my).sum(0))
        mx = (X > ZONE_X[0]) & (X < ZONE_X[1])
        ry = 100 * (1 - (P[-1] * A * mx).sum(1) / (P0[-1] * A * mx).sum(1))
        win = {f"{x0:g}": _red(P[-1], P0[-1], A, (X > x0) & (X < x0 + 500) & (Y < ZONE_Y))
               for x0 in np.arange(-250.0, 501.0, 25.0)}
        best = max(win, key=win.get)
        load = 2 * (P0[-1][zone] * A[zone]).sum()
        gain = 2 * ((P0[-1] - P[-1]) * A)
        parts = {"zone": zone, "downwind": (X > ZONE_X[1]), "upwind": (X < ZONE_X[0]),
                 "beside": (X > ZONE_X[0]) & (X < ZONE_X[1]) & (Y > ZONE_Y)}
        benefit = {k: float(100 * gain[m].sum() / load) for k, m in parts.items()}
        q = lambda m: _red(P[-1], P0[-1], A, m & zone)
        # zone-mean vertical profile of the PM decrease (last time) and the bottom of the mixed layer: going down
        # from the height of the largest decrease below the drones, the first height with less than 1 %
        V3 = ref["V"]
        P3, P30 = d["pm"][n - 1].astype(float), ref["pm"][n - 1].astype(float)
        prof_z = [float(100 * (1 - (P3[k][zone] * V3[k][zone]).sum() / (P30[k][zone] * V3[k][zone]).sum()))
                  for k in range(len(z))]
        below = [k for k in range(len(z)) if z[k] < ZH - 3]
        kmax = max(below, key=lambda k: prof_z[k])
        k = kmax
        while k > 0 and prof_z[k - 1] >= 1.0:
            k -= 1
        z_mix = float(z[k] - (z[k] - z[k - 1]) / 2) if k > 0 else 0.0
        jy = int(np.argmin(abs(y - SECTION_Y.get(kind, 25.0))))
        out["cases"][kind] = dict(
            ref=sp["ref"], drones=sp["drones"], model=sp["model"], t_min=((t - T_SPIN) / 60).tolist(),
            series=series, zone_last=series[-1], minutes=float((t[-1] - T_SPIN) / 60),
            zone_mean=float(np.trapezoid(series, t) / (t[-1] - t[0])),
            upwind_half=q(X < 250), downwind_half=q(X > 250), centre_strip=q(Y < 125), edge_strip=q(Y > 125),
            downwind_500m=_red(P[-1], P0[-1], A, down), max_ground=float(rmap.max()),
            best_window_x0=float(best), best_window=win[best], window=win,
            area_gt10_ha=float(2 * A[rmap > 10].sum() / 1e4), area_gt5_ha=float(2 * A[rmap > 5].sum() / 1e4),
            area_gt2_ha=float(2 * A[rmap > 2].sum() / 1e4),
            benefit_zone_equiv=benefit, benefit_total=float(sum(benefit.values())),
            mix_bottom_m=z_mix, profile_z=prof_z,
            profile_x=rx.tolist(), profile_y=ry.tolist(), map=np.round(rmap, 2).tolist(),
            section_y=float(y[jy]),
            section_ref=np.round(ref["pm"][n - 1][:, jy, :].astype(float), 4).tolist(),
            section_pm=np.round(d["pm"][n - 1][:, jy, :].astype(float), 4).tolist(),
            section_u=np.round(d["U"][0][:, jy, :].astype(float), 3).tolist(),
            section_w=np.round(d["U"][2][:, jy, :].astype(float), 3).tolist())
    # reference check: final model against the original model (no drones)
    if os.path.exists(os.path.join(RUNS, "none_noamb", "pm.npz")) and os.path.exists(os.path.join(RUNS, "none", "pm.npz")):
        a, b = np.load(os.path.join(RUNS, "none", "pm.npz")), np.load(os.path.join(RUNS, "none_noamb", "pm.npz"))
        n = min(len(a["t"]), len(b["t"]))
        A = a["V"][0] / (2 * a["z"][0])
        X, Y = np.meshgrid(a["x"], a["y"])
        zone = (X > 0) & (X < 500) & (Y < 250)
        out["reference_check"] = [100 * (1 - (a["pm"][m, 0][zone] * A[zone]).sum() / (b["pm"][m, 0][zone] * A[zone]).sum())
                                  for m in range(n)]
    json.dump(out, open(os.path.join(RUNS, "summary.json"), "w"), separators=(",", ":"))
    return out

def runinfo(kinds):
    """Steps, time step and CPU time of each case (from the solver logs) -> runs/<dir>/runinfo.json."""
    out = {}
    for kind in ["base"] + list(kinds):
        case = os.path.join(RUNS, kind)
        logs = [os.path.join(case, n) for n in ("log.solve.spin", "log.solve") if os.path.exists(os.path.join(case, n))]
        logs += sorted(glob.glob(os.path.join(case, "log.solve.resume*")))
        steps, dts, cpu = 0, [], 0.0
        for lg in logs:
            txt = open(lg).read()
            steps += len(re.findall(r"^Time = ", txt, re.M))
            dts += [float(v) for v in re.findall(r"^deltaT = ([0-9.e+-]+)", txt, re.M)]
            ex = re.findall(r"ExecutionTime = ([0-9.e+]+) s", txt)
            cpu += float(ex[-1]) if ex else 0.0
        nproc = 2 if kind == "base" else 1
        out[kind] = dict(steps=steps, dt_mean_s=float(np.mean(dts)) if dts else None,
                         dt_min_s=float(np.min(dts)) if dts else None, wall_s=cpu, cpu_core_h=cpu * nproc / 3600)
    m = re.search(r"nCells: *([0-9]+)", open(os.path.join(RUNS, "base", "log.blockMesh")).read())
    out["cells"] = int(m.group(1)) if m else None
    json.dump(out, open(os.path.join(RUNS, "runinfo.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "cellR":
        # disk-size study: argv = cellR <name> <cells per 50 m> <disk radius>
        nm, n, rd = sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
        off = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0      # shift so the drone sits on a cell centre
        case, nc = build_cell(nm, n=n, rd=rd, pos=(25.0 + off, 25.0 + off))
        print(nm, nc, "cells")
        run(case, "buoyantPimpleFoam", "solve")
        r = analyse_cell(nm)
        print(nm, round(r["decrease_30min_pct"], 1), [round(100 * x["z<2m"] / r["series"][0]["z<2m"], 1) for x in r["series"]])
    elif what == "pcell":
        # calm-air rhythm test: argv = pcell <kind ...>
        def one(kind):
            case, n = build_pcell("p" + kind, kind)
            run(case, "buoyantPimpleFoam", "solve")
            r = analyse_cell("p" + kind)
            return kind, n, r
        with ThreadPoolExecutor(max_workers=len(sys.argv[2:])) as ex:
            for kind, n, r in ex.map(one, sys.argv[2:]):
                print("p" + kind, n, round(r["decrease_30min_pct"], 2), [round(100 * x["z<2m"] / r["series"][0]["z<2m"], 1) for x in r["series"]])
    elif what == "cell":
        for nm, fine in (("cell", False), ("cell_fine", True)):
            case, n = build_cell(nm, fine, rd=7.0)
            print(nm, n, "cells")
            run(case, "buoyantPimpleFoam", "solve")
    elif what == "base":
        nproc = int(os.environ.get("NPROC", "4"))
        case, n, counts = build_base(nproc)
        print("base", n, "cells", counts)
        solve_parallel(case, nproc, "solve.spin")
    elif what == "branch":
        kinds = sys.argv[2:]
        nproc = int(os.environ.get("NPROC", "1"))
        with ThreadPoolExecutor(max_workers=len(kinds)) as ex:
            list(ex.map(lambda k: branch(k, nproc), kinds))
    elif what == "extract":
        for k in sys.argv[2:]:
            print(k, extract_fleet(k))
    elif what == "runinfo":
        print(json.dumps(runinfo(sys.argv[2:]), indent=1))
    elif what == "summary":
        r = summarise(SPECS[os.path.basename(RUNS)])
        for k, v in r["cases"].items():
            print(k, v["model"], f"{v['minutes']:.0f} min", "zone", round(v["zone_last"], 2), [round(a, 2) for a in v["series"]],
                  "down", round(v["downwind_500m"], 2), "max", round(v["max_ground"], 1), "best x0", v["best_window_x0"],
                  round(v["best_window"], 2), "mixed layer bottom", round(v["mix_bottom_m"], 1))
        print("reference check", r.get("reference_check"))
