#!/usr/bin/env python3
"""
Solution 1 -- drones as inversion mixers (bring cleaner air down from aloft).

Method (literature-standard pieces only):
  * Transient URANS, Boussinesq buoyancy (buoyantBoussinesqPimpleFoam),
    standard k-epsilon (Launder & Spalding 1974) with wall functions.
    Standard k-epsilon has no buoyancy damping of turbulence in stable air,
    so it OVER-predicts mixing -> results are an optimistic upper bound.
  * Drone = uniformly loaded actuator disk (momentum source, prescribed
    thrust = weight), all rotors lumped into one equivalent disk of equal total
    area (rotor wakes merge within a few diameters).
  * Axisymmetric "unit cell" (wedge): 100 drones on a 50 m grid over a
    500 x 500 m zone -> each drone owns a cylinder of equal area
    (R = 28.2 m) with a symmetry (slip) outer boundary, the classic unit-cell
    treatment of periodic arrays (cf. wind-farm / frost-fan array studies).
  * Night-time stable boundary layer: potential temperature rises linearly by
    dTheta over the inversion depth h, neutral above (residual layer).
  * Pollutant: passive scalar with Sc_t = Pr_t = 0.7. Because PM and
    potential temperature obey the same linear advection-diffusion equation
    with identical (zero-flux) boundary conditions, an initial PM profile that
    is a linear function of theta stays the same linear function of theta for
    all time. We prescribe PM decreasing linearly through the inversion
    (surface C_s -> r*C_s at z = h, uniform above) and recover PM from theta.
    Settling of PM10 (~5 mm/s, ~9 m per 30 min) is neglected.
  * Analytic check: max penetration of a momentum jet into linear
    stratification z_max ~ 2.9 (M/N^2)^(1/4) ... written M^(1/4) N^(-1/2)
    (Turner 1966; List 1982).
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
from foam import write, run, field, write_field, vtk_zone  # noqa: E402

RUNS = os.path.join(HERE, "runs", "mixing")
RHO, NU, THETA0 = 1.25, 1.42e-5, 285.0
WEDGE = 5.0  # degrees
R_CELL = math.sqrt(500 * 500 / 100 / math.pi)  # 28.2 m
H_DOM = 200.0
INV_H, D_THETA, PM_RATIO = 100.0, 5.0, 0.5    # inversion depth, strength, aloft/surface PM

CASES = {
    "heavy_z60": dict(T=50 * 9.81, Rd=math.sqrt(4 * 1.474 / math.pi), zh=60.0),
    "heavy_z90": dict(T=50 * 9.81, Rd=math.sqrt(4 * 1.474 / math.pi), zh=90.0),
    # quad: same thrust (momentum flux) on an enlarged 0.75 m disk -- far field of a
    # momentum jet depends on M = T/rho only (Turner 1966; List 1982); keeps dt tractable
    "quad_z30": dict(T=2.0 * 9.81, Rd=0.75, zh=30.0),
    "no_drone": dict(T=0.0, Rd=1.37, zh=60.0),
    # conservative case: stability-reduced turbulent transport of heat/PM,
    # Pr_t = Sc_t = 2.0 (stable-boundary-layer range ~1-3; e.g. Zilitinkevich et al. 2013)
    # grid-independence study (Celik et al. 2008): refinement ratio sqrt(2) around the base mesh
    "heavy_z60_gc": dict(T=50 * 9.81, Rd=math.sqrt(4 * 1.474 / math.pi), zh=60.0, refine=1 / math.sqrt(2)),
    "heavy_z60_gf": dict(T=50 * 9.81, Rd=math.sqrt(4 * 1.474 / math.pi), zh=60.0, refine=math.sqrt(2)),
    "heavy_z60_prt2": dict(T=50 * 9.81, Rd=math.sqrt(4 * 1.474 / math.pi), zh=60.0, Prt=2.0),
}
T_END = 1800.0


def theta(z):
    return THETA0 + D_THETA * np.minimum(z, INV_H) / INV_H


def pm_from_theta(th):
    """C/C_surface, linear in theta (exact mapping, see docstring)."""
    return 1.0 - (1.0 - PM_RATIO) * (th - THETA0) / D_THETA


def jet_penetration(T):
    M = T / RHO
    N = math.sqrt(9.81 / THETA0 * D_THETA / INV_H)
    return 2.9 * M ** 0.25 / math.sqrt(N)


def build(case, c):
    a = math.radians(WEDGE / 2)
    cs, sn = math.cos(a), math.sin(a)
    r0 = 0.02
    Rd, zh = c["Rd"], c["zh"]
    fm = c.get("refine", 1.0)           # mesh refinement factor (grid-independence study)
    dr = Rd / (4 * fm)                  # 4 cells across disk radius on the base mesh
    r_in = 4 * Rd
    n_in = int(math.ceil((r_in - r0) / dr))
    n_out = int(round(40 * fm))
    dz = dr
    z_lo, z_hi = max(zh - 4 * Rd - 2, 2.0), zh + 3 * Rd + 2
    n_z1 = int(round(30 * fm))
    n_z2 = int(math.ceil((z_hi - z_lo) / dz))
    n_z3 = int(round(60 * fm))
    nz = n_z1 + n_z2 + n_z3
    nr = n_in + n_out
    fr_in = (r_in - r0) / (R_CELL - r0)
    v = []
    for z in (0, H_DOM):
        v += [(r0 * cs, -r0 * sn, z), (R_CELL * cs, -R_CELL * sn, z), (R_CELL * cs, R_CELL * sn, z),
              (r0 * cs, r0 * sn, z)]
    vs = " ".join(f"({x:.6f} {y:.6f} {z:.6f})" for x, y, z in v)
    # bottom segment: expanding upward from fine near-ground cells; disk band uniform;
    # top segment expanding upward
    gz = (f"(({z_lo / H_DOM:.6f} {n_z1 / nz:.6f} 4) ({(z_hi - z_lo) / H_DOM:.6f} {n_z2 / nz:.6f} 1) "
          f"({(H_DOM - z_hi) / H_DOM:.6f} {n_z3 / nz:.6f} 6))")
    gr = f"(({fr_in:.6f} {n_in / nr:.6f} 1) ({1 - fr_in:.6f} {n_out / nr:.6f} 5))"
    write(case, "system/blockMeshDict", f"""
convertToMeters 1;
vertices ( {vs} );
blocks ( hex (0 1 2 3 4 5 6 7) ({nr} 1 {nz}) simpleGrading ({gr} 1 {gz}) );
boundary (
  ground {{ type wall; faces ( (0 3 2 1) ); }}
  top    {{ type wall; faces ( (4 5 6 7) ); }}
  outer  {{ type symmetry; faces ( (1 2 6 5) ); }}
  axis   {{ type symmetry; faces ( (0 4 7 3) ); }}
  front  {{ type wedge; faces ( (0 1 5 4) ); }}
  back   {{ type wedge; faces ( (3 7 6 2) ); }}
);
""")
    write(case, "system/topoSetDict", f"""
actions (
 {{ name disk; type cellSet; action new; source cylinderToCell;
   sourceInfo {{ p1 (0 0 {zh - dz}); p2 (0 0 {zh + dz}); radius {Rd}; }} }}
 {{ name diskZone; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set disk; }} }}
);
""")
    frac = WEDGE / 360.0
    fo = "" if c["T"] == 0 else f"""disk
{{
    type vectorSemiImplicitSource; active yes;
    vectorSemiImplicitSourceCoeffs
    {{
        timeStart 0; duration 1e9; selectionMode cellZone; cellZone diskZone; volumeMode absolute;
        injectionRateSuSp {{ U ((0 0 {-c['T'] / RHO * frac:.6g}) 0); }}
    }}
}}
"""
    write(case, "constant/fvOptions", fo)
    write(case, "constant/transportProperties",
          f"transportModel Newtonian;\nnu {NU};\nbeta {1 / THETA0:.6g};\nTRef {THETA0};\nPr 0.7;\nPrt {c.get('Prt', 0.7)};\n")
    write(case, "constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n",
          cls="uniformDimensionedVectorField")
    write(case, "constant/turbulenceProperties",
          "simulationType RAS;\nRAS { RASModel kEpsilon; turbulence on; printCoeffs off; }\n")
    write(case, "system/controlDict", f"""
application buoyantBoussinesqPimpleFoam; startFrom latestTime; stopAt endTime; endTime {T_END};
deltaT 0.01; writeControl adjustableRunTime; writeInterval 300; purgeWrite 0; writeFormat ascii;
writePrecision 8; runTimeModifiable false; adjustTimeStep yes; maxCo 1.5; maxDeltaT 1;
""")
    write(case, "system/fvSchemes", """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss limitedLinearV 1; div(phi,T) Gauss limitedLinear 1;
  div(phi,k) Gauss limitedLinear 1; div(phi,epsilon) Gauss limitedLinear 1;
  div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write(case, "system/fvSolution", """
solvers {
  p_rgh { solver GAMG; smoother DIC; tolerance 1e-7; relTol 0.01; }
  p_rghFinal { $p_rgh; relTol 0; }
  "(U|T|k|epsilon)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
  "(U|T|k|epsilon)Final" { $U; relTol 0; }
}
PIMPLE { momentumPredictor yes; nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0;
  pRefCell 0; pRefValue 0; }
""")
    wall = {"ground": None, "top": None}
    sym = {"outer": "        type symmetry;", "axis": "        type symmetry;",
           "front": "        type wedge;", "back": "        type wedge;"}

    def pf(walltxt):
        d = dict(sym)
        d["ground"] = walltxt
        d["top"] = walltxt
        return d
    write_field(case, "0", field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)",
                                 dict(sym, ground="        type noSlip;", top="        type slip;")))
    write_field(case, "0", field("p_rgh", "volScalarField", "[0 2 -2 0 0 0 0]", "0",
                                 pf("        type fixedFluxPressure; rho rhok; value uniform 0;")))
    write_field(case, "0", field("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0",
                                 pf("        type calculated; value uniform 0;")))
    write_field(case, "0", field("k", "volScalarField", "[0 2 -2 0 0 0 0]", "1e-5",
                                 dict(sym, ground="        type kqRWallFunction; value uniform 1e-5;",
                                      top="        type zeroGradient;")))
    write_field(case, "0", field("epsilon", "volScalarField", "[0 2 -3 0 0 0 0]", "1e-7",
                                 dict(sym, ground="        type epsilonWallFunction; value uniform 1e-7;",
                                      top="        type zeroGradient;")))
    write_field(case, "0", field("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "0",
                                 dict(sym, ground="        type nutkWallFunction; value uniform 0;",
                                      top="        type calculated; value uniform 0;")))
    write_field(case, "0", field("alphat", "volScalarField", "[0 2 -1 0 0 0 0]", "0",
                                 dict(sym, ground=f"        type alphatJayatillekeWallFunction; Prt {c.get('Prt', 0.7)}; value uniform 0;",
                                      top="        type calculated; value uniform 0;")))
    write_field(case, "0", field("T", "volScalarField", "[0 0 0 1 0 0 0]", f"{THETA0}",
                                 dict(sym, ground="        type zeroGradient;", top="        type zeroGradient;")))


def set_initial_theta(case):
    """Cell centres from a VTK export of the (pure hex) mesh -> nonuniform T."""
    g = vtk_zone(case, None, "VTK_mesh", fields=("T",))
    zc = np.asarray(g.cell_centers().points)[:, 2]
    th = theta(zc)
    p = os.path.join(case, "0", "T")
    s = open(p).read()
    body = "nonuniform List<scalar>\n%d\n(\n%s\n)" % (len(th), "\n".join(f"{x:.6f}" for x in th))
    s = s.replace(f"internalField   uniform {THETA0};", f"internalField   {body};")
    open(p, "w").write(s)
    return len(th)


def analyse(case):
    """PM in breathing / near-surface layers over the unit cell vs time."""
    out = []
    times = sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit()], key=float)
    import pyvista as pv
    import glob
    run(case, "foamToVTK -no-point-data -no-boundary -overwrite -name VTK_all -fields '(T U)'", "vtk.all")
    files = sorted(glob.glob(os.path.join(case, "VTK_all", "*", "internal.vtu")),
                   key=lambda f: int(f.split("_")[-1].split("/")[0]) if f.split("_")[-1].split("/")[0].isdigit() else 0)
    for f in files:
        g = pv.read(f).compute_cell_sizes(length=False, area=False, volume=True)
        cc = np.asarray(g.cell_centers().points)
        r = np.hypot(cc[:, 0], cc[:, 1])
        z = cc[:, 2]
        V = np.asarray(g.cell_data["Volume"])
        C = pm_from_theta(np.asarray(g.cell_data["T"]))
        rec = dict(file=f)
        for lab, zmax in (("z<2m", 2.0), ("z<10m", 10.0), ("z<30m", 30.0)):
            m = z < zmax
            rec[lab] = float((C[m] * V[m]).sum() / V[m].sum())
        out.append(rec)
    return out


def solve(name):
    c = CASES[name]
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    build(case, c)
    run(case, "blockMesh")
    run(case, "topoSet")
    n = set_initial_theta(case)
    run(case, "buoyantBoussinesqPimpleFoam", "solve")
    res = dict(case=name, cells=n, thrust_N=c["T"], hover_z=c["zh"], disk_R=c["Rd"],
               jet_penetration_theory_m=jet_penetration(c["T"]) if c["T"] else 0,
               series=analyse(case))
    json.dump(res, open(os.path.join(case, "result.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    names = sys.argv[1:] or list(CASES)
    os.makedirs(RUNS, exist_ok=True)
    with ThreadPoolExecutor(max_workers=int(os.environ.get('NPROC', '4'))) as ex:
        for r in ex.map(solve, names):
            print(r["case"], [(round(s["z<2m"], 3), round(s["z<10m"], 3)) for s in r["series"]])
