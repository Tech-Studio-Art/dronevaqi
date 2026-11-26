#!/usr/bin/env python3
"""
2-D steady RANS CFD for Solutions 2, 3a and 4a.

Common methodology (all literature-standard):
  * Steady RANS, standard k-epsilon, neutral atmospheric surface layer.
    Inlet profiles of U, k, epsilon from Richards & Hoxey (1993)
    (OpenFOAM atmBoundaryLayerInlet*), rough-wall function on the ground
    (nutkAtmRoughWallFunction), following the COST 732 / AIJ best-practice
    guidelines (Franke et al. 2007; Tominaga et al. 2008).
  * Pollutant: passive scalar, Sc_t = 0.7 (Tominaga & Stathopoulos 2007),
    solved as T of buoyantBoussinesqSimpleFoam with beta = 0 (pure tracer).
    Ambient (inflow) concentration normalised to 1. Fine (PM2.5) and coarse
    (2.5-10 um) fractions are separate runs where their sinks differ.
  * Porous devices (fence, fog net): Darcy-Forchheimer momentum sink with a
    pressure-loss coefficient K (dp = K 0.5 rho U^2), the standard porous-
    screen model in wind-engineering CFD (e.g. Santiago et al. 2007).
  * Collection inside a device: first-order sink:
      ESP fence  : Deutsch-Anderson, lambda = w * a  (w = effective migration
                   velocity, a = collecting area per volume) (White 1963).
      Fog net    : Langmuir/Rivera impaction on filaments, captured droplet
                   fraction per pass eta_d = Stk/(Stk + pi/2) (Rivera 2011),
                   lambda = -ln(1-eta_d) U/L. PM rides in activated fog drops.
  * Clean-air shower: supply jet of filtered air from the shelter roof.

Settling of PM10 (<= 5 mm/s) is neglected: residence times in these domains
are tens of seconds.
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

RHO, NU, MU = 1.25, 1.42e-5, 1.78e-5
RUNS = os.path.join(HERE, "runs", "2d")
Z0 = 0.05


# ------------------------------------------------------------- lattice mesher
MESH_F = 1.0   # global mesh refinement factor (grid-independence study); 1.0 = base mesh


def spacing(a, b, h):
    return max(1, int(math.ceil(abs(b - a) * MESH_F / h)))


def lattice_mesh(case, xb, zb, hx, hz, solid=(), special=None, depth=1.0):
    """Multi-block 2-D blockMesh on a lattice of x/z breakpoints.
    hx/hz: target cell size per interval. solid: set of (i, j) blocks removed.
    special: f(i, j, side) -> patch name or None for faces against solids/domain.
    Domain patches: inlet (x min), outlet (x max), ground (z min), top (z max)."""
    nx, nz = len(xb) - 1, len(zb) - 1
    vid = {}
    verts = []
    for k, y in enumerate((0.0, depth)):
        for j, z in enumerate(zb):
            for i, x in enumerate(xb):
                vid[(i, j, k)] = len(verts)
                verts.append((x, y, z))
    blocks, faces = [], {}

    def add(p, f):
        faces.setdefault(p, []).append(f)
    for j in range(nz):
        for i in range(nx):
            if (i, j) in solid:
                continue
            v = [vid[(i, j, 0)], vid[(i + 1, j, 0)], vid[(i + 1, j, 1)], vid[(i, j, 1)],
                 vid[(i, j + 1, 0)], vid[(i + 1, j + 1, 0)], vid[(i + 1, j + 1, 1)], vid[(i, j + 1, 1)]]
            blocks.append(f"hex ({' '.join(map(str, v))}) ({spacing(xb[i], xb[i + 1], hx[i])} 1 "
                          f"{spacing(zb[j], zb[j + 1], hz[j])}) simpleGrading (1 1 1)")
            # faces: x-, x+, z-, z+
            nb = {"xm": (i - 1, j), "xp": (i + 1, j), "zm": (i, j - 1), "zp": (i, j + 1)}
            fdef = {"xm": (v[0], v[4], v[7], v[3]), "xp": (v[1], v[2], v[6], v[5]),
                    "zm": (v[0], v[3], v[2], v[1]), "zp": (v[4], v[5], v[6], v[7])}
            for side, (ii, jj) in nb.items():
                outside = ii < 0 or ii >= nx or jj < 0 or jj >= nz
                if outside or (ii, jj) in solid:
                    p = special(i, j, side) if special else None
                    if p is None:
                        if outside:
                            p = {"xm": "inlet", "xp": "outlet", "zm": "ground", "zp": "top"}[side]
                        else:
                            p = "walls"
                    add(p, fdef[side])
            add("frontAndBack", (v[0], v[1], v[5], v[4]))
            add("frontAndBack", (v[3], v[7], v[6], v[2]))
    types = {"inlet": "patch", "outlet": "patch", "top": "patch", "ground": "wall", "walls": "wall",
             "frontAndBack": "empty"}
    bnd = "\n".join(f"  {p} {{ type {types.get(p, 'patch' if p != 'walls' else 'wall')}; faces ( "
                    + " ".join("(" + " ".join(map(str, f)) + ")" for f in fl) + " ); }"
                    for p, fl in faces.items())
    write(case, "system/blockMeshDict", "convertToMeters 1;\nvertices (\n"
          + "\n".join(f"({x} {y} {z})" for x, y, z in verts) + "\n);\nblocks (\n" + "\n".join(blocks)
          + "\n);\nedges ();\nboundary (\n" + bnd + "\n);\nmergePatchPairs ();\n")
    return list(faces)


# ------------------------------------------------------------- common setup
def common(case, patches, Uref, fvoptions="", supply=None, wall_patches=("walls",), n_iter=3000):
    write(case, "constant/fvOptions", fvoptions)
    write(case, "constant/transportProperties",
          f"transportModel Newtonian;\nnu {NU};\nbeta 0;\nTRef 0;\nPr 1.0;\nPrt 0.7;\n")
    write(case, "constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n",
          cls="uniformDimensionedVectorField")
    write(case, "constant/turbulenceProperties",
          "simulationType RAS;\nRAS { RASModel kEpsilon; turbulence on; printCoeffs off; }\n")
    write(case, "system/controlDict", f"""
application buoyantBoussinesqSimpleFoam; startFrom latestTime; stopAt endTime; endTime {n_iter};
deltaT 1; writeControl timeStep; writeInterval {n_iter}; purgeWrite 1; writeFormat ascii;
writePrecision 8; runTimeModifiable false;
""")
    write(case, "system/fvSchemes", """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U);
  div(phi,k) bounded Gauss upwind; div(phi,epsilon) bounded Gauss upwind;
  div(phi,T) bounded Gauss limitedLinear 1; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write(case, "system/fvSolution", """
solvers {
  p_rgh { solver GAMG; smoother GaussSeidel; tolerance 1e-8; relTol 0.05; }
  "(U|k|epsilon)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol 0.1; }
  T { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0.05; maxIter 200; }
}
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; pRefCell 0; pRefValue 0;
  residualControl { p_rgh 1e-6; U 1e-7; T 1e-8; "(k|epsilon)" 1e-7; } }
relaxationFactors { fields { p_rgh 0.3; } equations { U 0.7; k 0.5; epsilon 0.5; T 0.7; } }
""")
    abl = (f"flowDir (1 0 0); zDir (0 0 1); Uref {Uref}; Zref 10; z0 uniform {Z0}; zGround uniform 0; "
           f"kappa 0.41; Cmu 0.09;")
    U, P, K, E, NUT, AT, T = ({} for _ in range(7))
    for p in patches:
        if p == "inlet":
            U[p] = f"type atmBoundaryLayerInletVelocity; {abl} value uniform ({Uref} 0 0);"
            K[p] = f"type atmBoundaryLayerInletK; {abl} value uniform 0.01;"
            E[p] = f"type atmBoundaryLayerInletEpsilon; {abl} value uniform 0.001;"
            P[p] = "type fixedFluxPressure; value uniform 0;"
            T[p] = "type fixedValue; value uniform 1;"
            NUT[p] = AT[p] = "type calculated; value uniform 0;"
        elif p == "outlet":
            U[p] = "type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0);"
            P[p] = "type fixedValue; value uniform 0;"
            K[p] = "type inletOutlet; inletValue uniform 1e-4; value uniform 1e-4;"
            E[p] = "type inletOutlet; inletValue uniform 1e-5; value uniform 1e-5;"
            T[p] = "type inletOutlet; inletValue uniform 1; value uniform 1;"
            NUT[p] = AT[p] = "type calculated; value uniform 0;"
        elif p in ("top", "sides"):
            U[p] = "type slip;"
            P[p] = K[p] = E[p] = T[p] = "type zeroGradient;"
            NUT[p] = AT[p] = "type calculated; value uniform 0;"
        elif p == "frontAndBack":
            U[p] = P[p] = K[p] = E[p] = T[p] = NUT[p] = AT[p] = "type empty;"
        elif supply and p in supply:
            w, c = supply[p]
            U[p] = f"type fixedValue; value uniform ({w[0]} 0 {w[1]});"
            P[p] = "type fixedFluxPressure; value uniform 0;"
            K[p] = "type fixedValue; value uniform 0.001;"
            E[p] = "type fixedValue; value uniform 0.0005;"
            T[p] = f"type fixedValue; value uniform {c};"
            NUT[p] = AT[p] = "type calculated; value uniform 0;"
        else:  # ground / walls
            U[p] = "type noSlip;"
            P[p] = "type fixedFluxPressure; value uniform 0;"
            K[p] = "type kqRWallFunction; value uniform 0.01;"
            E[p] = "type epsilonWallFunction; value uniform 0.001;"
            T[p] = "type zeroGradient;"
            NUT[p] = (f"type nutkAtmRoughWallFunction; z0 uniform {Z0}; value uniform 0;" if p == "ground"
                      else "type nutkWallFunction; value uniform 0;")
            AT[p] = "type alphatJayatillekeWallFunction; Prt 0.7; value uniform 0;"
    for name, cls, dim, init, d in (
            ("U", "volVectorField", "[0 1 -1 0 0 0 0]", f"({Uref} 0 0)", U),
            ("p_rgh", "volScalarField", "[0 2 -2 0 0 0 0]", "0", P),
            ("k", "volScalarField", "[0 2 -2 0 0 0 0]", "0.01", K),
            ("epsilon", "volScalarField", "[0 2 -3 0 0 0 0]", "0.001", E),
            ("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "0", NUT),
            ("alphat", "volScalarField", "[0 2 -1 0 0 0 0]", "0", AT),
            ("T", "volScalarField", "[0 0 0 1 0 0 0]", "1", T)):
        write_field(case, "0", field(name, cls, dim, init, {p: "        " + b for p, b in d.items()}))
    pp = {p: "        type calculated; value uniform 0;" if p != "frontAndBack" else "        type empty;"
          for p in patches}
    write_field(case, "0", field("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0", pp))


def porous_fv(zone, K, L, lam=0.0):
    f = K / L
    s = f"""porous
{{
    type explicitPorositySource; active yes;
    explicitPorositySourceCoeffs
    {{
        selectionMode cellZone; cellZone {zone}; type DarcyForchheimer;
        DarcyForchheimerCoeffs
        {{
            d d [0 -2 0 0 0 0 0] (0 0 0);
            f f [0 -1 0 0 0 0 0] ({f:.6g} {100 * f:.6g} {100 * f:.6g});
            coordinateSystem {{ origin (0 0 0); e1 (1 0 0); e2 (0 1 0); }}
        }}
    }}
}}
"""
    if lam > 0:
        s += f"""capture
{{
    type scalarSemiImplicitSource; active yes;
    scalarSemiImplicitSourceCoeffs
    {{
        timeStart 0; duration 1e9; selectionMode cellZone; cellZone {zone}; volumeMode specific;
        injectionRateSuSp {{ T (0 {-lam:.6g}); }}
    }}
}}
"""
    return s


def box_zone(case, name, p0, p1):
    write(case, "system/topoSetDict", f"""
actions (
 {{ name {name}; type cellSet; action new; source boxToCell; sourceInfo {{ box ({p0[0]} -1 {p0[1]}) ({p1[0]} 2 {p1[1]}); }} }}
 {{ name {name}Zone; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set {name}; }} }}
);
""")


def load(case):
    g = vtk_zone(case, None, "VTK_full", fields=("U", "T"))
    cc = np.asarray(g.cell_centers().points)
    return dict(x=cc[:, 0], z=cc[:, 2], V=np.asarray(g.cell_data["Volume"]),
                T=np.asarray(g.cell_data["T"]), U=np.asarray(g.cell_data["U"]))


def region_mean(d, x0, x1, z0, z1):
    m = (d["x"] > x0) & (d["x"] < x1) & (d["z"] > z0) & (d["z"] < z1)
    return float((d["T"][m] * d["V"][m]).sum() / d["V"][m].sum())


# ============================================================ Solution 2
def shelter(Uref, w_supply, wind_from="front", c_supply=0.05):
    """Bus shelter cross-section: depth 2.5 m, roof at 2.5-2.6 m, back wall 0.1 m.
    wind_from='front': wind blows into the open side; 'back': over the back wall."""
    name = f"shelter_U{Uref}_w{w_supply}_{wind_from}"
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    if wind_from == "back":   # back wall upwind (x 0-0.1), opening faces downwind (+x)
        xb = [-40, -6, -1, 0, 0.1, 0.3, 2.3, 2.6, 4, 10, 80]
        wall_i, roof_i, sup_i = {3}, {3, 4, 5, 6}, {5}
        inside = (0.1, 2.6)
    else:                     # opening faces upwind, back wall at x 2.5-2.6
        xb = [-40, -10, -4, 0, 0.3, 2.3, 2.5, 2.6, 4, 10, 80]
        wall_i, roof_i, sup_i = {6}, {3, 4, 5, 6}, {4}
        inside = (0.0, 2.5)
    zb = [0, 0.5, 1.8, 2.5, 2.6, 4, 10, 40]
    hx = [{0: 3, 1: 1, 2: 0.25, 9: 1.5, 8: 0.6}.get(i, 0.05) for i in range(len(xb) - 1)]
    hz = [0.05, 0.05, 0.05, 0.05, 0.1, 0.4, 1.5]
    solid = {(i, 3) for i in roof_i} | {(i, j) for i in wall_i for j in (0, 1, 2)}

    def special(i, j, side):
        if side == "zp" and j == 2 and i in sup_i:
            return "supply"
        return None
    patches = lattice_mesh(case, xb, zb, hx, hz, solid, special)
    supply = {"supply": ((0, -w_supply), c_supply)} if w_supply > 0 else None
    if w_supply == 0:
        # closed roof panel instead of supply
        txt = open(os.path.join(case, "system/blockMeshDict")).read().replace("supply { type patch;", "supply { type wall;")
        open(os.path.join(case, "system/blockMeshDict"), "w").write(txt)
    common(case, patches, Uref, supply=supply)
    run(case, "blockMesh")
    run(case, "buoyantBoussinesqSimpleFoam", "solve")
    d = load(case)
    bz = region_mean(d, inside[0], inside[1], 0.5, 1.8)
    sup_len = 2.0
    Q = w_supply * sup_len           # m3/s per metre of shelter length
    # blower power: filter (HEPA-class 250 Pa) + ducting/outlet 100 Pa, fan eff 0.5
    P = Q * 350 / 0.5
    return dict(case=name, Uref=Uref, w_supply=w_supply, wind_from=wind_from,
                breathing_zone_C=bz, reduction_pct=100 * (1 - bz), Q_per_m=Q, fan_W_per_m=P)


# ============================================================ Solution 3a
ESP = dict(H=6.0, L=0.5, K=1.5, a=20.0, w25=0.04, wc=0.12)


def esp_fence(Uref, size):
    name = f"espfence_U{Uref}_{size}"
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    H, L = ESP["H"], ESP["L"]
    xb = [-60, -10, -2, 0, L, 3, 15, 60, 150]
    zb = [0, 2, H, H + 2, 15, 60]
    hx = [3, 1, 0.25, 0.1, 0.25, 0.5, 1.5, 3]
    hz = [0.1, 0.2, 0.25, 0.75, 3]
    patches = lattice_mesh(case, xb, zb, hx, hz)
    lam = ESP["a"] * (ESP["w25"] if size == "PM25" else ESP["wc"])
    box_zone(case, "fence", (0, 0), (L, H))
    common(case, patches, Uref, fvoptions=porous_fv("fenceZone", ESP["K"], L, lam))
    run(case, "blockMesh")
    run(case, "topoSet")
    run(case, "buoyantBoussinesqSimpleFoam", "solve")
    d = load(case)
    m = (d["x"] > 0) & (d["x"] < L) & (d["z"] < H)
    # mass removal per metre of fence (m3/s ambient-equivalent) = lambda * int(C dV) / depth
    cadr = lam * float((d["T"][m] * d["V"][m]).sum())
    u_through = float((d["U"][m, 0] * d["V"][m]).sum() / d["V"][m].sum())
    prof = {f"x={x}m": region_mean(d, x - 1.0, x + 1.0, 1.0, 2.0) for x in (2, 5, 10, 20, 40, 80)}
    return dict(case=name, Uref=Uref, size=size, lam=lam, u_through=u_through,
                frac_of_free_wind=u_through / (Uref * math.log((H / 2 + Z0) / Z0) / math.log((10 + Z0) / Z0)),
                single_pass_eff=1 - math.exp(-lam * L / max(u_through, 1e-6)),
                CADR_per_m=cadr, pedestrian_C=prof)


# ============================================================ Solution 4a
FOG = dict(H=4.0, zb=2.0, L=0.2, shade=0.35, d_fil=1.0e-3, d_drop=10e-6, LWC=0.2e-3, f_act=0.5)


def fog_C0(s):
    """Pressure-drop coefficient of a Raschel mesh (Idel'chik form, used by Rivera 2011)."""
    return 1.3 * s + (s / (1 - s)) ** 2


def fog_net(Uref):
    name = f"fognet_U{Uref}"
    case = os.path.join(RUNS, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    H, z0n, L = FOG["H"], FOG["zb"], FOG["L"]
    xb = [-60, -10, -2, 0, L, 3, 15, 60, 150]
    zb = [0, z0n, z0n + H, z0n + H + 2, 15, 60]
    hx = [3, 1, 0.25, 0.05, 0.25, 0.5, 1.5, 3]
    hz = [0.2, 0.2, 0.25, 0.75, 3]
    patches = lattice_mesh(case, xb, zb, hx, hz)
    box_zone(case, "net", (0, z0n), (L, z0n + H))
    K = fog_C0(FOG["shade"])
    # pass 1: flow only, to get the through-velocity
    common(case, patches, Uref, fvoptions=porous_fv("netZone", K, L), n_iter=2500)
    run(case, "blockMesh")
    run(case, "topoSet")
    run(case, "buoyantBoussinesqSimpleFoam", "solve.flow")
    d = load(case)
    m = (d["x"] > 0) & (d["x"] < L) & (d["z"] > z0n) & (d["z"] < z0n + H)
    u = float((d["U"][m, 0] * d["V"][m]).sum() / d["V"][m].sum())
    Stk = 1000 * FOG["d_drop"] ** 2 * u / (9 * MU * FOG["d_fil"])
    eta_d = Stk / (Stk + math.pi / 2)
    # droplets meet filaments only on the solid (shade) fraction of the mesh
    eta_pass = FOG["shade"] * eta_d
    lam = -math.log(1 - eta_pass) * u / L
    common(case, patches, Uref, fvoptions=porous_fv("netZone", K, L, lam), n_iter=4000)
    run(case, "buoyantBoussinesqSimpleFoam", "solve.drops")
    d = load(case)
    cap = lam * float((d["T"][m] * d["V"][m]).sum())     # m3/s of foggy air "emptied" per m of net
    u10_free = Uref * math.log((z0n + H / 2 + Z0) / Z0) / math.log((10 + Z0) / Z0)
    eta_ae = u / u10_free * 1.0
    water_L_per_m2_day = cap * FOG["LWC"] / H * 86400  # kg = L
    pm_cadr = FOG["f_act"] * cap
    prof = {f"x={x}m": 1 - FOG["f_act"] * (1 - region_mean(d, x - 1.0, x + 1.0, 1.0, 2.0))
            for x in (2, 5, 10, 20, 40)}
    return dict(case=name, Uref=Uref, C0=K, u_through=u, frac_free_wind=eta_ae, Stk=Stk, eta_d=eta_d,
                eta_total=eta_ae * eta_pass, water_L_m2_day=water_L_per_m2_day,
                PM_CADR_per_m_net=pm_cadr, PM_CADR_per_m2_net=pm_cadr / H,
                pedestrian_PM_C=prof)


JOBS = ([("shelter", (u, w, wf)) for u in (0.5, 1.0, 2.0) for w in (0.0, 0.25, 0.5, 1.0)
         for wf in ("front", "back")]
        + [("esp", (u, s)) for u in (1.0, 2.0) for s in ("PM25", "coarse")]
        + [("fog", (u,)) for u in (1.0, 2.0)])


def dispatch(job):
    kind, args = job
    try:
        return {"shelter": shelter, "esp": esp_fence, "fog": fog_net}[kind](*args)
    except Exception as e:  # keep the sweep going, record failure
        return dict(case=f"{kind}{args}", error=str(e))


if __name__ == "__main__":
    os.makedirs(RUNS, exist_ok=True)
    sel = sys.argv[1:]
    jobs = [j for j in JOBS if not sel or j[0] in sel]
    nproc = int(os.environ.get("NPROC", "4"))
    with ThreadPoolExecutor(max_workers=nproc) as ex:
        res = list(ex.map(dispatch, jobs))
    tag = "_".join(sel) or "all"
    json.dump(res, open(os.path.join(RUNS, f"summary_{tag}.json"), "w"), indent=1, default=float)
    for r in res:
        print(r)
