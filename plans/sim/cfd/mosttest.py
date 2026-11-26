#!/usr/bin/env python3
"""
Validation (step 1b-ii): stable atmospheric surface layer, single column.

Purpose: test whether a turbulence model reproduces the measured flux-profile relations of
the stable surface layer (Monin-Obukhov similarity). Measured relations (field data):
  Businger et al. (1971):   phi_m = 1 + 4.7 zeta,  phi_h = 0.74 + 4.7 zeta
  Dyer (1974):              phi_m = phi_h = 1 + 5 zeta
  Hogstrom (1988):          phi_m = 1 + 6 zeta,    phi_h = 0.95 + 7.8 zeta
with zeta = z/L. A model without buoyancy effects on turbulence gives phi_m ~ 1 at all zeta.

Set-up (GABLS1-type forcing, no Coriolis): 1-D column 400 m high, z0 = 0.1 m, steady driving
force for u* ~ 0.3 m/s, initial theta 265 K below 100 m and +0.01 K/m above, surface cooled
at 0.25 K/h for 9 h. Local similarity (Nieuwstadt 1984) gives phi_m and zeta from local
fluxes: u*_l^2 = nu_t dU/dz, w'theta' = -(alpha_t) dtheta/dz, Lambda = -u*_l^3 theta0/(kappa g w'theta').

Models:
  std  : buoyantBoussinesqPimpleFoam, standard k-epsilon (no buoyancy terms in k, epsilon)
  buoy : buoyantPimpleFoam, buoyantKEpsilon (buoyancy production/destruction in k and epsilon),
         Boussinesq equation of state and a large specific heat, so that T is the potential
         temperature (adiabatic lapse rate 0.1 K/km, negligible): the potential-temperature
         form used in atmospheric boundary-layer CFD.
"""
import json
import math
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from foam import write, run, field, write_field  # noqa: E402

RUNS = os.path.join(HERE, "runs", "validation", os.environ.get("MOST_DIR", "most"))
H, TH0 = 400.0, 265.0
N = int(os.environ.get("MOST_N", "128"))
Z0 = float(os.environ.get("MOST_Z0", "0.1"))
GRADING = float(os.environ.get("MOST_GRADING", "8"))
USTAR = 0.3
FX = USTAR ** 2 / H              # m/s2 driving force (stress balance u*^2 = fx H)
RHO0, NU, KAPPA, G = 1.3, 1.4e-5, 0.41, 9.81
T_SPIN = 24 * 3600.0            # neutral spin-up before the surface cooling starts
T_END = T_SPIN + float(os.environ.get("MOST_COOL_H", "9")) * 3600.0
CP_BIG = 1.0e5


def theta0(z):
    return TH0 + np.maximum(z - 100.0, 0) * 0.01


def mesh(case):
    write(case, "system/blockMeshDict", f"""
convertToMeters 1;
vertices ( (0 0 0) (300 0 0) (300 1 0) (0 1 0) (0 0 {H}) (300 0 {H}) (300 1 {H}) (0 1 {H}) );
blocks ( hex (0 1 2 3 4 5 6 7) (3 1 {N}) simpleGrading (1 1 {GRADING}) );
boundary (
  ground {{ type wall; faces ( (0 3 2 1) ); }}
  top    {{ type wall; faces ( (4 5 6 7) ); }}
  left   {{ type cyclic; neighbourPatch right; faces ( (0 4 7 3) ); }}
  right  {{ type cyclic; neighbourPatch left; faces ( (1 2 6 5) ); }}
  frontAndBack {{ type empty; faces ( (0 1 5 4) (3 7 6 2) ); }}
);
""")


def common_dicts(case, model, app):
    write(case, "constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n", cls="uniformDimensionedVectorField")
    write(case, "constant/turbulenceProperties",
          f"simulationType RAS;\nRAS {{ RASModel {model}; turbulence on; printCoeffs on; "
          f"buoyantKEpsilonCoeffs {{ Cg {os.environ.get('MOST_CG', '1.0')}; }} }}\n")
    write(case, "system/controlDict", f"""
application {app}; startFrom latestTime; stopAt endTime; endTime {T_END};
deltaT 1; writeControl adjustableRunTime; writeInterval 3600; purgeWrite 0; writeFormat ascii;
writePrecision 10; runTimeModifiable false; adjustTimeStep yes; maxCo 1; maxDeltaT 5;
libs ("libatmosphericModels.so");
""")
    write(case, "system/fvSchemes", """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default Gauss upwind; div(phi,U) Gauss upwind; div(phi,T) Gauss upwind; div(phi,h) Gauss upwind;
  div(phi,K) Gauss linear; div(phi,k) Gauss upwind; div(phi,epsilon) Gauss upwind;
  div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write(case, "system/fvSolution", """
solvers {
  "(p_rgh|p_rghFinal)" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; maxIter 2000; }
  "(U|T|h|k|epsilon)(|Final)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }
  "rho.*" { solver diagonal; }
}
PIMPLE { momentumPredictor yes; nOuterCorrectors 2; nCorrectors 2; nNonOrthogonalCorrectors 0;
  pRefCell 0; pRefValue 0; }
relaxationFactors { equations { ".*" 1; } }
""")


def fields(case, comp):
    z = (np.arange(N) + 0.5) / N  # placeholder; real centres set after blockMesh
    U = "        type noSlip;"
    walls = dict(ground=None, top=None)
    emp = {"left": "        type cyclic;", "right": "        type cyclic;", "frontAndBack": "        type empty;"}
    write_field(case, "0", field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(5 0 0)",
                                 dict(emp, ground="        type noSlip;", top="        type slip;")))
    pdim = "[1 -1 -2 0 0 0 0]" if comp else "[0 2 -2 0 0 0 0]"
    pval = "1e5" if comp else "0"
    write_field(case, "0", field("p_rgh", "volScalarField", pdim, pval,
                                 dict(emp, ground="        type fixedFluxPressure; value uniform %s;" % pval,
                                      top="        type fixedFluxPressure; value uniform %s;" % pval)))
    write_field(case, "0", field("p", "volScalarField", pdim, pval,
                                 dict(emp, ground="        type calculated; value uniform %s;" % pval,
                                      top="        type calculated; value uniform %s;" % pval)))
    k0, e0 = USTAR ** 2 / math.sqrt(0.09), USTAR ** 3 / (KAPPA * 10)
    write_field(case, "0", field("k", "volScalarField", "[0 2 -2 0 0 0 0]", k0,
                                 dict(emp, ground=f"        type kqRWallFunction; value uniform {k0};", top="        type zeroGradient;")))
    write_field(case, "0", field("epsilon", "volScalarField", "[0 2 -3 0 0 0 0]", e0,
                                 dict(emp, ground=f"        type epsilonWallFunction; value uniform {e0};", top="        type zeroGradient;")))
    write_field(case, "0", field("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "0",
                                 dict(emp, ground=f"        type nutkAtmRoughWallFunction; z0 uniform {Z0}; value uniform 0;",
                                      top="        type calculated; value uniform 0;")))
    if comp:
        at = "        type compressible::alphatJayatillekeWallFunction; Prt 0.85; value uniform 0;"
        adim = "[1 -1 -1 0 0 0 0]"
    else:
        at = "        type alphatJayatillekeWallFunction; Prt 0.85; value uniform 0;"
        adim = "[0 2 -1 0 0 0 0]"
    write_field(case, "0", field("alphat", "volScalarField", adim, "0",
                                 dict(emp, ground=at, top="        type calculated; value uniform 0;")))
    write_field(case, "0", field("T", "volScalarField", "[0 0 0 1 0 0 0]", TH0,
                                 dict(emp, ground=f"        type uniformFixedValue; uniformValue table ( (0 {TH0}) ({T_SPIN} {TH0}) ({T_END} {TH0 - 0.25 * (T_END - T_SPIN) / 3600}) ); value uniform {TH0};",
                                      top="        type zeroGradient;")))


def set_profiles(case):
    """Nonuniform initial T (GABLS1 profile) and log-law U from the actual cell centres."""
    pts = open(os.path.join(case, "constant/polyMesh/points")).read()
    import re
    zs = sorted(set(round(float(m.split()[2]), 6) for m in re.findall(r"\(([-\d.e+ ]+)\)", pts) if len(m.split()) == 3))
    zc1 = np.array([(zs[i] + zs[i + 1]) / 2 for i in range(len(zs) - 1)])
    zc = np.repeat(zc1, 3)          # 3 cells per level (blockMesh orders x fastest)
    th = theta0(zc)
    u = np.maximum(USTAR / KAPPA * np.log((zc + Z0) / Z0), 0.1)
    for name, vals, fmt in (("T", th, "{:.6f}"), ("U", u, "({:.6f} 0 0)")):
        p = os.path.join(case, "0", name)
        s = open(p).read()
        kind = "List<scalar>" if name == "T" else "List<vector>"
        body = f"nonuniform {kind}\n{len(vals)}\n(\n" + "\n".join(fmt.format(v) for v in vals) + "\n)"
        s = re.sub(r"internalField\s+uniform [^;]+;", f"internalField   {body};", s, count=1)
        open(p, "w").write(s)
    return zc


def build(variant):
    case = os.path.join(RUNS, variant)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    mesh(case)
    comp = variant == "buoy"
    if comp:
        common_dicts(case, "buoyantKEpsilon", "buoyantPimpleFoam")
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
        # semi-implicit driving force in compressible form (rho*fx per unit volume)
        write(case, "constant/fvOptions", f"""
drive {{ type vectorSemiImplicitSource; active yes;
  vectorSemiImplicitSourceCoeffs {{ timeStart 0; duration 1e9; selectionMode all; volumeMode specific;
    injectionRateSuSp {{ U (({RHO0 * FX:.8g} 0 0) 0); }} }} }}
""")
    else:
        common_dicts(case, "kEpsilon", "buoyantBoussinesqPimpleFoam")
        write(case, "constant/transportProperties",
              f"transportModel Newtonian;\nnu {NU};\nbeta {1 / TH0:.8g};\nTRef {TH0};\nPr 0.7;\nPrt 0.85;\n")
        write(case, "constant/fvOptions", f"""
drive {{ type vectorSemiImplicitSource; active yes;
  vectorSemiImplicitSourceCoeffs {{ timeStart 0; duration 1e9; selectionMode all; volumeMode specific;
    injectionRateSuSp {{ U (({FX:.8g} 0 0) 0); }} }} }}
""")
    fields(case, comp)
    run(case, "blockMesh")
    set_profiles(case)
    return case


def analyse(case, comp):
    from foam import vtk_zone
    import glob
    import pyvista as pv
    run(case, f"foamToVTK -time '{T_SPIN + 3600:g}:' -no-point-data -no-boundary -overwrite -name VTK_all "
              f"-fields '(U T nut alphat)'", "vtk.all")
    times = sorted(float(d) for d in os.listdir(case) if d.replace(".", "").replace("e", "").isdigit() and float(d) > T_SPIN + 1)
    out = []
    files = sorted(glob.glob(os.path.join(case, "VTK_all", "*", "internal.vtu")),
                   key=lambda f: int(f.split("_")[-1].split("/")[0]) if f.split("_")[-1].split("/")[0].isdigit() else 0)
    for f, tt in zip(files, times):
        g = pv.read(f)
        zr = np.round(np.asarray(g.cell_centers().points)[:, 2], 6)
        z, inv = np.unique(zr, return_inverse=True)          # 3 identical cells per level -> mean
        avg = lambda a: np.bincount(inv, weights=a) / np.bincount(inv)
        U = avg(np.asarray(g.cell_data["U"])[:, 0])
        T = avg(np.asarray(g.cell_data["T"]))
        nut = avg(np.asarray(g.cell_data["nut"]))
        at = avg(np.asarray(g.cell_data["alphat"])) / (RHO0 if comp else 1.0)
        dU = np.gradient(U, z)
        dT = np.gradient(T, z)
        ust = np.sqrt(np.maximum(nut * dU, 1e-12))
        wth = -at * dT
        with np.errstate(divide="ignore", invalid="ignore"):
            Lam = -ust ** 3 * TH0 / (KAPPA * G * wth)
            zeta = z / Lam
            phim = KAPPA * z * dU / ust
            phih = KAPPA * z * ust * dT / np.where(wth != 0, -wth, np.nan)
        sel = (z > 3) & (z < 60) & np.isfinite(zeta) & (zeta > 0)
        out.append(dict(file=f, hours_cooling=(tt - T_SPIN) / 3600, z=z[sel].tolist(), zeta=zeta[sel].tolist(), phim=phim[sel].tolist(),
                        phih=phih[sel].tolist(), ustar_surface=float(ust[1]),
                        sbl_height=float(z[np.argmax(ust < 0.05 * ust[1])]) if (ust < 0.05 * ust[1]).any() else None))
    return out


def main(variants):
    res = {}
    for v in variants:
        case = build(v)
        app = "buoyantPimpleFoam" if v == "buoy" else "buoyantBoussinesqPimpleFoam"
        run(case, app, "solve")
        res[v] = analyse(case, v == "buoy")
    os.makedirs(RUNS, exist_ok=True)
    old = json.load(open(os.path.join(RUNS, "summary.json"))) if os.path.exists(os.path.join(RUNS, "summary.json")) else {}
    old.update(res)
    json.dump(old, open(os.path.join(RUNS, "summary.json"), "w"))
    for v, r in res.items():
        last = r[-1]
        print(v, "u*_s", round(last["ustar_surface"], 3), "h", last["sbl_height"],
              [(round(a, 2), round(b, 2)) for a, b in zip(last["zeta"][::6], last["phim"][::6])])


if __name__ == "__main__":
    main(sys.argv[1:] or ["std", "buoy"])


def _read_internal(path, vec=False):
    import re
    s = open(path).read()
    seg = s[s.index("internalField"):]
    if "nonuniform" not in seg.split(";")[0]:
        v = seg.split("uniform")[1].split(";")[0].strip()
        return None, (float(v.strip("()").split()[0]) if vec else float(v))
    n = int(re.search(r"\n(\d+)\n\(", seg).group(1))
    body = seg[seg.index("(") + 1:]
    if vec:
        vals = re.findall(r"\(([^)]+)\)", body)[:n]
        return np.array([float(v.split()[0]) for v in vals]), None
    return np.array([float(x) for x in body.split(")")[0].split()[:n]]), None


def analyse_budget(case, comp):
    """Fluxes from the solver's own budgets (exact for a horizontally homogeneous column):
       tau(z)   = integral_z^H (fx - dU/dt) dz'      (stress-free top)
       wth(z)   = integral_z^H d(theta)/dt dz'      (no heat flux at the top)
    Gradients from the hour-mean profiles (3-point smoothing)."""
    import re
    pts = open(os.path.join(case, "constant/polyMesh/points")).read()
    zs = sorted(set(round(float(p.split()[2]), 6) for p in re.findall(r"\(([-\d.e+ ]+)\)", pts) if len(p.split()) == 3))
    zf = np.array(zs)
    z = 0.5 * (zf[1:] + zf[:-1])
    dz = np.diff(zf)
    times = sorted((float(d), d) for d in os.listdir(case) if re.fullmatch(r"[0-9.e+]+", d))
    times = [t for t in times if t[0] >= T_SPIN]
    sm = lambda a: np.convolve(np.pad(a, 1, mode="edge"), np.ones(3) / 3, mode="valid")
    out = []
    for (t1, d1), (t2, d2) in zip(times[:-1], times[1:]):
        U1 = _read_internal(os.path.join(case, d1, "U"), True)[0].reshape(-1, 3).mean(1)
        U2 = _read_internal(os.path.join(case, d2, "U"), True)[0].reshape(-1, 3).mean(1)
        T1 = _read_internal(os.path.join(case, d1, "T"))[0].reshape(-1, 3).mean(1)
        T2 = _read_internal(os.path.join(case, d2, "T"))[0].reshape(-1, 3).mean(1)
        dt = t2 - t1
        tend_u, tend_t = (U2 - U1) / dt, (T2 - T1) / dt
        # integrate from the top: value at cell centre = sum over cells above + half own cell
        cum = lambda f: np.cumsum((f * dz)[::-1])[::-1] - 0.5 * f * dz
        tau = cum(FX - tend_u)
        wth = cum(tend_t)
        Um, Tm = sm(0.5 * (U1 + U2)), sm(0.5 * (T1 + T2))
        dU, dT = np.gradient(Um, z), np.gradient(Tm, z)
        ust = np.sqrt(np.maximum(tau, 1e-12))
        with np.errstate(divide="ignore", invalid="ignore"):
            Lam = -ust ** 3 * TH0 / (KAPPA * G * wth)
            zeta = z / Lam
            phim = KAPPA * z * dU / ust
            phih = KAPPA * z * dT * ust / (-wth)
        sel = (z > 2) & (z < 80) & np.isfinite(zeta) & (zeta > 0) & (wth < 0)
        out.append(dict(hours_cooling=0.5 * (t1 + t2 - 2 * T_SPIN) / 3600, z=z[sel].tolist(), zeta=zeta[sel].tolist(),
                        phim=phim[sel].tolist(), phih=phih[sel].tolist(), ustar_surface=float(np.sqrt(max(tau[0], 0))),
                        wth_surface=float(wth[0])))
    return out
