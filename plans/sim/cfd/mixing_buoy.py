#!/usr/bin/env python3
"""
Step 1c: inversion mixing with a buoyancy-corrected turbulence model.

Same unit cell, mesh, drone actuator disk, inversion and PM profile as mixing_axi.py.
Changes:
  * solver buoyantPimpleFoam with turbulence model buoyantKEpsilon (buoyancy production /
    destruction of turbulence, i.e. stable air damps the eddies);
  * Boussinesq equation of state and a large specific heat (1e5 J/kg/K): T is the potential
    temperature (adiabatic lapse rate 0.1 K/km, negligible) - the potential-temperature form
    used in atmospheric boundary-layer CFD, checked in the stable-surface-layer test
    (mosttest.py);
  * buoyancy coefficient Cg set from that test (see runs/validation/most_*).
PM is still read from theta (same transport equation, Pr_t = Sc_t).
"""
import json
import math
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mixing_axi as mx  # noqa: E402
from foam import write, run, field, write_field  # noqa: E402

RUNS = os.path.join(HERE, "runs", "mixing_buoy")
CP_BIG = 1.0e5
RHO0 = 1.25
CG = float(os.environ.get("MIX_CG", "1.0"))
PRT = 0.7

CASES = {k: dict(v) for k, v in mx.CASES.items() if k in ("heavy_z60", "heavy_z90", "quad_z30", "no_drone")}
CASES["heavy_z60_coarse"] = dict(mx.CASES["heavy_z60_gc"])


def convert(case, c):
    """Turn a mixing_axi case into a buoyantPimpleFoam / buoyantKEpsilon case."""
    os.remove(os.path.join(case, "constant/transportProperties"))
    write(case, "constant/thermophysicalProperties", f"""
thermoType {{ type heRhoThermo; mixture pureMixture; transport const; thermo hConst;
  equationOfState Boussinesq; specie specie; energy sensibleEnthalpy; }}
dpdt no;
mixture {{
  specie {{ molWeight 28.96; }}
  equationOfState {{ rho0 {RHO0}; T0 {mx.THETA0}; beta {1 / mx.THETA0:.8g}; }}
  thermodynamics {{ Cp {CP_BIG}; Hf 0; }}
  transport {{ mu {mx.NU * RHO0:.6g}; Pr 0.7; }}
}}
""")
    write(case, "constant/turbulenceProperties",
          f"simulationType RAS;\nRAS {{ RASModel buoyantKEpsilon; turbulence on; printCoeffs on; "
          f"buoyantKEpsilonCoeffs {{ Cg {CG}; Prt {PRT}; }} }}\n")
    # actuator disk: compressible momentum equation -> force in N (not kinematic)
    fv = os.path.join(case, "constant/fvOptions")
    s = open(fv).read()
    if "injectionRateSuSp" in s:
        s = re.sub(r"U \(\(0 0 (-?[0-9.e+-]+)\) 0\)",
                   lambda m: f"U ((0 0 {float(m.group(1)) * mx.RHO:.6g}) 0)", s)
        open(fv, "w").write(s)
    cd = os.path.join(case, "system/controlDict")
    s = open(cd).read().replace("application buoyantBoussinesqPimpleFoam;", "application buoyantPimpleFoam;")
    open(cd, "w").write(s)
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
    sym = {"outer": "        type symmetry;", "axis": "        type symmetry;",
           "front": "        type wedge;", "back": "        type wedge;"}

    def pf(walltxt):
        return dict(sym, ground=walltxt, top=walltxt)
    write_field(case, "0", field("p_rgh", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5",
                                 pf("        type fixedFluxPressure; value uniform 1e5;")))
    write_field(case, "0", field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5",
                                 pf("        type calculated; value uniform 1e5;")))
    write_field(case, "0", field("alphat", "volScalarField", "[1 -1 -1 0 0 0 0]", "0",
                                 dict(sym, ground=f"        type compressible::alphatJayatillekeWallFunction; Prt {PRT}; value uniform 0;",
                                      top="        type calculated; value uniform 0;")))


def solve(name):
    c = CASES[name]
    case = os.path.join(RUNS, name + os.environ.get("MIX_TAG", ""))
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    mx.build(case, c)
    convert(case, c)
    run(case, "blockMesh")
    run(case, "topoSet")
    n = mx.set_initial_theta(case)
    run(case, "buoyantPimpleFoam", "solve")
    res = dict(case=name + os.environ.get("MIX_TAG", ""), Cg=CG, cells=n, thrust_N=c["T"], hover_z=c["zh"], series=mx.analyse(case))
    json.dump(res, open(os.path.join(case, "result.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    names = sys.argv[1:] or ["heavy_z60", "quad_z30", "no_drone", "heavy_z90"]
    os.makedirs(RUNS, exist_ok=True)
    with ThreadPoolExecutor(max_workers=int(os.environ.get("NPROC", "2"))) as ex:
        for r in ex.map(solve, names):
            s = r["series"]
            print(r["case"], "Cg", r["Cg"], [round(100 * x["z<2m"] / s[0]["z<2m"], 1) for x in s])
