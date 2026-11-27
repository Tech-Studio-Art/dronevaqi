/*---------------------------------------------------------------------------*\
    Registers buoyantKEpsilonR as a compressible (fluidThermo) RAS model.
    Load with  libs ("libbuoyantKEpsilonR.so");  in system/controlDict.
    License: GNU GPL v3 or later (derived from OpenFOAM v1912).
\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels.H"
#include "buoyantKEpsilonR.H"

makeRASModel(buoyantKEpsilonR);
