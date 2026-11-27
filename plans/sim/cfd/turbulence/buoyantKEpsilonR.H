/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
    Modified 2026 for the Drone V AQI study (Tech Studio .Art): buoyantKEpsilonR
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM (derived work).

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::buoyantKEpsilonR

Description
    buoyantKEpsilon (OpenFOAM v1912) with one change in the epsilon equation.

    In buoyantKEpsilon the buoyancy term of the epsilon equation is
        C1*C3*G_b*epsilon/k  with  C3 = tanh(|v|/|u|)
    (v: velocity along gravity, u: velocity normal to gravity). In horizontal
    flow C3 -> 0. In unstable stratification G_b > 0 then feeds k but not
    epsilon, so k/epsilon (and nut) can grow without limit.

    buoyantKEpsilonR keeps C3 = tanh(|v|/|u|) in stable stratification (the
    calibrated case) and uses C3 = C3u (default 1) in unstable stratification,
    as recommended by Rodi (1987), J. Geophys. Res. 92(C5), 5305-5328.

        buoyantKEpsilonRCoeffs { Cg 1.0; C3u 1.0; }

SourceFiles
    buoyantKEpsilonR.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantKEpsilonR_H
#define buoyantKEpsilonR_H

#include "kEpsilon.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class buoyantKEpsilonR
:
    public kEpsilon<BasicTurbulenceModel>
{
    buoyantKEpsilonR(const buoyantKEpsilonR&) = delete;
    void operator=(const buoyantKEpsilonR&) = delete;

protected:

        dimensionedScalar Cg_;
        dimensionedScalar C3u_;

        tmp<volScalarField> Gcoef() const;

        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("buoyantKEpsilonR");

        buoyantKEpsilonR
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );

    virtual ~buoyantKEpsilonR() = default;

        virtual bool read();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "buoyantKEpsilonR.C"
#endif

#endif
