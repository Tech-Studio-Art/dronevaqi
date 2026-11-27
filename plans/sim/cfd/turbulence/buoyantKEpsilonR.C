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
    This file is part of OpenFOAM (derived work), GNU GPL v3 or later.
\*---------------------------------------------------------------------------*/

#include "buoyantKEpsilonR.H"
#include "gravityMeshObject.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
buoyantKEpsilonR<BasicTurbulenceModel>::buoyantKEpsilonR
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName, type
    ),
    Cg_(dimensioned<scalar>::lookupOrAddToDict("Cg", this->coeffDict_, 1.0)),
    C3u_(dimensioned<scalar>::lookupOrAddToDict("C3u", this->coeffDict_, 1.0))
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


template<class BasicTurbulenceModel>
bool buoyantKEpsilonR<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        Cg_.readIfPresent(this->coeffDict());
        C3u_.readIfPresent(this->coeffDict());
        return true;
    }
    return false;
}


// Same as buoyantKEpsilon: Gcoef > 0 in stable, < 0 in unstable stratification
template<class BasicTurbulenceModel>
tmp<volScalarField> buoyantKEpsilonR<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(this->mesh_.time());

    return
        (Cg_*this->Cmu_)*this->alpha_*this->k_*(g & fvc::grad(this->rho_))
       /(this->epsilon_ + this->epsilonMin_);
}


// Same as buoyantKEpsilon
template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> buoyantKEpsilonR<BasicTurbulenceModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(this->mesh_.time());

    if (mag(g.value()) > SMALL)
    {
        return -fvm::SuSp(Gcoef(), this->k_);
    }
    return kEpsilon<BasicTurbulenceModel>::kSource();
}


// Changed: C3 = tanh(|v|/|u|) in stable stratification, C3 = C3u in unstable
template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> buoyantKEpsilonR<BasicTurbulenceModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(this->mesh_.time());

    if (mag(g.value()) > SMALL)
    {
        vector gHat(g.value()/mag(g.value()));

        volScalarField v(gHat & this->U_);
        volScalarField u
        (
            mag(this->U_ - gHat*v)
          + dimensionedScalar("SMALL", dimVelocity, SMALL)
        );

        const volScalarField G(Gcoef());
        const volScalarField C3(pos0(G)*tanh(mag(v)/u) + neg(G)*C3u_.value());

        return -fvm::SuSp(this->C1_*C3*G, this->epsilon_);
    }
    return kEpsilon<BasicTurbulenceModel>::epsilonSource();
}

} // End namespace RASModels
} // End namespace Foam
