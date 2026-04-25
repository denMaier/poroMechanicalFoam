/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef varSatPoroMechanicalLawTerms_H
#define varSatPoroMechanicalLawTerms_H

#include "volFields.H"
#include "surfaceFields.H"
#include "tmp.H"

namespace Foam
{
namespace varSatPoroMechanicalLawTerms
{
    tmp<volScalarField> mixtureDensity
    (
        const dimensionedScalar& saturatedDensity,
        const volScalarField& S,
        const volScalarField& n,
        const dimensionedScalar& rhoWater
    );

    tmp<volSymmTensorField> initialEffectiveStress
    (
        const volSymmTensorField& totalStress,
        const dimensionedScalar& biotCoeff,
        const volScalarField& chi,
        const volScalarField& p
    );

    tmp<surfaceSymmTensorField> initialEffectiveStress
    (
        const surfaceSymmTensorField& totalStress,
        const dimensionedScalar& biotCoeff,
        const surfaceScalarField& chi,
        const surfaceScalarField& p
    );

    tmp<volSymmTensorField> totalStress
    (
        const volSymmTensorField& effectiveStress,
        const dimensionedScalar& biotCoeff,
        const volScalarField& chi,
        const volScalarField& p
    );

    tmp<surfaceSymmTensorField> totalStress
    (
        const surfaceSymmTensorField& effectiveStress,
        const dimensionedScalar& biotCoeff,
        const surfaceScalarField& chi,
        const surfaceScalarField& p
    );
}
}

#endif

// ************************************************************************* //
