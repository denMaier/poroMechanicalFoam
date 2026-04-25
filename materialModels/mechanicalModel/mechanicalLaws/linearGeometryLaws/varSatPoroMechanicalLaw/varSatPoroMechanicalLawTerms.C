/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "varSatPoroMechanicalLawTerms.H"

Foam::tmp<Foam::volScalarField>
Foam::varSatPoroMechanicalLawTerms::mixtureDensity
(
    const dimensionedScalar& saturatedDensity,
    const volScalarField& S,
    const volScalarField& n,
    const dimensionedScalar& rhoWater
)
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                S.mesh().time().timeName(),
                S.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            saturatedDensity - (scalar(1) - S)*n*rhoWater,
            "zeroGradient"
        )
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::varSatPoroMechanicalLawTerms::initialEffectiveStress
(
    const volSymmTensorField& totalStress,
    const dimensionedScalar& biotCoeff,
    const volScalarField& chi,
    const volScalarField& p
)
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            totalStress + biotCoeff*chi*p*symmTensor(I)
        )
    );
}

Foam::tmp<Foam::surfaceSymmTensorField>
Foam::varSatPoroMechanicalLawTerms::initialEffectiveStress
(
    const surfaceSymmTensorField& totalStress,
    const dimensionedScalar& biotCoeff,
    const surfaceScalarField& chi,
    const surfaceScalarField& p
)
{
    return tmp<surfaceSymmTensorField>
    (
        new surfaceSymmTensorField
        (
            totalStress + biotCoeff*chi*p*symmTensor(I)
        )
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::varSatPoroMechanicalLawTerms::totalStress
(
    const volSymmTensorField& effectiveStress,
    const dimensionedScalar& biotCoeff,
    const volScalarField& chi,
    const volScalarField& p
)
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            effectiveStress - biotCoeff*chi*p*symmTensor(I)
        )
    );
}

Foam::tmp<Foam::surfaceSymmTensorField>
Foam::varSatPoroMechanicalLawTerms::totalStress
(
    const surfaceSymmTensorField& effectiveStress,
    const dimensionedScalar& biotCoeff,
    const surfaceScalarField& chi,
    const surfaceScalarField& p
)
{
    return tmp<surfaceSymmTensorField>
    (
        new surfaceSymmTensorField
        (
            effectiveStress - biotCoeff*chi*p*symmTensor(I)
        )
    );
}

// ************************************************************************* //
