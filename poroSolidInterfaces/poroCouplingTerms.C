/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "poroCouplingTerms.H"
#include "fvc.H"

Foam::tmp<Foam::volScalarField> Foam::poroCouplingTerms::nDot
(
    const volScalarField& couplingCoeff,
    const volVectorField& U
)
{
    const tmp<volScalarField> tDivU(fvc::div(U));

    tmp<volScalarField> tnDot
    (
        new volScalarField
        (
            IOobject
            (
                "nDot",
                U.mesh().time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            couplingCoeff*tDivU()
        )
    );

    return tnDot;
}

Foam::tmp<Foam::volScalarField> Foam::poroCouplingTerms::fixedStressStabil
(
    const volScalarField& biotCoeff,
    const volScalarField& impK
)
{
    tmp<volScalarField> tStabil
    (
        new volScalarField
        (
            "fixedStressSplitCoeff",
            min
            (
                pow(biotCoeff, 2)
               /max
                (
                    impK,
                    dimensionedScalar("small", impK.dimensions(), VSMALL)
                ),
                dimensionedScalar("great", dimless/impK.dimensions(), VGREAT)
            )
        )
    );

    return tStabil;
}

Foam::tmp<Foam::surfaceVectorField>
Foam::poroCouplingTerms::relativeAccelerationFlux
(
    const surfaceScalarField& kf,
    const volVectorField& a,
    const dimensionedScalar& magG
)
{
    tmp<surfaceVectorField> tq
    (
        new surfaceVectorField
        (
            "q_relAcc",
            kf*fvc::interpolate(a/magG)
        )
    );

    return tq;
}

Foam::tmp<Foam::volScalarField> Foam::poroCouplingTerms::explicitCouplingSource
(
    const volScalarField& nDot,
    const surfaceVectorField& qRelAcc
)
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            nDot - fvc::div(nDot.mesh().Sf() & qRelAcc)
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::poroCouplingTerms::implicitCouplingRate
(
    const volScalarField& implicitCoeff,
    const dimensionedScalar& deltaT
)
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "implicitCouplingRate",
            implicitCoeff/deltaT
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::poroCouplingTerms::explicitCouplingRate
(
    const volScalarField& explicitSource
)
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "explicitCouplingRate",
            -explicitSource
        )
    );
}

// ************************************************************************* //
