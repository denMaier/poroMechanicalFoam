/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "poroCouplingTerms.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvmSup.H"

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

void Foam::poroCouplingTerms::updateCouplingFields
(
    const volScalarField& explicitCouplingCoeff,
    const volScalarField& stabilCouplingCoeff,
    const volScalarField& impK,
    const volVectorField& U,
    autoPtr<volScalarField>& nDotField,
    autoPtr<volScalarField>& fixedStressStabilField
)
{
    if(!nDotField.valid())
    {
        tmp<volScalarField> tnDot(nDot(explicitCouplingCoeff, U));
        nDotField.reset(tnDot.ptr());

        if
        (
            !nDotField().mesh().objectRegistry::foundObject<volScalarField>
            (
                nDotField().name()
            )
        )
        {
            nDotField().mesh().objectRegistry::checkIn(nDotField());
        }
    }
    else
    {
        nDotField.ref() = nDot(explicitCouplingCoeff, U);
    }

    const tmp<volScalarField> tStabil
    (
        fixedStressStabil(stabilCouplingCoeff, impK)
    );
    fixedStressStabilField.reset(new volScalarField(tStabil()));
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

void Foam::poroCouplingTerms::addCouplingSource
(
    fvMatrix<scalar>& eqn,
    const volScalarField& pField,
    const volScalarField& implicitCoupling,
    const volScalarField& explicitCoupling,
    const dimensionedScalar& deltaT
)
{
    const tmp<volScalarField> tImplicitRate
    (
        implicitCouplingRate(implicitCoupling, deltaT)
    );

    eqn += fvm::SuSp(tImplicitRate(), pField);

    const tmp<volScalarField> tExplicitRate
    (
        explicitCouplingRate(explicitCoupling)
    );

    eqn += fvm::Su(tExplicitRate(), pField);
}

// ************************************************************************* //
