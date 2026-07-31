/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "poroCouplingTerms.H"
#include "ddtScheme.H"
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

Foam::tmp<Foam::fvMatrix<Foam::scalar>>
Foam::poroCouplingTerms::implicitCouplingMatrix
(
    const volScalarField& implicitCoeff,
    const volScalarField& pField,
    Istream& ddtSchemeData
)
{
    return fv::ddtScheme<scalar>::New
    (
        pField.mesh(),
        ddtSchemeData
    ).ref().fvmDdt(implicitCoeff, pField);
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
    const volScalarField& pRef,
    const fvMatrix<scalar>& implicitDdt,
    const volScalarField& explicitCoupling
)
{
    if (implicitDdt.dimensions() != eqn.dimensions())
    {
        FatalErrorInFunction
            << "Fixed-stress ddt matrix dimensions "
            << implicitDdt.dimensions()
            << " do not match pressure equation dimensions "
            << eqn.dimensions()
            << exit(FatalError);
    }

    // A ddt matrix is diagonal, so use its diagonal directly instead of
    // allocating an A() field and temporary Sp/Su matrices. Its source is
    // deliberately ignored: all physical-time history cancels from the
    // fixed-stress difference, leaving A_t*(p - pRef).
    const scalarField& implicitDiag = implicitDdt.diag();

    // This matrix is returned through fvOptions() on the RHS of the
    // pressure equation, so it is subtracted when `lhs == fvOptions(...)`
    // is assembled.  Build the opposite sign here to add
    // A_t*(p - pRef) as positive storage in the final pEqn, where A_t is
    // the current-time coefficient of ddt(L,p) and pRef is the pressure at
    // the last coupling-term assembly (see header).
    eqn.diag() -= implicitDiag;
    eqn.source() -= implicitDiag*pRef.primitiveField();

    const tmp<volScalarField> tExplicitRate
    (
        explicitCouplingRate(explicitCoupling)
    );

    eqn += fvm::Su(tExplicitRate(), pField);
}

// ************************************************************************* //
