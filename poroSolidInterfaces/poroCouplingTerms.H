/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef poroCouplingTerms_H
#define poroCouplingTerms_H

#include "volFields.H"
#include "surfaceFields.H"
#include "tmp.H"
#include "autoPtr.H"

namespace Foam
{
template<class Type>
class fvMatrix;

namespace poroCouplingTerms
{
    tmp<volScalarField> nDot
    (
        const volScalarField& couplingCoeff,
        const volVectorField& U
    );

    tmp<volScalarField> fixedStressStabil
    (
        const volScalarField& biotCoeff,
        const volScalarField& impK
    );

    void updateCouplingFields
    (
        const volScalarField& explicitCouplingCoeff,
        const volScalarField& stabilCouplingCoeff,
        const volScalarField& impK,
        const volVectorField& U,
        autoPtr<volScalarField>& nDotField,
        autoPtr<volScalarField>& fixedStressStabilField
    );

    tmp<surfaceVectorField> relativeAccelerationFlux
    (
        const surfaceScalarField& kf,
        const volVectorField& a,
        const dimensionedScalar& magG
    );

    tmp<volScalarField> explicitCouplingSource
    (
        const volScalarField& nDot,
        const surfaceVectorField& qRelAcc
    );

    tmp<volScalarField> implicitCouplingRate
    (
        const volScalarField& implicitCoeff,
        const dimensionedScalar& deltaT
    );

    tmp<volScalarField> explicitCouplingRate
    (
        const volScalarField& explicitSource
    );

    void addCouplingSource
    (
        fvMatrix<scalar>& eqn,
        const volScalarField& pField,
        const volScalarField& implicitCoupling,
        const volScalarField& explicitCoupling,
        const dimensionedScalar& deltaT
    );
}
}

#endif

// ************************************************************************* //
