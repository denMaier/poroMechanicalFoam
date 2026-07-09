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

    //- Add the fixed-stress and explicit deformation coupling contributions
    //  to the pressure equation. pRef is the pressure reference the
    //  fixed-stress term relaxes against: it must be the pressure at which
    //  the coupling terms (implicitCoupling/explicitCoupling) were
    //  assembled, i.e. a snapshot taken once per coupling iteration. Do NOT
    //  pass pField.prevIter() here: that slot is re-stored every fluid
    //  iteration for relaxation, which would make the fixed-stress term
    //  vanish as the fluid sub-loop converges.
    void addCouplingSource
    (
        fvMatrix<scalar>& eqn,
        const volScalarField& pField,
        const volScalarField& pRef,
        const volScalarField& implicitCoupling,
        const volScalarField& explicitCoupling,
        const dimensionedScalar& deltaT
    );
}
}

#endif

// ************************************************************************* //
