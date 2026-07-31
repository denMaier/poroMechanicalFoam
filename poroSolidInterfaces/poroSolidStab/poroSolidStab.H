/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef poroSolidStab_H
#define poroSolidStab_H

#include "poroSolid.H"

namespace Foam
{
namespace poroSolidInteractions
{

class poroSolidStab
:
    public poroSolid
{
    //- Checkerboard stabilization treatment.
    const word stabilizationType_;

    //- Dimensionless scale multiplying h^2/impK.
    const scalar stabFactor_;

    //- Cell zones whose touching faces receive no artificial stabilization.
    const wordList excludedCellZones_;

    //- Avoid repeating the static zone/face summary every coupling iteration.
    bool exclusionSummaryPrinted_;

    //- Face coefficient stabFactor*h^2/impK on the fluid mesh.
    autoPtr<surfaceScalarField> checkerboardStabilCoeff_;

    //- Mark internal and coupled faces touching any configured cell zone.
    boolList excludedStabilizationFaces
    (
        const fvMesh& mesh,
        label& excludedCellCount
    ) const;

    //- Update the face coefficient from 1/impK on the fluid mesh.
    void updateCheckerboardStabilCoeff(const volScalarField& rImpK);

protected:

    //- Assemble standard coupling terms and the checkerboard coefficient.
    virtual void assembleCouplingTerms();

    //- Clear all per-coupling-iteration fields.
    virtual void clearCouplingTerms();

public:
    TypeName("poroSolidStab");

    poroSolidStab
    (
        Time& runTime,
        const word& region = "poroFluid"
    );

    virtual void addPressureEquationTerms(fvMatrix<scalar>& eqn) const;
};

}
}

#endif

// ************************************************************************* //
