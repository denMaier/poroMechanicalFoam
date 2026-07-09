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
    Switch pressureStabilization_;

    dimensionedScalar pressureStabilizationTargetK_;

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
