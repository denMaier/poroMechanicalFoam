/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef poroPressureUnits_H
#define poroPressureUnits_H

#include "dimensionedScalar.H"
#include "objectRegistry.H"

namespace Foam
{
namespace poroPressureUnits
{
    dimensionedScalar pressureScale
    (
        const dimensionSet& pressureFieldDimensions,
        const objectRegistry& registry
    );
}
}

#endif

// ************************************************************************* //
