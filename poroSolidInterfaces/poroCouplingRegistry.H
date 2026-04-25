/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef poroCouplingRegistry_H
#define poroCouplingRegistry_H

#include "HashSet.H"
#include "objectRegistry.H"
#include "volFields.H"

namespace Foam
{
namespace poroCouplingRegistry
{
    bool ensureVolScalarFieldRegistered
    (
        objectRegistry& targetRegistry,
        volScalarField& field,
        const word& ownerRegistryName,
        HashSet<word>& borrowedFieldNames
    );
}
}

#endif

// ************************************************************************* //
