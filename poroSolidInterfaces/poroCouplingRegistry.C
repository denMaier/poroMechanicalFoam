/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "poroCouplingRegistry.H"

bool Foam::poroCouplingRegistry::ensureVolScalarFieldRegistered
(
    objectRegistry& targetRegistry,
    volScalarField& field,
    const word& ownerRegistryName,
    HashSet<word>& borrowedFieldNames
)
{
    if(targetRegistry.foundObject<volScalarField>(field.name()))
    {
        const volScalarField& registeredField =
            targetRegistry.lookupObject<volScalarField>(field.name());

        if(&registeredField != &field)
        {
            FatalErrorInFunction
                << "Target registry already contains a different volScalarField named '"
                << field.name() << "' while trying to expose the field from '"
                << ownerRegistryName << "'" << nl
                << "Existing object db: " << registeredField.db().name() << nl
                << "Requested object db: " << field.db().name()
                << exit(FatalError);
        }

        return false;
    }

    targetRegistry.checkIn(field);
    borrowedFieldNames.insert(field.name());
    return true;
}

// ************************************************************************* //
