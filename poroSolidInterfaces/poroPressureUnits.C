/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "poroPressureUnits.H"
#include "UniformDimensionedField.H"
#include "volFields.H"

Foam::dimensionedScalar Foam::poroPressureUnits::pressureScale
(
    const dimensionSet& pressureFieldDimensions,
    const objectRegistry& registry
)
{
    if(pressureFieldDimensions == dimPressure)
    {
        return dimensionedScalar("gammaw", dimless, 1.0);
    }

    if(pressureFieldDimensions == dimLength)
    {
        if(!registry.foundObject<UniformDimensionedField<vector>>("gamma_water"))
        {
            FatalErrorInFunction
                << "Head-based variably saturated coupling requires a "
                << "'gamma_water' uniform field once the conversion factor "
                << "is actually needed"
                << exit(FatalError);
        }

        const UniformDimensionedField<vector>& gammaW =
            registry.lookupObject<UniformDimensionedField<vector>>("gamma_water");

        return dimensionedScalar
        (
            "gammaw",
            gammaW.dimensions(),
            mag(gammaW.value())
        );
    }

    FatalErrorInFunction
        << "pore pressure field is neither in head nor in pressure dimensions"
        << exit(FatalError);

    return dimensionedScalar("gammaw", dimless, 1.0);
}

// ************************************************************************* //
