#include "MassBalanceTerms.H"
#include "volFields.H"

Foam::scalarField Foam::MassBalanceTerms::residualValues
(
    const volScalarField& massBalance,
    const bool relative
)
{
    scalarField values(mag(massBalance.primitiveField()));

    if (relative)
    {
        const dimensionedScalar dimensionedSmall
        (
            "",
            massBalance.dimensions(),
            SMALL
        );

        values =
            mag(values)
          / max(massBalance.oldTime(), dimensionedSmall)().primitiveField();
    }

    return values;
}
