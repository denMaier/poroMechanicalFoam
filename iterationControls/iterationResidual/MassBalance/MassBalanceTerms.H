#ifndef MassBalanceTerms_H
#define MassBalanceTerms_H

#include "scalarField.H"
#include "volFieldsFwd.H"

namespace Foam
{
namespace MassBalanceTerms
{
    scalarField residualValues
    (
        const volScalarField& massBalance,
        const bool relative
    );
}
}

#endif
