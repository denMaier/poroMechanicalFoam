/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "poroSolidStab.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "poroPressureUnits.H"

namespace Foam
{
namespace poroSolidInteractions
{
    defineTypeNameAndDebug(poroSolidStab, 0);
    addToRunTimeSelectionTable(poroSolidInterface, poroSolidStab, dictionary);

    poroSolidStab::poroSolidStab
    (
        Time& runTime,
        const word& region
    )
    :
        poroSolid(typeName, runTime, region),
        pressureStabilization_
        (
            interactionProperties().lookupOrDefault<Switch>
            (
                "pressureStabilization",
                true
            )
        ),
        pressureStabilizationTargetK_
        (
            interactionProperties().lookupOrDefault<dimensionedScalar>
            (
                "pressureStabilizationTargetK",
                dimensionedScalar
                (
                    "pressureStabilizationTargetK",
                    dimLength/dimTime,
                    0.0
                )
            )
        )
    {
        if (pressureStabilization_)
        {
            Info<< "Pressure stabilization enabled in poroSolidStab" << nl
                << "Pressure stabilization target k: "
                << pressureStabilizationTargetK_ << endl;
        }
    }

    void poroSolidStab::addPressureEquationTerms
    (
        fvMatrix<scalar>& eqn
    ) const
    {
        poroSolid::addPressureEquationTerms(eqn);

        if (!pressureStabilization_)
        {
            return;
        }

        const tmp<surfaceScalarField> tkf
        (
            poroFluid().relativeAccelerationConductivity()
        );

        const dimensionedScalar pressureScale
        (
            poroPressureUnits::pressureScale
            (
                pField().dimensions(),
                poroFluid()
            )
        );

        const tmp<surfaceScalarField> tAddedByScale
        (
            max
            (
                pressureStabilizationTargetK_/pressureScale - tkf()/pressureScale,
                dimensionedScalar
                (
                    "zero",
                    tkf().dimensions()/pressureScale.dimensions(),
                    0.0
                )
            )
        );

        // Deferred correction: the implicit/explicit pair must cancel when
        // the fluid loop converges, so prevIter() is the right reference
        // here — the fluid model stores it at the top of every fluid
        // iteration, before this term is assembled through fvOptions. Do
        // not replace it with pCouplingRef(): that snapshot is frozen per
        // coupling iteration and would leave a spurious residual flux.
        eqn -= fvm::laplacian(tAddedByScale(), pField());
        eqn += fvc::laplacian(tAddedByScale(), pField().prevIter());
    }
}
}

// ************************************************************************* //
