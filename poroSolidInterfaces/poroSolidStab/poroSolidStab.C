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
#include "mechanicalModel.H"
#include "syncTools.H"

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
        stabilizationType_
        (
            interactionProperties().lookupOrDefault<word>
            (
                "stabilizationType",
                "implicit"
            )
        ),
        stabFactor_
        (
            interactionProperties().lookupOrDefault<scalar>
            (
                "stabFactor",
                1.0
            )
        ),
        excludedCellZones_
        (
            interactionProperties().lookupOrDefault<wordList>
            (
                "excludedCellZones",
                wordList()
            )
        ),
        exclusionSummaryPrinted_(false),
        checkerboardStabilCoeff_()
    {
        if
        (
            stabilizationType_ != "implicit"
         && stabilizationType_ != "explicit"
         && stabilizationType_ != "pressureJump"
        )
        {
            FatalIOErrorInFunction(interactionProperties())
                << "Unknown stabilizationType '" << stabilizationType_ << "'. "
                << "Valid choices are 'implicit', 'explicit', and "
                << "'pressureJump'."
                << exit(FatalIOError);
        }

        if (stabFactor_ < 0.0)
        {
            FatalIOErrorInFunction(interactionProperties())
                << "stabFactor must be non-negative, but is " << stabFactor_
                << exit(FatalIOError);
        }

        Info<< "Checkerboard pressure stabilization enabled in poroSolidStab" << nl
            << "    stabilizationType: " << stabilizationType_ << nl
            << "    stabFactor: " << stabFactor_ << nl
            << "    excludedCellZones: " << excludedCellZones_ << endl;
    }

    boolList poroSolidStab::excludedStabilizationFaces
    (
        const fvMesh& mesh,
        label& excludedCellCount
    ) const
    {
        boolList excludedCells(mesh.nCells(), false);

        forAll(excludedCellZones_, zoneNameI)
        {
            const word& zoneName = excludedCellZones_[zoneNameI];
            const label zoneI = mesh.cellZones().findZoneID(zoneName);

            if (zoneI < 0)
            {
                FatalErrorInFunction
                    << "Cannot exclude stabilization around cell zone '"
                    << zoneName << "' because it does not exist on mesh '"
                    << mesh.name() << "'. Available cell zones are "
                    << mesh.cellZones().names()
                    << exit(FatalError);
            }

            const cellZone& zone = mesh.cellZones()[zoneI];
            forAll(zone, zoneCellI)
            {
                excludedCells[zone[zoneCellI]] = true;
            }
        }

        excludedCellCount = 0;
        forAll(excludedCells, cellI)
        {
            excludedCellCount += excludedCells[cellI] ? 1 : 0;
        }

        boolList excludedFaces(mesh.nFaces(), false);
        const labelUList& owner = mesh.faceOwner();
        const labelUList& neighbour = mesh.faceNeighbour();

        for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
        {
            excludedFaces[faceI] =
                excludedCells[owner[faceI]]
             || excludedCells[neighbour[faceI]];
        }

        // Coupled faces need the exclusion state from the cell on the other
        // side.  swapBoundaryCellList handles processor and cyclic patches
        // without turning the exact Boolean face rule into an interpolation.
        const boolList neighbourExcluded
        (
            syncTools::swapBoundaryCellList(mesh, excludedCells)
        );

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];
            if (!patch.coupled())
            {
                continue;
            }

            forAll(patch, patchFaceI)
            {
                const label faceI = patch.start() + patchFaceI;
                const label boundaryFaceI = faceI - mesh.nInternalFaces();
                excludedFaces[faceI] =
                    excludedCells[owner[faceI]]
                 || neighbourExcluded[boundaryFaceI];
            }
        }

        return excludedFaces;
    }

    void poroSolidStab::updateCheckerboardStabilCoeff
    (
        const volScalarField& rImpK
    )
    {
        const fvMesh& mesh = poroFluidMesh();
        tmp<surfaceScalarField> tCoeff
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "checkerboardStabilCoeff",
                    runTime().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
               stabFactor_
               *sqr(1.0/mesh.deltaCoeffs())
               *linearInterpolate(rImpK)
            )
        );

        // Do not permit the numerical stabilization to create an external
        // pressure flux. Coupled boundaries retain the internal face
        // treatment; all physical boundaries receive zero coefficient.
        forAll(tCoeff().boundaryField(), patchI)
        {
            if (!tCoeff().boundaryField()[patchI].coupled())
            {
                tCoeff.ref().boundaryFieldRef()[patchI] = 0.0;
            }
        }

        if (!excludedCellZones_.empty())
        {
            label localExcludedCellCount = 0;
            const boolList excludedFaces
            (
                excludedStabilizationFaces(mesh, localExcludedCellCount)
            );
            scalarField& internalCoeff = tCoeff.ref().primitiveFieldRef();

            forAll(internalCoeff, faceI)
            {
                if (excludedFaces[faceI])
                {
                    internalCoeff[faceI] = 0.0;
                }
            }

            scalar localExcludedFaceCount = 0.0;
            for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
            {
                localExcludedFaceCount += excludedFaces[faceI] ? 1.0 : 0.0;
            }

            forAll(tCoeff().boundaryField(), patchI)
            {
                scalarField& patchCoeff =
                    tCoeff.ref().boundaryFieldRef()[patchI];
                const polyPatch& patch = mesh.boundaryMesh()[patchI];

                forAll(patchCoeff, patchFaceI)
                {
                    const label faceI = patch.start() + patchFaceI;
                    if (excludedFaces[faceI])
                    {
                        patchCoeff[patchFaceI] = 0.0;
                        // A coupled connection exists on both sides, so each
                        // local representation contributes half a global face.
                        localExcludedFaceCount += 0.5;
                    }
                }
            }

            if (!exclusionSummaryPrinted_)
            {
                const label globalExcludedCellCount = returnReduce
                (
                    localExcludedCellCount,
                    sumOp<label>()
                );
                const scalar globalExcludedFaceCount = returnReduce
                (
                    localExcludedFaceCount,
                    sumOp<scalar>()
                );

                Info<< "Stabilization face exclusion: zones "
                    << excludedCellZones_ << ", cells "
                    << globalExcludedCellCount << ", internal/coupled faces "
                    << label(globalExcludedFaceCount + 0.5) << endl;
                exclusionSummaryPrinted_ = true;
            }
        }

        if (!checkerboardStabilCoeff_.valid())
        {
            checkerboardStabilCoeff_.reset(tCoeff.ptr());
        }
        else
        {
            checkerboardStabilCoeff_.ref() = tCoeff();
        }
    }

    void poroSolidStab::assembleCouplingTerms()
    {
        poroSolid::assembleCouplingTerms();

        const volScalarField impKSolid(solid().mechanical().impK());
        const dimensionedScalar minImpK
        (
            "minImpK",
            impKSolid.dimensions(),
            VSMALL
        );
        const volScalarField rImpKSolid
        (
            IOobject
            (
                "checkerboardRImpK",
                runTime().timeName(),
                solidMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            1.0/max(impKSolid, minImpK)
        );

        if (sharedMesh())
        {
            updateCheckerboardStabilCoeff(rImpKSolid);
        }
        else
        {
            // Map compliance, not stiffness. Linear interpolation of 1/impK
            // is equivalent to a weighted harmonic face interpolation of
            // impK and avoids taking a reciprocal after mesh mapping.
            const tmp<volScalarField> trImpKFluid
            (
                solidToPoroFluid().mapTgtToSrc(rImpKSolid)
            );
            updateCheckerboardStabilCoeff(trImpKFluid());
        }
    }

    void poroSolidStab::clearCouplingTerms()
    {
        poroSolid::clearCouplingTerms();
        checkerboardStabilCoeff_.clear();
    }

    void poroSolidStab::addPressureEquationTerms
    (
        fvMatrix<scalar>& eqn
    ) const
    {
        poroSolid::addPressureEquationTerms(eqn);

        if (!checkerboardStabilCoeff_.valid())
        {
            FatalErrorInFunction
                << "Checkerboard stabilization coefficient is not initialized"
                << exit(FatalError);
        }

        const surfaceScalarField& coeff = checkerboardStabilCoeff_();
        const fvMesh& mesh = pField().mesh();

        // The coupled mass balance is based on the finite displacement
        // increment over the current time step, not on a selectable ddt(D)
        // scheme.  Use the matching Euler pressure increments here.
        const volScalarField pDotPredictor
        (
            IOobject
            (
                "checkerboardPredictorPressureRate",
                runTime().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (pCouplingRef() - pField().oldTime())/runTime().deltaT()
        );
        const tmp<volVectorField> tGradPredictorRate
        (
            fvc::grad(pDotPredictor)
        );
        const tmp<volScalarField> tWideRate
        (
            fvc::div
            (
                coeff
               *(
                    mesh.Sf()
                  & fvc::interpolate(tGradPredictorRate())
                )
            )
        );

        if (stabilizationType_ == "pressureJump")
        {
            const tmp<surfaceScalarField> tCoeffByDt
            (
                coeff/runTime().deltaT()
            );

            // Paper-aligned pressure-jump stabilization.  Keep the complete
            // pressure-increment diffusion operator inside the flow solve;
            // unlike the compact-minus-wide variants, this introduces no
            // outer-iteration-lagged wide-stencil term.
            eqn += fvm::laplacian(tCoeffByDt(), pField());
            eqn -= fvc::laplacian(tCoeffByDt(), pField().oldTime());
        }
        else if (stabilizationType_ == "implicit")
        {
            const tmp<surfaceScalarField> tCoeffByDt
            (
                coeff/runTime().deltaT()
            );

            // Implicit compact-minus-wide stabilization of the temporal
            // pressure increment.  The coupling snapshot is frozen while
            // the fluid iterates, so the wide predictor is time-consistent.
            // addPressureEquationTerms() builds the fvOptions matrix on the
            // right-hand side of `lhs == fvOptions(...)`.  Put
            // compact-minus-wide into that option; OpenFOAM subtracts it
            // during equation assembly, yielding the stabilizing
            // wide-minus-compact contribution on the final LHS.
            eqn += fvm::laplacian(tCoeffByDt(), pField());
            eqn -= fvc::laplacian(tCoeffByDt(), pField().oldTime());
            eqn -= tWideRate();
        }
        else
        {
            // Fully lagged compact-minus-wide replacement.  Both stencils
            // act on the temporal pressure increment that produced the
            // incoming solid predictor.  This fvOption contains
            // compact-minus-wide and is subtracted from the final pressure
            // equation, just like the implicit form above.
            eqn += fvc::laplacian(coeff, pDotPredictor);
            eqn -= tWideRate();
        }
    }
}
}

// ************************************************************************* //
