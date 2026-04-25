/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "poroSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "mechanicalModel.H"
#include "iterationControl.H"
#include "poroCouplingTerms.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace poroSolidInteractions
    {
        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(poroSolid, 0);
        addToRunTimeSelectionTable(poroSolidInterface, poroSolid, dictionary);

        // * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //
        void poroSolid::initializeSolidHydraulicFields()
        {
            if(sharedMesh())
            {
                ensureSharedSolidFieldRegistered(poroFluidRef().p(), poroFluidMesh().name());
                ensureSharedSolidFieldRegistered(poroFluidRef().p_rgh(), poroFluidMesh().name());
            }
            else
            {
                if(!pSolidMesh_.valid())
                {
                    pSolidMesh_.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "p",
                                runTime().timeName(),
                                solidMesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            solidToPoroFluid().mapSrcToTgt(poroFluid().p())()
                        )
                    );
                }

                if(!pRghSolidMesh_.valid())
                {
                    pRghSolidMesh_.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "p_rgh",
                                runTime().timeName(),
                                solidMesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            solidToPoroFluid().mapSrcToTgt(poroFluid().p_rgh())()
                        )
                    );
                }
            }
        }

        void poroSolid::syncSolidHydraulicFields()
        {
            if(sharedMesh())
            {
                return;
            }

            initializeSolidHydraulicFields();
            mapPressuresToSolidMesh(pSolidMesh_, pRghSolidMesh_);
        }

        void poroSolid::prepareCouplingLoop()
        {
            Info << "Preparing poroSolid solver" << endl;
        }

        void poroSolid::assembleCouplingTerms()
        {
            if(sharedMesh())
            {
                const volScalarField impK(solid().mechanical().impK());
                const tmp<surfaceScalarField> tkf(poroFluid().relativeAccelerationConductivity());
                const tmp<volVectorField> tSolidA(fvc::d2dt2(solid().D()));

                updateCouplingTerms(b(), impK, solid().U(), nDot_, fixedStressStabil_);
                q_relAcc_.reset(q_relAcc(tkf(), tSolidA()).ptr());
            }
            else
            {
                Info << "Mapping fields to poroFluid mesh";
                const volScalarField impK(solid().mechanical().impK());
                const tmp<surfaceScalarField> tkf(poroFluid().relativeAccelerationConductivity());
                tmp<volVectorField> UFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().U());
                const tmp<volVectorField> tSolidA(fvc::d2dt2(solid().D()));
                tmp<volVectorField> aFluidMesh = solidToPoroFluid().mapTgtToSrc(tSolidA());
                tmp<volScalarField> tmpImpK(solidToPoroFluid().mapTgtToSrc(impK));

                updateCouplingTerms(b(), tmpImpK(), UFluidMesh(), nDot_, fixedStressStabil_);
                q_relAcc_.reset(q_relAcc(tkf(), aFluidMesh()).ptr());

                tmpImpK.clear();
                aFluidMesh.clear();
                UFluidMesh.clear();
            }
        }

        void poroSolid::afterFluidSolve()
        {
            q_relAcc_.clear();
        }

        void poroSolid::clearCouplingTerms()
        {
            q_relAcc_.clear();
            nDot_.clear();
            fixedStressStabil_.clear();
        }

        //- Darcy-type flux induced by solid-matrix acceleration.
        //  With hydraulic conductivity K [m/s], the reduced saturated
        //  fluid-momentum balance gives q_acc = -(K/|g|) a_s in addition to the
        //  standard pressure-driven Darcy flux. This helper returns the
        //  positive face-vector quantity K/|g| a_s; explicitCouplingDtoP()
        //  later inserts it as -div(q_relAcc), so the assembled flux
        //  contribution is the required -q_relAcc.
        tmp<surfaceVectorField> poroSolid::q_relAcc(const surfaceScalarField& kf, const volVectorField& a)
        {
            return poroCouplingTerms::relativeAccelerationFlux(kf, a, mag(poroFluid().g()));
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        poroSolid::poroSolid(
            Time &runTime,
            const word &region)
            : poroSolidInterface(typeName, runTime, region),
              fixedStressStabil_(), // Pointer to stabilization term
              nDot_(), // Pointer to porosity change per unit time (exchange term)
              q_relAcc_(), // Pointer to rel. Acceleration term (exchange term)
              pRghSolidMesh_(), // p_rgh on solid Mesh (used in buoyancy less mode)
              pSolidMesh_() // p on solid Mesh (includes the buoyancy force on the solid)
        {
            // Sanity check if an appropriate poroSolid solver is used for the selected poroFluid solver
            if(poroFluid().name() == "varSatPoroFluid" || poroFluid().name() == "varSatPoroFluidHead")
            {
                Warning() << "'varSatPoroFluid/Head' should be used with 'varSatPoroSolid'" << endl;
            }
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        //- explicit coupling terms to pressure equation
        const tmp<volScalarField> poroSolid::explicitCouplingDtoP() const
        {
                if(!nDot_.valid())
                {
                    FatalErrorInFunction
                        << "Explicit coupling term nDot is not initialized"
                        << exit(FatalError);
                }
                if(q_relAcc_.valid())
                {
                    // See q_relAcc() above: the physical acceleration-induced
                    // Darcy contribution is -q_relAcc, hence the source enters
                    // continuity as -div(q_relAcc).
                    return tmp<volScalarField>
                    (
                        poroCouplingTerms::explicitCouplingSource(nDot_(), q_relAcc_()).ptr()
                    );
                }

                return tmp<volScalarField>(new volScalarField(nDot_()));
            }
            //- implicit coupling terms to pressure equation (get multiplied with pField)
            const tmp<volScalarField>  poroSolid::implicitCouplingDtoP() const
            {
                if(!fixedStressStabil_.valid())
                {
                    FatalErrorInFunction
                        << "Implicit coupling term fixedStressStabil is not initialized"
                        << exit(FatalError);
                }
            tmp<volScalarField> tSp(
                new volScalarField(fixedStressStabil_())
                );
            return tSp;
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace poroSolidInteractions

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
