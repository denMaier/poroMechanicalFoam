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

#include "varSatPoroSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "mechanicalModel.H"
#include "iterationControl.H"
#include "poroCouplingTerms.H"
#include "poroPressureUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace poroSolidInteractions
    {
        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(varSatPoroSolid, 0);
        addToRunTimeSelectionTable(poroSolidInterface, varSatPoroSolid, dictionary);

        // * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //
        void varSatPoroSolid::initializePressureUnits()
        {
            if(magGammaWInitialized_)
            {
                return;
            }

            magGammaW_ = poroPressureUnits::pressureScale(pField().dimensions(), poroFluid());

            magGammaWInitialized_ = true;
        }

        const volScalarField& varSatPoroSolid::saturationField() const
        {
            if(!poroFluid().foundObject<volScalarField>("S"))
            {
                FatalErrorInFunction
                    << "variably saturated coupling requires a saturation field 'S' "
                    << "once saturation is actually used"
                    << exit(FatalError);
            }

            return poroFluid().lookupObject<volScalarField>("S");
        }

        void varSatPoroSolid::initializeSolidHydraulicFields()
        {
            if(sharedMesh())
            {
                ensureSharedSolidFieldRegistered(poroFluidRef().p(), poroFluidMesh().name());
                ensureSharedSolidFieldRegistered(poroFluidRef().p_rgh(), poroFluidMesh().name());
                ensureSharedSolidFieldRegistered(const_cast<volScalarField&>(saturationField()), poroFluidMesh().name());
                ensureSharedSolidFieldRegistered(const_cast<volScalarField&>(poroFluidRef().n()), poroFluidMesh().name());
            }
            else
            {
                const volScalarField& S = saturationField();

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

                if(!SSolidMesh_.valid())
                {
                    SSolidMesh_.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "S",
                                runTime().timeName(),
                                solidMesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            solidToPoroFluid().mapSrcToTgt(S)(),
                            "zeroGradient"
                        )
                    );
                }

                if(!nSolidMesh_.valid())
                {
                    nSolidMesh_.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "n",
                                runTime().timeName(),
                                solidMesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            solidToPoroFluid().mapSrcToTgt(poroFluidRef().n())(),
                            "zeroGradient"
                        )
                    );
                }
            }
        }

        void varSatPoroSolid::syncSolidHydraulicFields()
        {
            if(sharedMesh())
            {
                return;
            }

            initializeSolidHydraulicFields();
            mapPressuresToSolidMesh(pSolidMesh_, pRghSolidMesh_);
            SSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(saturationField());
        }

        void varSatPoroSolid::syncSolidPorosityField()
        {
            if(sharedMesh())
            {
                return;
            }

            initializeSolidHydraulicFields();
            nSolidMesh_.ref() = solidToPoroFluid().mapSrcToTgt(poroFluidRef().n())();
        }

        void varSatPoroSolid::prepareCouplingLoop()
        {
            Info << "Preparing varSatPoroSolid solver" << endl;
            initializePressureUnits();
        }

        void varSatPoroSolid::assembleCouplingTerms()
        {
            const volScalarField& SFluidMesh = saturationField();

            if(sharedMesh())
            {
                const volScalarField impK(solid().mechanical().impK()/magGammaW_);
                tmp<volScalarField> bishopBiot(SFluidMesh*b());
                const tmp<surfaceScalarField> tkf(poroFluid().relativeAccelerationConductivity());
                const tmp<volVectorField> tSolidA(fvc::d2dt2(solid().D()));

                // The explicit deformation-to-flow term is saturation-weighted,
                // but the fixed-stress stabilization intentionally keeps the
                // canonical Biot-based form through updateCouplingTerms().
                // This term is used for split-scheme robustness, not as a
                // direct physical coupling law, so we do not weaken it by S.
                updateCouplingTerms(bishopBiot(), impK, solid().U(), nDot_, fixedStressStabil_);
                q_relAcc_.reset(q_relAcc(tkf(), tSolidA()).ptr());
            }
            else
            {
                Info << "Mapping fields to poroFluid mesh" << endl;
                const volScalarField impK(solid().mechanical().impK()/magGammaW_);
                const tmp<surfaceScalarField> tkf(poroFluid().relativeAccelerationConductivity());
                tmp<volVectorField> UFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().U());
                const tmp<volVectorField> tSolidA(fvc::d2dt2(solid().D()));
                tmp<volVectorField> aFluidMesh = solidToPoroFluid().mapTgtToSrc(tSolidA());
                tmp<volScalarField> tmpImpK(solidToPoroFluid().mapTgtToSrc(impK)());
                tmp<volScalarField> bishopBiot(SFluidMesh*b());

                // The explicit deformation-to-flow term is saturation-weighted,
                // but the fixed-stress stabilization intentionally keeps the
                // canonical Biot-based form through updateCouplingTerms().
                // This term is used for split-scheme robustness, not as a
                // direct physical coupling law, so we do not weaken it by S.
                updateCouplingTerms(bishopBiot(), tmpImpK(), UFluidMesh(), nDot_, fixedStressStabil_);
                q_relAcc_.reset(q_relAcc(tkf(), aFluidMesh()).ptr());

                tmpImpK.clear();
                aFluidMesh.clear();
                UFluidMesh.clear();
                bishopBiot.clear();
            }
        }

        void varSatPoroSolid::afterPorosityUpdate()
        {
            syncSolidPorosityField();
        }

        void varSatPoroSolid::afterFluidSolve()
        {
            q_relAcc_.clear();
        }

        void varSatPoroSolid::beforeSolidSolve()
        {
            if(couplingControl().index() > 0)
            {
                solidRef().recalculateRho();
            }
        }

        void varSatPoroSolid::clearCouplingTerms()
        {
            q_relAcc_.clear();
            nDot_.clear();
            fixedStressStabil_.clear();
        }

        void varSatPoroSolid::writeAdditionalFields(const Time&)
        {
            solid().rho().write();
        }

        //- Darcy-type flux induced by solid-matrix acceleration.
        //  For the liquid phase the reduced momentum balance gives
        //  q_acc = -(K_eff/|g|) a_s, where K_eff already contains the
        //  saturation-dependent mobility. This helper returns the positive
        //  face-vector quantity K_eff/|g| a_s; explicitCouplingDtoP() later
        //  inserts it as -div(q_relAcc), so the assembled flux contribution is
        //  the required -q_relAcc.
        tmp<surfaceVectorField> varSatPoroSolid::q_relAcc(const surfaceScalarField& kf, const volVectorField& a)
        {
            return poroCouplingTerms::relativeAccelerationFlux(kf, a, mag(poroFluid().g()));
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
        varSatPoroSolid::varSatPoroSolid(
            Time &runTime,
            const word &region)
            : poroSolidInterface(typeName, runTime, region),
              fixedStressStabil_(), // Pointer to stabilization term
              nDot_(), // Pointer to porosity change per unit time (exchange term)
              q_relAcc_(), // Pointer to rel. Acceleration term (exchange term)
              pRghSolidMesh_(), // p_rgh on solid Mesh (used in buoyancy less mode)
              pSolidMesh_(), // p on solid Mesh (includes the buoyancy force on the solid)
              SSolidMesh_(), // saturation on solid Mesh 
              nSolidMesh_(), // porosity on solid Mesh 
              magGammaW_(dimensionedScalar("gammaw",dimless,1.0)), // this to facilitate pressureHead and p_rgh based solvers in one class
              magGammaWInitialized_(false)
        {
            // Sanity check if an appropriate poroSolid solver is used for the selected poroFluid solver
            if(poroFluid().name() == "poroFluid")
            {
                FatalErrorInFunction()
                    << "The selected poroFluid model '" << poroFluid().name()
                    << "' is not suitable for variably saturated calculations." << nl
                    << "Use a variably saturated poroFluid model such as varSatPoroFluid or varSatPoroFluidHead."
                    << exit(FatalError);
            }
        }

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        //- explicit coupling terms to pressure equation
        const tmp<volScalarField> varSatPoroSolid::explicitCouplingDtoP() const
        {
            if(!nDot_.valid())
            {
                FatalErrorInFunction
                    << "Explicit coupling term nDot is not initialized"
                    << exit(FatalError);
            }
            if(q_relAcc_.valid())
            {
                // See q_relAcc() above: the physical acceleration-induced Darcy
                // contribution is -q_relAcc, hence the source enters continuity
                // as -div(q_relAcc).
                return tmp<volScalarField>
                (
                    poroCouplingTerms::explicitCouplingSource(nDot_(), q_relAcc_()).ptr()
                );
            }

            return tmp<volScalarField>(new volScalarField(nDot_()));
        }
        //- implicit coupling terms to pressure equation (get multiplied with pField)
        const tmp<volScalarField>  varSatPoroSolid::implicitCouplingDtoP() const
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
