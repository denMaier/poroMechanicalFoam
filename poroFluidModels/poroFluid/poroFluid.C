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

#include "poroFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "iterationControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

        namespace poroFluidModels
        {
            // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

            defineTypeNameAndDebug(poroFluid, 0);
            addToRunTimeSelectionTable(poroFluidModel, poroFluid, dictionary);
            

            // * * * * * * * * * * * * * Private Memeber Functions * * * * * * * * * * * //

            poroHydraulicModel &poroFluid::poroHydraulic()
            {
                if (!poroHydPtr_.valid())
                {
                    poroHydPtr_.reset
                    (
                        new poroHydraulicModel(p_rgh(), g())
                    );
                }

                return poroHydPtr_.ref();
            }

            // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

            poroFluid::poroFluid(
                Time &runTime,
                const word &region,
                const bool sharedmesh)
                : poroFluidModel(typeName, runTime, "p_rgh", region, sharedmesh),
                  poroHydPtr_(),
                  Ss_
                  (
                    (
                      IOobject
                      (
                          "Ss", //Call it effective to make it easier for BC to find it.
                          runTime.timeName(),
                          p_rgh().db(),
                          IOobject::NO_READ,
                          IOobject::AUTO_WRITE
                      ),
                      poroHydraulic().Ss(n(),p_rgh())
                    )
                  ),
                  kbyGammaf_
                  (
                      IOobject
                      (
                          "kEffbyGammaf", //Call it effective to make it easier for BC to find it.
                          runTime.timeName(),
                          p_rgh().db(),
                          IOobject::NO_READ,
                          IOobject::AUTO_WRITE
                      ),
                      poroHydraulic().kf() / poroHydraulic().magGamma()
                  )
            {
                pRGHisRequired();
                //Info << runTime().constant() << nl << runTime().caseConstant() <<endl;
                phi() = (-kbyGammaf_ * mesh().magSf() * fvc::snGrad(p_rgh()));
                Info << "Done generating poroFluidModel" << endl;

                if(debug)
                {
                    n().write();
                    const tmp<volScalarField> tk(poroHydraulic().k());
                    volScalarField kDiagnostic(tk());
                    kDiagnostic.writeOpt() = IOobject::AUTO_WRITE;
                    kDiagnostic.write();
                    Ss_.write();
                }
            }

            // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

            bool poroFluid::evolve()
            {

                // Initialize performance tracking

                SolverPerformance<scalar> solverPerfp;
                SolverPerformance<scalar>::debug=0;

	            iterCtrl().reset();

                Info << "Evolving fluid model: " << this->type() << endl;

                do // NonOrhtogonal Correctors
                {
                    // Store prev fields for relaxation
                    p_rgh().storePrevIter();

                    if (debug)
                    {
                        Info << "Assembling pEqn"
                             << endl;
                    }

                    if(poroHydraulic().updatesSs())
                    {
                        Ss_ = poroHydraulic().Ss(n(),p_rgh());
                    }

                    fvScalarMatrix pEqn(
                        fvm::ddt(Ss_, p_rgh())                //- Storage
                        - fvm::laplacian(kbyGammaf_, p_rgh())  //- Pressure Fluxes
                        == fvOptions()(Ss_,p_rgh())
                    );
                    
                    fvOptions().constrain(pEqn);

                    solverPerfp = pEqn.solve();

                    p_rgh().relax();
                
                    fvOptions().correct(p_rgh());
                    
                    phi() = pEqn.flux();
                    hydraulicGradient() = fvc::grad(p_rgh())/poroHydraulic().magGamma();

                    if(poroHydraulic().updatesK())
                    {
                        kbyGammaf_ = poroHydraulic().kf() / poroHydraulic().magGamma();
                    }

                    p_rgh().correctBoundaryConditions();

                } while (iterCtrl().loop());

                iterCtrl().write();

                Info << "poroFluid evolved" << endl;

                p() = p_rgh()+poroHydraulic().p_Hyd();

                pDot()= fvc::ddt(p());
                return 0;
            }

            scalar poroFluid::newDeltaT()
            {
                const tmp<volScalarField> tDenom
                (
                    Ss_ * poroHydraulic().magGamma()
                );
                const dimensionedScalar minDenom
                (
                    "minDenom",
                    tDenom().dimensions(),
                    SMALL
                );
                const scalar maxCv
                (
                    max
                    (
                        poroHydraulic().k()
                      / max(tDenom(), minDenom)
                    ).value()
                );
                Info << "max c_v: " << maxCv << endl;

                if (maxCv <= SMALL)
                {
                    WarningInFunction
                        << "Degenerate hydraulic diffusivity estimate detected. "
                        << "Keeping current deltaT = " << runTime().deltaTValue()
                        << endl;
                    return runTime().deltaTValue();
                }

                const scalar minDeltaX
                (
                    min(mesh().V() / fvc::surfaceSum(mag(mesh().Sf())))
                );
                const scalar wishedTimeStep
                (
                    CoNumber_ * pow(minDeltaX, 2) / maxCv
                );
                const scalar maxDeltaTFact
                (
                    wishedTimeStep / max(runTime().deltaTValue(), SMALL)
                );
                const scalar deltaTFact
                (
                    min(min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2)
                );

                return deltaTFact * runTime().deltaTValue();
            }

            void poroFluid::writeFields(const Time &runTime)
            {
                poroFluidModel::writeFields(runTime);
            }

            void poroFluid::end()
            {
                poroHydraulic().write();
                poroFluidModel::end();
            }

            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        } // End namespace poroFluidModels
} // End namespace Foam

// ************************************************************************* //
