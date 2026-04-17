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

Class
    poroFluid

Description
    Saturated groundwater-flow solver using the excess-pressure formulation
    `p_rgh`, i.e. pressure relative to the initial hydrostatic state.

    This is the simplest concrete `poroFluidModel` implementation in the
    hierarchy. It solves a linear diffusion/storage equation for
    single-phase saturated pore-fluid flow and delegates constitutive data to
    `poroHydraulicModel`.

    Main responsibilities:
    - keep `p_rgh` as the primary unknown exposed through `pField()`
    - assemble the linear storage + Darcy-flux equation
    - update the total pore pressure `p = p_rgh + p_Hyd`
    - publish flux and hydraulic-gradient fields for coupling and boundary
      conditions

    The solver is selected at run time via the `poroFluidModel` factory and is
    therefore the saturated counterpart to the variably saturated Richards
    solvers in `variableSaturated/`.

Author
    Denis Maier, BAW.   All rights reserved.

SourceFiles
    poroFluid.C

\*---------------------------------------------------------------------------*/

#ifndef poroFluid_H
#define poroFluid_H

#include "poroFluidModel.H"
#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

        namespace poroFluidModels
        {
            /*---------------------------------------------------------------------------*\
                                    Class poroFluid Declaration
            \*---------------------------------------------------------------------------*/

            /**
             * Saturated pore-fluid solver based on excess pressure `p_rgh`.
             *
             * The class extends `poroFluidModel` with a lazily owned
             * `poroHydraulicModel`, a storage field `Ss_`, and a face-based
             * conductivity field `kbyGammaf_`. Together these fields define the
             * linear saturated-flow equation solved in `evolve()`.
             *
             * Interface notes:
             * - `pField()` aliases `p_rgh()` so generic coupling code can work
             *   with the active primary unknown.
             * - `p()` stores the reconstructed total pore pressure required by
             *   post-processing and coupled solid mechanics.
             * - `phi()` and `hydraulicGradient()` are refreshed every solve for
             *   downstream boundary conditions and source models.
             */
            class poroFluid
                : public poroFluidModel
            {
                
                //- Constitutive model for saturated conductivity, storage and
                //  hydrostatic pressure reconstruction.
                mutable autoPtr<poroHydraulicModel> poroHydPtr_;

                //- Specific storage evaluated from the hydraulic model.
                volScalarField Ss_;

                //- Face-interpolated hydraulic conductivity divided by the
                //  specific weight of the pore fluid, used in the Darcy term.
                surfaceScalarField kbyGammaf_;

                // Private Member Functions
                //- Return the initialized hydraulic model instance.
                poroHydraulicModel& poroHydraulic();

                //- Disallow default bitwise copy construct
                poroFluid(const poroFluid &);

                //- Disallow default bitwise assignment
                void operator=(const poroFluid &);

            public:
                //- Runtime type information
                TypeName("poroFluid");
                // Declare run-time constructor selection table

                // Constructors

                //- Construct from components
                poroFluid(
                    Time &runTime,
                    const word &region = dynamicFvMesh::defaultRegion,
                    const bool sharedMesh = false);

                // Destructor

                virtual ~poroFluid()
                {
                }

                // Selector

                // Member Functions

                // Access
                
                //- Alias the generic primary field handle to the excess
                //  pressure solved by this class.
                virtual volScalarField& pField()
                {
                    return p_rgh();
                }
                
                //- Const overload of the `p_rgh` primary-field alias.
                virtual const volScalarField& pField() const
                {
                    return p_rgh();
                }

                //- Return the reference porosity field from the hydraulic model.
                virtual tmp<volScalarField> n0()
                {
                    return poroHydraulic().n0();
                }

                // Special Access
                //- Const access to the saturated hydraulic constitutive model.
                const poroHydraulicModel& poroHydraulic() const
                {
                    return poroHydPtr_();
                }

                // Edit

                //- Advance the saturated-flow equation for one time step.
                virtual bool evolve();

                //- Estimate a stable next time step from a diffusivity scale.
                virtual scalar newDeltaT();

                //- No separate nonlinear material residual is required because
                //  the governing equation is linear in the solved state.
                virtual scalar materialResidual()
                {
                    return 0.0;
                }

                //- Report the outer iteration counter managed by iterationControl.
                virtual label outerIters()
                {
                    return iterCtrl().index();
                }

                //- Delegate generic field writing to the common base class.
                virtual void writeFields(const Time &runTime);

                //- Persist hydraulic-model state before the simulation ends.
                virtual void end();
            };

            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        } // End namespace poroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
