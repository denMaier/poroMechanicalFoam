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
    variablySaturatedPoroFluid

Description
    Abstract intermediate base for variably saturated pore-fluid solvers based on
    Richards equation.

    The class factors out everything shared by the concrete `p_rgh` and
    `pHead` formulations:
    - saturation-dependent constitutive data via `varSatPoroHydraulicModel`
    - the pore-fluid saturation field `S`
    - optional mass-balance residual bookkeeping
    - runtime selection of the nonlinear/linearization strategy through
      `richardsLinearization`

    The model consistently treats `S` as pore-fluid saturation, i.e. the
    occupied fraction of pore volume by the mobile pore fluid. This is distinct
    from pure water saturation when the constitutive model accounts for trapped
    gas or similar effects. As a consequence, fluxes and storage terms are based
    on effective hydraulic properties rather than simply saturated values.

    Concrete subclasses only need to decide which pressure variable is solved
    (`p_rgh` or `pHead`) and how the corresponding equation is assembled.

Author
    Denis Maier, BAW.  All rights reserved.

SourceFiles
    variablySaturatedPoroFluid.C

\*---------------------------------------------------------------------------*/

#ifndef variablySaturatedPoroFluid_H
#define variablySaturatedPoroFluid_H

#include "poroFluidModel.H"
#include "varSatPoroHydraulicModel.H"
#include "richardsLinearization.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace poroFluidModels
    {
        /*---------------------------------------------------------------------------*\
                        Class variablySaturatedPoroFluid Declaration
\*---------------------------------------------------------------------------*/

        /**
         * Shared base class for Richards-equation pore-fluid solvers.
         *
         * `variablySaturatedPoroFluid` sits between the generic
         * `poroFluidModel` interface and the concrete variably saturated
         * formulations. It owns the constitutive model, the saturation field and
         * the linearization strategy, while the children implement the actual
         * equation form and choose the primary field exposed through `pField()`.
         */
        class variablySaturatedPoroFluid
            : public poroFluidModel
        {
            // Private data

            //- Disallow default bitwise copy construct
            variablySaturatedPoroFluid(const variablySaturatedPoroFluid &);

            //- Disallow default bitwise assignment
            void operator=(const variablySaturatedPoroFluid &);

        protected:
            // Protected Data

            //- Variably saturated constitutive model for saturation, storage,
            //  conductivity and hydrostatic reconstruction.
            mutable autoPtr<varSatPoroHydraulicModel> poroHydPtr_;

            //- Runtime-selected linearization/nonlinear update strategy for
            //  Richards equation.
            mutable autoPtr<richardsLinearization> richardsLinearizationPtr_;

            //- Saturation of the mobile pore fluid occupying the pore space.
            autoPtr<volScalarField> SPtr_;

            //- Optional field storing a mass-balance residual for diagnostics.
            autoPtr<volScalarField> MassBalancePtr_;

            // Protected Member Functions
            //- Construct the concrete constitutive model for the active child.
            virtual void makePoroHydraulic() = 0;

            //- Allocate the saturation field on first use.
            void makeS();

            //- Copy a fully assembled saturation field into the registered S
            //  object without renaming or re-registering it in the registry.
            void setS(const tmp<volScalarField>& sat)
            {
                volScalarField& SRef = S();
                const volScalarField& satRef = sat();

                SRef.internalFieldRef() = satRef.internalField();

                forAll(SRef.boundaryFieldRef(), patchI)
                {
                    SRef.boundaryFieldRef()[patchI] = satRef.boundaryField()[patchI];
                }

                SRef.correctBoundaryConditions();
            }

            //- Recompute saturation from pressure and write it into the registered
            //- saturation field.
            void updateS(const volScalarField& pressure)
            {
                setS(poroHydraulic().S(pressure));
            }

            //- Construct the runtime-selected Richards linearization scheme.
            void makeRichardsLinearization();

            //- Non-const access to the lazily created pore-fluid saturation.
            volScalarField& S()
            {
                if(!SPtr_.valid())
                {   
                    makeS();
                }
                return SPtr_.ref();
            }

            //- Non-const access to the lazily created constitutive model.
            varSatPoroHydraulicModel& poroHydraulic()
            {
                if(!poroHydPtr_.valid())
                {   
                    makePoroHydraulic();
                }
                return poroHydPtr_.ref();
            }

            //- Non-const access to the lazily created linearization object.
            richardsLinearization& linearization()
            {
                if(!richardsLinearizationPtr_.valid())
                {   
                    makeRichardsLinearization();
                }
                return richardsLinearizationPtr_.ref();
            }


        public:
            //- Runtime type information
            TypeName("variablySaturatedPoroFluid");
            // Declare run-time constructor selection table

            // Constructors

            //- Construct from components
            variablySaturatedPoroFluid(
                const word &type,
                Time &runTime,
                const word &fieldName,
                const word &region = dynamicFvMesh::defaultRegion,
                const bool sharedMesh = false);

            // Destructor

            virtual ~variablySaturatedPoroFluid()
            {
            }

            // Selector

            // Member Functions

            // Special Acess

            //- Const access to the pore-fluid saturation field.
            virtual const volScalarField& S() const override
            {
                if (!SPtr_.valid())
                {
                    const_cast<variablySaturatedPoroFluid&>(*this).makeS();
                }

                return SPtr_();
            }

            //- Return the reference porosity from the constitutive model.
            virtual tmp<volScalarField> n0()
            {
                return poroHydraulic().n0();
            }

            //- Const access to the variably saturated constitutive model.
            const varSatPoroHydraulicModel& poroHydraulic() const
            {
                if (!poroHydPtr_.valid())
                {
                    const_cast<variablySaturatedPoroFluid&>(*this).makePoroHydraulic();
                }

                return poroHydPtr_();
            }

            // Edit

            // Reminder: This still needs to be implemented by children
            //- Evolve the fluid model
            //virtual bool evolve() = 0;

            // Reminder: This still needs to be implemented by children
            //- Return the desired new time-step size
            //virtual scalar newDeltaT() = 0;

            //- Allocate diagnostics needed for mass-balance monitoring.
            void checkMassBalance();

            //- Return the accumulated material residual if diagnostics are active.
            virtual scalar materialResidual()
            {
                return MassBalancePtr_.valid()
                           ? gSum(MassBalancePtr_())
                           : 0.0;
            }

            //- Return number of top-level iterations.
            virtual label outerIters()
            {
                notImplemented("outerIters()");
                return iterCtrl().index(); // for now, to keep the compiler happy
            }

            // Reminder: This still needs to be implemented by children
            //virtual void writeFields(const Time &runTime) = 0;

            //- Persist hydraulic-model state before the simulation ends.
            virtual void end();
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace variablySaturatedPoroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
