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
    poroSolid

Description
    Coupling implementation for fully saturated poromechanics.

    This class specializes Foam::poroSolidInterface for the saturated case,
    where pore pressure alone drives the hydraulic-to-mechanical coupling.

    Responsibilities:
    - build the deformation-to-flow source terms used by
      Foam::fv::poroSolidToFluidCouplingSource
    - expose pore pressure fields on the solid side, either by registering the
      fluid fields directly for shared meshes or by storing mapped copies on
      the solid mesh
    - execute the staggered pressure-displacement outer loop:
      1. assemble coupling terms from the current solid state
      2. evolve the poroFluid model
      3. map hydraulic fields back to the solid mesh if required
      4. evolve the solid model

    The flow-to-solid coupling itself is consumed by the porous mechanical laws
    on the solid side, which look up `p` and `p_rgh` in the solid mesh object
    registry.

Author
    Denis Maier, BAW. All rights reserved.
    Héctor Mentenegro, BAW. All rights reserved.

SourceFiles
    poroSolid.C

\*---------------------------------------------------------------------------*/

#ifndef poroSolid_H
#define poroSolid_H

#include "poroSolidInterface.H"
#include "poroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace poroSolidInteractions
    {

        /*---------------------------------------------------------------------------*\
                        Class poroSolid Declaration
        \*---------------------------------------------------------------------------*/

        /**
         * \brief Fully saturated poro-solid coupling driver.
         *
         * The class extends Foam::poroSolidInterface with the fields required
         * by the saturated formulation. It owns transient coupling fields that
         * are assembled before each fluid solve and exposes them through
         * explicitCouplingDtoP() and implicitCouplingDtoP().
         *
         * For shared-mesh runs, the fluid pressure fields are registered
         * directly on the solid mesh registry. For non-shared runs, mapped
         * copies are stored in pSolidMesh_ and pRghSolidMesh_ so the solid-side
         * constitutive laws see the same names independent of mesh layout.
         */
        class poroSolid
            : public poroSolidInterface
        {
            // Private data

            //- Fixed-stress stabilization coefficient for the flow equation.
            autoPtr<volScalarField> fixedStressStabil_;
            //- Volumetric strain-rate contribution added explicitly to flow.
            autoPtr<volScalarField> nDot_;
            //- Relative acceleration correction for Darcy fluxes (currently unused).
            autoPtr<surfaceVectorField> q_relAcc_;
            //- Pressure head on the solid mesh for non-shared mesh coupling.
            autoPtr<volScalarField> pRghSolidMesh_;
            //- Absolute pore pressure on the solid mesh for non-shared mesh coupling.
            autoPtr<volScalarField> pSolidMesh_;

            /**
             * \brief Estimate the additional flux caused by relative acceleration.
             *
             * The term is kept for completeness but is currently not assembled
             * into the active coupling operators.
             */
            tmp<surfaceVectorField> q_relAcc(const surfaceScalarField& kf, const volVectorField& a);

            //- Disallow default bitwise copy construct
            poroSolid(const poroSolid &);

            //- Disallow default bitwise assignment
            void operator=(const poroSolid &);

        protected:
            // Protected member functions

            virtual void prepareCouplingLoop();

            //- Lazily register or allocate hydraulic fields required on the solid side.
            virtual void initializeSolidHydraulicFields();

            //- Assemble saturated coupling terms on the fluid mesh.
            virtual void assembleCouplingTerms();

            //- Synchronize the hydraulic state required on the solid side.
            virtual void syncSolidHydraulicFields();

            //- Cleanup per-iteration transient terms after the fluid solve.
            virtual void afterFluidSolve();

            //- Cleanup transient coupling fields after the outer loop.
            virtual void clearCouplingTerms();

        public:
            //- Runtime type information
            TypeName("poroSolid");

            // Constructors

            /**
             * \brief Construct the saturated poro-solid coupling model.
             *
             * Besides base-class construction, this prepares the saturated
             * coupling object. The solid-side hydraulic fields are then
             * registered or allocated lazily before the first coupled solve.
             */
            poroSolid(
                Time &runTime,
                const word &region = "poroFluid");

            // Destructor

            virtual ~poroSolid()
            {
            }

            // Member Functions

            // Edit

            /**
             * \brief Advance one coupled pressure-displacement time step.
             *
             * The method performs the staggered outer iteration managed by
             * couplingControl(). It assembles deformation-to-flow terms,
             * evolves the poroFluid model, maps hydraulic fields back to the
             * solid side if necessary, evolves the solid model, and finally
             * updates porosity for the completed time step.
             */
            // Access

            /**
             * \brief Return the explicit deformation source for the flow solve.
             *
             * The returned field is the currently assembled `nDot` term, i.e.
             * the Biot-weighted volumetric strain-rate contribution.
             */
            virtual const tmp<volScalarField>  explicitCouplingDtoP() const;

            /**
             * \brief Return the implicit fixed-stress stabilization coefficient.
             *
             * The fvOption multiplies this field with the current pressure field
             * when assembling the stabilized flow equation.
             */
            virtual const tmp<volScalarField>  implicitCouplingDtoP() const; 

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // namespace solidModels

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
