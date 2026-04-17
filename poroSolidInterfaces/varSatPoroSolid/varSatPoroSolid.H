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
    varSatPoroSolid

Description
    Coupling implementation for variably saturated poromechanics.

    This specialization augments the saturated coupling workflow with the
    additional state required by the partially saturated formulation:
    saturation, porosity on the solid mesh, and pressure-unit conversion for
    pressure-head based fluid solvers.

    The class follows the same staggered outer loop as Foam::poroSolid, but it
    scales the deformation-to-flow coupling with saturation (Bishop-type
    weighting) and keeps the solid-side copies of `S` and `n` synchronized when
    separate meshes are used.

Author
    Denis Maier, BAW. All rights reserved.
    Héctor Mentenegro, BAW. All rights reserved.

SourceFiles
    varSatPoroSolid.C

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroSolid_H
#define varSatPoroSolid_H

#include "poroSolidInterface.H"
#include "varSatPoroHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    namespace poroSolidInteractions
    {

        /*---------------------------------------------------------------------------*\
                        Class varSatPoroSolid Declaration
        \*---------------------------------------------------------------------------*/

        /**
         * \brief Variably saturated poro-solid coupling driver.
         *
         * The class extends Foam::poroSolidInterface with the extra fields
         * required by the variably saturated formulation:
         * - mapped or shared access to `p` and `p_rgh`
         * - mapped or shared access to saturation `S`
         * - mapped or shared access to porosity `n`
         * - a conversion factor between pressure-head and pressure units
         *
         * These fields allow the solid constitutive model and density update
         * logic to consume a consistent hydraulic state regardless of whether
         * the fluid and solid run on the same mesh.
         */
        class varSatPoroSolid
            : public poroSolidInterface
        {
            // Private data

            //- Fixed-stress stabilization coefficient for the flow equation.
            autoPtr<volScalarField> fixedStressStabil_;
            //- Saturation-weighted volumetric strain-rate source for flow.
            autoPtr<volScalarField> nDot_;
            //- Relative acceleration correction for Darcy fluxes (currently unused).
            autoPtr<surfaceVectorField> q_relAcc_;
            //- Pressure head on the solid mesh for non-shared mesh coupling.
            autoPtr<volScalarField> pRghSolidMesh_;
            //- Absolute pore pressure on the solid mesh for non-shared mesh coupling.
            autoPtr<volScalarField> pSolidMesh_;
            //- Saturation field mirrored onto the solid mesh.
            autoPtr<volScalarField> SSolidMesh_;
            //- Porosity field mirrored onto the solid mesh for density updates.
            autoPtr<volScalarField> nSolidMesh_;
            //- Pressure-head to pressure conversion factor for head-based solvers.
            dimensionedScalar magGammaW_;
            //- Has the pressure-head conversion factor been initialized?
            bool magGammaWInitialized_;

            //- Lazily initialize the pressure-head conversion factor when needed.
            void initializePressureUnits();

            //- Return the saturation field once it actually exists in the poroFluid registry.
            const volScalarField& saturationField() const;

            //- Synchronize the porosity field required on the solid side.
            void syncSolidPorosityField();

            /**
             * \brief Estimate the additional flux caused by relative acceleration.
             *
             * The term is currently kept as an optional extension path and is
             * not assembled into the active coupling operators.
             */
            tmp<surfaceVectorField> q_relAcc(const surfaceScalarField& kf, const volVectorField& U);

            // Private Member Functions

            //- Disallow default bitwise copy construct
            varSatPoroSolid(const varSatPoroSolid &);

            //- Disallow default bitwise assignment
            void operator=(const varSatPoroSolid &);

        protected:

            virtual void prepareCouplingLoop();

            //- Lazily register or allocate hydraulic fields required on the solid side.
            virtual void initializeSolidHydraulicFields();

            //- Assemble variably saturated coupling terms on the fluid mesh.
            virtual void assembleCouplingTerms();

            //- Synchronize pressure and saturation fields required on the solid side.
            virtual void syncSolidHydraulicFields();

            //- Keep mirrored solid-side porosity in sync after porosity updates.
            virtual void afterPorosityUpdate();

            //- Cleanup per-iteration transient terms after the fluid solve.
            virtual void afterFluidSolve();

            //- Update solid-side density state after fluid-side hydraulic refresh.
            virtual void beforeSolidSolve();

            //- Cleanup transient coupling fields after the outer loop.
            virtual void clearCouplingTerms();

            //- Write solid density in addition to the standard outputs.
            virtual void writeAdditionalFields(const Time& runTime);

        public:
            //- Runtime type information
            TypeName("varSatPoroSolid");

            // Constructors

            /**
             * \brief Construct the variably saturated poro-solid coupling model.
             *
             * Construction validates that a variably saturated poroFluid model
             * is in use and prepares the lazy initialization paths for the
             * solid-side hydraulic state and any pressure-head conversion.
             */
            varSatPoroSolid(
                Time &runTime,
                const word &region = "poroFluid");

            // Destructor

            virtual ~varSatPoroSolid()
            {
            }

            // Member Functions

            // Edit

            /**
             * \brief Advance one coupled pressure-displacement time step.
             *
             * Compared to the saturated model, the deformation-to-flow source is
             * weighted by saturation and the mirrored solid-side saturation and
             * porosity fields are refreshed during the coupling loop so the
             * constitutive law and density updates see the latest hydraulic
             * state.
             */
            // Access

            /**
             * \brief Return the explicit deformation source for the flow solve.
             *
             * This is the assembled `nDot` field based on the saturation-
             * weighted Biot coefficient and the current solid velocity field.
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
