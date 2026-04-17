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
    varSatPoroFluid

Description
    Richards-equation solver for variably saturated flow using the excess
    pressure `p_rgh` as the primary unknown.

    Relative to `varSatPoroFluidHead`, this formulation absorbs the hydrostatic
    reference state into the field decomposition and solves only for deviations
    from that state. The class therefore reconstructs the total pore pressure
    as `p = p_rgh + p_Hyd` and evaluates saturation-dependent constitutive data
    from that total pressure.

    The solver combines:
    - `variablySaturatedPoroFluid` for shared Richards-equation infrastructure
    - `varSatPoroHydraulicModel` for constitutive laws
    - a runtime-selected `richardsLinearization` scheme for nonlinear updates

Author
    Denis Maier, BAW.  All rights reserved.

SourceFiles
    varSatPoroFluid.C

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroFluid_H
#define varSatPoroFluid_H

#include "variablySaturatedPoroFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace poroFluidModels
    {
        /*---------------------------------------------------------------------------*\
                        Class varSatPoroFluid Declaration
\*---------------------------------------------------------------------------*/

        /**
         * Variably saturated Richards solver in excess-pressure form.
         *
         * This concrete class specializes `variablySaturatedPoroFluid` by
         * choosing `p_rgh` as the primary unknown and by assembling the
         * corresponding storage, saturation and Darcy-flux terms in `evolve()`.
         */
        class varSatPoroFluid
            : public variablySaturatedPoroFluid
        {
            // Private data

            //- Compressibility / storage contribution used in the pressure equation.
            volScalarField Ss_;
            
            //- Face-effective conductivity divided by specific weight for the
            //  `p_rgh` Darcy term.
            autoPtr<surfaceScalarField> kEffbyGammafPtr_;

            //- Flag indicating whether the chosen ddt scheme is steady-state.
            bool steadyState_;

            // Private Member Functions

            //- Disallow default bitwise copy construct
            varSatPoroFluid(const varSatPoroFluid &);

            //- Disallow default bitwise assignment
            void operator=(const varSatPoroFluid &);

            //- Allocate the face conductivity field on first use.
            void makeKEffbyGammaf();

            //- Construct the variably saturated constitutive model for the
            //  `p_rgh` formulation.
            virtual void makePoroHydraulic();

        protected:

            //- Non-const access to the effective face conductivity field.
            surfaceScalarField &kEffbyGammaf();

        public:
            //- Runtime type information
            TypeName("varSatPoroFluid");
            // Declare run-time constructor selection table

            // Constructors

            //- Construct from components
            varSatPoroFluid(
                Time &runTime,
                const word &region = dynamicFvMesh::defaultRegion,
                const bool sharedMesh = false);

            // Destructor

            virtual ~varSatPoroFluid()
            {
            }

            // Selector

            // Member Functions

            // Standard Access
            //- Expose `p_rgh` as the generic primary field.
            virtual volScalarField &pField()
            {
                return p_rgh();
            }
                
            //- Const overload of the `p_rgh` primary-field alias.
            virtual const volScalarField &pField() const
            {
                return p_rgh();
            }

            // Edit

            //- Advance one Richards-equation time step in `p_rgh` form.
            virtual bool evolve();

            //- Return the face-effective conductivity needed by q_relAcc.
            virtual tmp<surfaceScalarField> relativeAccelerationConductivity() const
            {
                return poroHydraulic().kEfff(p());
            }

            //- Estimate the next time step from the current nonlinear state.
            virtual scalar newDeltaT();

            //- Write derived diagnostic and post-processing fields.
            virtual void writeFields(const Time &runTime);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace varSatPoroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
