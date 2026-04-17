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
    varSatPoroFluidHead

Description
    Richards-equation solver for variably saturated flow using hydraulic head
    `pHead` as the primary unknown.

    Compared with `varSatPoroFluid`, this formulation keeps gravity explicit in
    the governing equation. The solved field is hydraulic head in units of water
    column, and total pore pressure is reconstructed afterwards as
    `p = gamma * pHead`.

    This class is the head-based concrete child of
    `variablySaturatedPoroFluid`. It reuses the same constitutive and
    linearization infrastructure as the `p_rgh` version, but assembles the
    equation in a form that is natural for head-based boundary conditions and
    post-processing.

Author
    Denis Maier, BAW.  All rights reserved.

SourceFiles
    varSatPoroFluidHead.C

\*---------------------------------------------------------------------------*/

#ifndef varSatPoroFluidHead_H
#define varSatPoroFluidHead_H

#include "variablySaturatedPoroFluid.H"
#include "richardsLinearization.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
        namespace poroFluidModels
        {
        /*---------------------------------------------------------------------------*\
                        Class varSatPoroFluidHead Declaration
        \*---------------------------------------------------------------------------*/

        /**
         * Variably saturated Richards solver in hydraulic-head form.
         *
         * The class exposes `pHead_` as the primary field, reconstructs total
         * pressure after convergence, and carries a face-effective conductivity
         * field appropriate for the head-based Darcy formulation.
         */
        class varSatPoroFluidHead
            : public variablySaturatedPoroFluid
        {
            // Private data
            //- Primary unknown: hydraulic head of the pore fluid.
            volScalarField pHead_;

            //- Compressibility / storage contribution used in the pressure equation.
            volScalarField Ss_;

            //- Face-effective hydraulic conductivity for the head formulation.
            autoPtr<surfaceScalarField> kEfffPtr_;

            //- Flag indicating whether the chosen ddt scheme is steady-state.
            bool steadyState_;

            // Private Member Functions

            //- Disallow default bitwise copy construct
            varSatPoroFluidHead(const varSatPoroFluidHead &);

            //- Disallow default bitwise assignment
            void operator=(const varSatPoroFluidHead &);

            //- Allocate the face conductivity field on first use.
            void makeKEfff();

            //- Construct the variably saturated constitutive model for the
            //  `pHead` formulation.
            virtual void makePoroHydraulic();

        protected:
            
            //- Non-const access to the effective face conductivity field.
            surfaceScalarField &kEfff();

        public:
            //- Runtime type information
            TypeName("varSatPoroFluidHead");

            // Declare run-time constructor selection table

            // Constructors

            //- Construct from components
            varSatPoroFluidHead(
                    Time &runTime,
                    const word &region = dynamicFvMesh::defaultRegion,
                    const bool sharedMesh = false);

            // Destructor

            virtual ~varSatPoroFluidHead()
            {
            }

            // Selector

            // Member Functions

            // Access

            //- Expose hydraulic head as the generic primary field.
            virtual volScalarField &pField()
            {
                return pHead_;
            }
                
            //- Const overload of the hydraulic-head primary-field alias.
            virtual const volScalarField &pField() const
            {
                return pHead_;
            }

            // Edit

            //- Advance one Richards-equation time step in hydraulic-head form.
            virtual bool evolve();

            //- Estimate the next time step from the current nonlinear state.
            virtual scalar newDeltaT();

            //- Write derived diagnostic and post-processing fields.
            virtual void writeFields(const Time &runTime);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        } // End namespace poroFluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
