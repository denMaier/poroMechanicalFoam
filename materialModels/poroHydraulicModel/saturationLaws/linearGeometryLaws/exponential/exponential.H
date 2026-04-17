/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::saturationLaws::exponential

Description
    Exponential SWCC.

    The exponential coefficients may be supplied as fields, so both cell and
    boundary-face evaluations use local coefficient values rather than uniform
    scalars.

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "saturationLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class exponential Declaration
\*---------------------------------------------------------------------------*/

        class exponential
            : public saturationLaw
        {
            // Private data

            dictionary& exponentialCoeffs_;
            //- Optional field-backed exponential shape parameter.
            const volScalarField exp_;
            //- Optional field-backed saturation scaling coefficient.
            const volScalarField coeffS_;
            //- Optional field-backed relative-conductivity scaling coefficient.
            const volScalarField coeffk_;
            //- Optional field-backed maximum saturation.
            const volScalarField S_0;

            // Private Member Functions

            scalar H
            (
                const scalar p,
                const scalar expScalar,
                const scalar coeffSScalar
            ) const;

        public:
            //- Runtime type information
            TypeName("exponential");

            // Constructors

            //- Construct from components
            exponential(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~exponential()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            virtual scalar pStar(const label cellI) const;

            virtual scalar C(const scalar p, const label cellI);

            virtual scalar S(const scalar p, const label cellI);

            virtual scalar S
            (
                const scalar p,
                const label patchI,
                const label faceI
            );

            virtual scalar kr(const scalar p, const label cellI);

            virtual scalar kr
            (
                const scalar p,
                const label patchI,
                const label faceI
            );

            //- No-op legacy hook; this model has no custom write-time state.
            virtual void write() const{};

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace saturationLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
