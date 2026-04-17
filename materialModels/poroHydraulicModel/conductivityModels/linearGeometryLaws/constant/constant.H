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
    Foam::conductivityModels::constant

Description
    constant k.  

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "conductivityModel.H"
#include "volFields.H"
#include "Function1.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace conductivityModels
    {

        /*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

        class constant
            : public conductivityModel
        {
            // Private data
            //- Legacy local name storage. The base class already keeps the
            //  selected runtime type; this member remains only to preserve the
            //  existing class layout.
            word name_;
            //- Uniform fallback conductivity from the dictionary.
            const dimensionedScalar kScalar_;
            //- Optional conductivity field read from disk when a spatially
            //  varying constant k is provided by the case.
            autoPtr<volScalarField> k_;         

        public:
            //- Runtime type information
            TypeName("constant");

            // Constructors

            //- Construct from components
            constant(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~constant()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            //- Return the historic user-facing model name.
            virtual word name() { return "constant"; }

            virtual scalar k(const label cellI);

            virtual scalar k
            (
                const label patchI,
                const label faceI
            );
            
            //- This model is static once configured, so callers do not need to
            //  reassemble k in response to evolving hydraulic state.
            virtual bool updatesK() const
            {
                return false;
            }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        };
    }; // End namespace conductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
