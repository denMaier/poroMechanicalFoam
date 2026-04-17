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
    Foam::storageLaws::storageCoeff

Description
    Single coefficient, that can be read either from file or from the dictionary.

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    storageCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef storageCoeff_H
#define storageCoeff_H

#include "storageLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace storageLaws
    {

        /*---------------------------------------------------------------------------*\
                           Class storageCoeff Declaration
\*---------------------------------------------------------------------------*/

        class storageCoeff
            : public storageLaw
        {
            // Private data
            dictionary& storageCoeffCoeffs_;
            //- Legacy local copy of the selected type name. The base class also
            //  stores the runtime-selected name; this member is kept only to
            //  preserve the current class layout and behavior.
            word name_;
            //- Uniform fallback storage coefficient from the dictionary.
            dimensionedScalar SsScalar_;
            //- Optional per-cell storage field read from disk when present.
            autoPtr<volScalarField> SsField_;

        public:
            //- Runtime type information
            TypeName("storageCoeff");

            // Constructors

            //- Construct from components
            storageCoeff(
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField);

            // Destructor

            ~storageCoeff()
            {
            }

            // Member Functions

            //- Read transportProperties dictionary

            //- Return the stored model name without changing the historic API.
            virtual word name() { return name_; }

            virtual scalar Ss
            (
                const scalar n,
                const scalar p,
                const label cellI
            );

            //- This model is static once configured, so callers do not need to
            //  reassemble Ss in response to evolving hydraulic state.
            virtual bool updatesSs()  const 
            {
                return false;
            }

        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace storageLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
