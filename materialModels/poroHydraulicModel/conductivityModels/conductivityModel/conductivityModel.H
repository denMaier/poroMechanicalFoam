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

Namespace
    Foam::conductivityModels

Description
    A namespace for various conductivityModel implementations.

Class
    Foam::conductivityModel

Description
    An abstract base class for conductivityModels

    Constitutive relations for soils: k(n) or k(i) or other functions

Author
    Denis Maier, BAW. All rights reserved.

SourceFiles
    conductivityModel.C
    newconductivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef conductivityModel_H
#define conductivityModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                Class conductivityModel Declaration
\*---------------------------------------------------------------------------*/

    class conductivityModel
    {


        //- Run-time-selected conductivity-model type name.
        word name_;
        //- Zone sub-dictionary backing this conductivity-model instance.
        //  Derived laws may write default entries into it so the final
        //  `_withDefaultValues` output captures the effective configuration.
        dictionary &poroHydraulicProperties_;
        //- Reference to the fluid mesh on which conductivity is evaluated.
        const fvMesh &mesh_;
        //- Registry used by models that depend on auxiliary fields.
        const objectRegistry &db_;

        // Protected data
    protected:

        //- Disallow copy construct
        conductivityModel(const conductivityModel &);

        //- Disallow default bitwise assignment
        void operator=(const conductivityModel &);

    public:
        //- Runtime type information
        TypeName("conductivityModel");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable(
            autoPtr,
            conductivityModel,
            dictionary,
            (
                const word &name,
                dictionary &poroHydraulicProperties,
                const volScalarField &pField),
            (name, poroHydraulicProperties, pField));

        // Selectors

        //- Return a reference to the selected soil model
        static autoPtr<conductivityModel> New(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField);

        // Constructors

        //- Construct from components
        conductivityModel(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField);

        // Destructor

        virtual ~conductivityModel()
        {
        }

        // Member Functions

        const objectRegistry &db() const
        {
            return db_;
        }

        const fvMesh &mesh() const
        {
            return mesh_;
        }

        //- Read transportProperties dictionary

        virtual word name()
        {
            return name_;
        }

        //- Return conductivity for an internal cell.
        virtual scalar k(const label cellI) = 0;

        //- Return conductivity for a boundary face.
        virtual scalar k
        (
            const label patchI,
            const label faceI
        ) = 0;

        virtual bool updatesK() const
        {
            return true;
        } 
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
