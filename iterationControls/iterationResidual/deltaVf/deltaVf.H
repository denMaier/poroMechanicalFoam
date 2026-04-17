/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::deltaVf

Description
    Using the change in any field as the convergence critera.

Author
    Denis Maier, BAW. All rights reserved.

SeeAlso
    fvSolution

SourceFiles
    deltaVf.C

\*---------------------------------------------------------------------------*/

#ifndef delta_H
#define delta_H

#include "iterationResidual.H"
#include "scalar.H"
#include "surfaceFields.H"
#include "wordList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class deltaVf Declaration
\*---------------------------------------------------------------------------*/
class deltaVf
:
    public iterationResidual
{
    // Private Member Functions
        Switch writeField_;
        word variableName_;
        scalar sqrt32_;
        refPtr<objectRegistry> db_;
        refPtr<fvMesh> mesh_;
        refPtr<dimensionSet> dimensions_;
        word type_;
        autoPtr<volScalarField> vf_;
        autoPtr<volScalarField> deltaVf_;
        autoPtr<surfaceScalarField> sf_;
        autoPtr<surfaceScalarField> deltaSf_;
        bool prevIterStored_;

        bool isSurfaceType() const
        {
            return type_ == "surfaceScalar" || type_ == "surfaceVector";
        }

        volScalarField& vf()
        {
            if (!db_.valid())
            {
                firstLookup();
            }
            if(!vf_.valid())
            {
                makeVfScalar();
            }
            return vf_.ref();
        }

        surfaceScalarField& sf()
        {
            if (!db_.valid())
            {
                firstLookup();
            }
            if(!sf_.valid())
            {
                makeSfScalar();
            }
            return sf_.ref();
        }

        void firstLookup();
        void lookupAndMakeScalar();
        void makeVfScalar();
        void makeDeltaVf();
        void makeSfScalar();
        void makeDeltaSf();

        //- No copy construct
        deltaVf(const deltaVf&) = delete;

public:
    //- runTime type information
        TypeName("delta");

        static bool fieldExists
        (
            const Time& runTime,
            const word& variableName
        );

        static wordList candidateFieldNames(const Time& runTime);

    // Constructors
        deltaVf
        (
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField
        );


    //- Destructor
    ~deltaVf();


    // Member Functions

        //- Execute specified function names
        virtual scalar calcResidual();

        virtual word name() const
        {
            return relative_
                ? "rel. " + iterationResidual::name()
                : iterationResidual::name();
        }

        virtual void reset();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
