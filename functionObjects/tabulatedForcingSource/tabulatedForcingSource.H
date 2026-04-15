/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::fv::tabulatedForcingSource

Group
    grpFvOptionsSources

Description
    Solid-body 6-DoF acceleration source

Usage
    Example usage:
    \verbatim
    SBM
    {
        type            tabulatedForcingSource;
        active          yes;

        timeDataFileName "constant/acceleration-terms.dat";
    }
    \endverbatim

SourceFiles
    tabulatedForcingSource.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedAccelerationSource_H
#define tabulatedAccelerationSource_H

#include "cellSetOption.H"
//#include "tabulated6DoFAcceleration.H"
#include "Function1.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class tabulatedForcingSource Declaration
\*---------------------------------------------------------------------------*/

class tabulatedForcingSource
:
    public cellSetOption
{
protected:

    // Protected data

        autoPtr<Function1<vector>> timeForce_;
        //- Velocity field name, default = U
        word UName_;

private:

    // Private Member Functions

        //- No copy construct
        tabulatedForcingSource
        (
            const tabulatedForcingSource&
        ) = delete;

        //- No copy assignment
        void operator=(const tabulatedForcingSource&) = delete;


        //- Source term to momentum equation
        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

public:

    //- Runtime type information
    TypeName("tabulatedForcingSource");


    // Constructors

        //- Construct from components
        tabulatedForcingSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~tabulatedForcingSource() = default;


    // Member Functions

        //- Source term to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tabulatedForcingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
