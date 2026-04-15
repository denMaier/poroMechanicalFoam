/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "tabulatedForcingSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "uniformDimensionedFields.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::tabulatedForcingSource::addSup
(
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    vector forceDensity
    (
        timeForce_().value
        (
            mesh_.time().timeOutputValue()
        ) / this->V()
    );

    volVectorField::Internal Su(
         IOobject(
              "forces",
              mesh_.time().constant(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE),
        mesh_,
        dimensionedVector("rhoTimesA", dimForce/dimVolume, vector::zero)
    );

    UIndirectList<vector>(Su, cells_) = forceDensity;
    eqn -=  Su;
}


// ************************************************************************* //
