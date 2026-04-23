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

\*---------------------------------------------------------------------------*/

#include "LinearSolverRes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace
{
    Foam::scalar parseLinearSolverTolerance
    (
        const Foam::word& fieldName,
        const Foam::ITstream& stream
    )
    {
        if (!stream.size())
        {
            FatalErrorInFunction
                << "linearSolver convergence entry for field '" << fieldName
                << "' requires at least the keyword 'linearSolver'."
                << Foam::exit(Foam::FatalError);
        }

        if (!stream.first().isWord() || stream.first().wordToken() != "linearSolver")
        {
            FatalErrorInFunction
                << "linearSolver convergence entry for field '" << fieldName
                << "' must start with the keyword 'linearSolver'."
                << Foam::exit(Foam::FatalError);
        }

        if (stream.size() == 1)
        {
            return -1.0;
        }

        if (stream.size() != 2)
        {
            FatalErrorInFunction
                << "linearSolver convergence entry for field '" << fieldName
                << "' must be 'linearSolver <tolerance|show>'."
                << Foam::exit(Foam::FatalError);
        }

        if (stream.last().isNumber())
        {
            return stream.last().number();
        }

        if (stream.last().isWord() && stream.last().wordToken() == "show")
        {
            return -1.0;
        }

        FatalErrorInFunction
            << "linearSolver convergence entry for field '" << fieldName
            << "' must end with a number or the word show."
            << Foam::exit(Foam::FatalError);

        return -1.0;
    }

    template<class GeoField, class ValueType>
    Foam::scalar initialResidualForField
    (
        const Foam::fvMesh& mesh,
        const Foam::word& fieldName
    )
    {
        if (!mesh.foundObject<GeoField>(fieldName))
        {
            return -1.0;
        }

        const Foam::dictionary& solverDict = mesh.data().solverPerformanceDict();

        if (!solverDict.found(fieldName))
        {
            FatalErrorInFunction
                << "Field '" << fieldName << "' exists on mesh '"
                << mesh.name() << "' but has no solverPerformance entry in the "
                << "current outer iteration." << Foam::nl
                << "Make sure the field is actually solved before convergence "
                << "is checked."
                << Foam::exit(Foam::FatalError);
        }

        const Foam::List<Foam::SolverPerformance<ValueType>> sp
        (
            solverDict.lookup(fieldName)
        );

        if (!sp.size())
        {
            FatalErrorInFunction
                << "solverPerformance entry for field '" << fieldName
                << "' on mesh '" << mesh.name() << "' is empty."
                << Foam::exit(Foam::FatalError);
        }

        return Foam::cmptMax(sp.first().initialResidual());
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam{
     defineTypeNameAndDebug(LinearSolverRes, 0);

    addToRunTimeSelectionTable(
            iterationResidual,
            LinearSolverRes,
            dictionary);
}

void Foam::LinearSolverRes::findMesh()
{
    HashTable<const fvMesh*> meshes = runTime().lookupClass<const fvMesh>();

    forAllConstIters(meshes, meshIter)
    {
        const fvMesh& regionMesh = *(meshIter.val());

        if
        (
            regionMesh.foundObject<volScalarField>(fieldName_)
         || regionMesh.foundObject<volVectorField>(fieldName_)
         || regionMesh.foundObject<volSymmTensorField>(fieldName_)
         || regionMesh.foundObject<volTensorField>(fieldName_)
        )
        {
            meshPtr_ = &regionMesh;
            return;
        }
    }

    FatalErrorInFunction
        << "No supported volume field named '" << fieldName_
        << "' was found for linearSolver convergence." << nl
        << "linearSolver convergence in OpenFOAM v2412 is field-based and uses "
        << "solverPerformanceDict entries keyed by solved volume fields."
        << exit(FatalError);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LinearSolverRes::LinearSolverRes
(
    const Time& runTime,
    const word name,
    const ITstream stream,
    const bool writeField
)
:
    iterationResidual
    (
        runTime,
        name,
        parseLinearSolverTolerance(name, stream),
        "linearSolver",
        false
    ),
    fieldName_(name),
    meshPtr_(nullptr)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LinearSolverRes::~LinearSolverRes()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::LinearSolverRes::calcResidual()
{
    if (!meshPtr_)
    {
        findMesh();
    }

    if (meshPtr_->foundObject<volScalarField>(fieldName_))
    {
        residual_ = initialResidualForField<volScalarField, scalar>
        (
            *meshPtr_,
            fieldName_
        );
    }
    else if (meshPtr_->foundObject<volVectorField>(fieldName_))
    {
        residual_ = initialResidualForField<volVectorField, vector>
        (
            *meshPtr_,
            fieldName_
        );
    }
    else if (meshPtr_->foundObject<volSymmTensorField>(fieldName_))
    {
        residual_ = initialResidualForField<volSymmTensorField, symmTensor>
        (
            *meshPtr_,
            fieldName_
        );
    }
    else if (meshPtr_->foundObject<volTensorField>(fieldName_))
    {
        residual_ = initialResidualForField<volTensorField, tensor>
        (
            *meshPtr_,
            fieldName_
        );
    }
    else
    {
        FatalErrorInFunction
            << "Field '" << fieldName_ << "' disappeared from mesh '"
            << meshPtr_->name() << "' before the linearSolver residual could "
            << "be evaluated."
            << exit(FatalError);
    }

    return residual_;
}


// * * * * * * * * * * * * * * Ostream operation  * * * * * * * * * * * * * * //

// ************************************************************************* //
