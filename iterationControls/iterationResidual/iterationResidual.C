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

#include "iterationResidual.H"
#include "deltaVf.H"
#include "LinearSolverRes.H"

namespace
{
    Foam::scalar parseTolerance
    (
        const Foam::word& name,
        const Foam::ITstream& stream
    )
    {
        if (!stream.size())
        {
            return -1.0;
        }

        if (stream.last().isNumber())
        {
            return stream.last().number();
        }

        if (stream.last().isWord() && stream.last().wordToken() == "show")
        {
            return -1.0;
        }

        Foam::FatalErrorInFunction
            << "last entry for residual " << name
            << " must be either a number or the word show"
            << ::Foam::abort(Foam::FatalError);

        return -1.0;
    }

    bool parseRelativeFlag(const Foam::ITstream& stream)
    {
        return stream.size() && stream.found(Foam::token(Foam::word("rel")));
    }

    Foam::word parseOperationName
    (
        const Foam::word& name,
        const Foam::ITstream& stream
    )
    {
        if (!stream.size())
        {
            return "max";
        }

        if (stream.size() < 2)
        {
            Foam::FatalErrorInFunction
                << "Residual " << name
                << " requires '<operation> <tolerance|show>'." << Foam::nl
                << "For example: " << name << " max 1e-6"
                << ::Foam::abort(Foam::FatalError);
        }

        return stream.first().wordToken();
    }
}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(iterationResidual, 0);
    defineRunTimeSelectionTable(iterationResidual, dictionary);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iterationResidual::iterationResidual
(
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField
)
:       runTime_(runTime),
        name_(name),
        tolerance_(parseTolerance(name, stream)),
        operationName_(parseOperationName(name, stream)),
        operation_(),
        residual_(GREAT),
        relative_(parseRelativeFlag(stream))
{
    relative_
         ? Info << "rel. "
         : Info << " ";
    Info  << name_
          << " Tol.: " <<  tolerance_ << endl;
    operation_.reset(residualOperation::New(operationName_));
}

Foam::iterationResidual::iterationResidual
(
    const Time& runTime,
    const word name,
    const scalar tolerance,
    const word& operationName,
    const bool relative
)
:   runTime_(runTime),
    name_(name),
    tolerance_(tolerance),
    operationName_(operationName),
    operation_(),
    residual_(GREAT),
    relative_(relative)
{
    relative_
         ? Info << "rel. "
         : Info << " ";
    Info  << name_
          << " Tol.: " <<  tolerance_ << endl;

}

Foam::autoPtr<Foam::iterationResidual> Foam::iterationResidual::New(
            const Time& runTime,
            const word name,
            const ITstream stream,
            const bool writeField)
{
    Info << "Initializing iteration control: " << nl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(name);

#else
    dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(name);

    auto* ctorPtr =
        cstrIter == dictionaryConstructorTablePtr_->end()
      ? nullptr
      : cstrIter();
#endif

    if (ctorPtr)
    {
        return autoPtr<Foam::iterationResidual>
        (
            ctorPtr(runTime, name, stream, writeField)
        );
    }

    if
    (
        name == "LinearSolver"
     || name == "linearSolver"
     || (stream.size() && stream.first().isWord() && stream.first().wordToken() == "linearSolver")
    )
    {
        return autoPtr<Foam::iterationResidual>
        (
            new Foam::LinearSolverRes
            (
                runTime,
                name,
                stream,
                writeField
            )
        );
    }

    if (Foam::deltaVf::fieldExists(runTime, name))
    {
        if (stream.size() < 2)
        {
            FatalErrorInFunction
                << "Field-based convergence entry '" << name << "' requires "
                << "'<operation> <tolerance|show>'." << nl
                << "For example: \"" << name << "\" L2 1e-3"
                << exit(FatalError);
        }

        return autoPtr<Foam::iterationResidual>
        (
            new Foam::deltaVf(runTime, name, stream, writeField)
        );
    }

    FatalErrorInFunction
        << "Unknown convergence entry '" << name << "'." << nl
        << "Expected either a registered residual type or an existing field "
        << "name." << nl
        << "Registered residual types include: linearSolver, MassBalance" << nl
        << "Trackable registered fields currently visible in the mesh registries:" << nl
        << Foam::deltaVf::candidateFieldNames(runTime) << nl << nl
        << "Only registered vol/surface scalar/vector/tensor fields can be used "
        << "for delta-based convergence checks." << nl
        << "They must also persist across the outer-iteration loop so prevIter() "
        << "state is meaningful."
        << exit(FatalError);

    return autoPtr<Foam::iterationResidual>(nullptr);

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Ostream operation  * * * * * * * * * * * * * * //

// ************************************************************************* //
