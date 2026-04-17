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

#include "deltaVf.H"
#include "DynamicList.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam{
     defineTypeNameAndDebug(deltaVf, 0);
}

namespace
{
    bool isTrackableField(const Foam::objectRegistry& registry, const Foam::word& name)
    {
        return
        (
            registry.foundObject<Foam::volScalarField>(name)
         || registry.foundObject<Foam::volVectorField>(name)
         || registry.foundObject<Foam::volSymmTensorField>(name)
         || registry.foundObject<Foam::volTensorField>(name)
         || registry.foundObject<Foam::surfaceScalarField>(name)
         || registry.foundObject<Foam::surfaceVectorField>(name)
        );
    }
}

bool Foam::deltaVf::fieldExists
(
    const Time& runTime,
    const word& variableName
)
{
    HashTable<const fvMesh*> meshes = runTime.lookupClass<const fvMesh>();

    forAllConstIters(meshes, meshIter)
    {
        const fvMesh& regionMesh = *(meshIter.val());

        if
        (
            regionMesh.objectRegistry::foundObject<volScalarField>(variableName)
         || regionMesh.objectRegistry::foundObject<volVectorField>(variableName)
         || regionMesh.objectRegistry::foundObject<volSymmTensorField>(variableName)
         || regionMesh.objectRegistry::foundObject<volTensorField>(variableName)
         || regionMesh.objectRegistry::foundObject<surfaceScalarField>(variableName)
         || regionMesh.objectRegistry::foundObject<surfaceVectorField>(variableName)
        )
        {
            return true;
        }
    }

    return false;
}

Foam::wordList Foam::deltaVf::candidateFieldNames(const Time& runTime)
{
    DynamicList<word> uniqueNames;
    HashTable<const fvMesh*> meshes = runTime.lookupClass<const fvMesh>();

    forAllConstIters(meshes, meshIter)
    {
        const fvMesh& regionMesh = *(meshIter.val());
        const wordList registryNames(regionMesh.sortedNames());

        forAll(registryNames, iName)
        {
            const word& name = registryNames[iName];

            if (!isTrackableField(regionMesh, name))
            {
                continue;
            }

            bool alreadyAdded = false;
            forAll(uniqueNames, iUnique)
            {
                if (uniqueNames[iUnique] == name)
                {
                    alreadyAdded = true;
                    break;
                }
            }

            if (!alreadyAdded)
            {
                uniqueNames.append(name);
            }
        }
    }

    return wordList(uniqueNames);
}

void Foam::deltaVf::firstLookup()
{
    HashTable<const fvMesh*> meshes = runTime().lookupClass<const fvMesh>();
    forAllConstIters(meshes, meshIter)
    {
        bool found = false;
        const fvMesh& regionMesh = *(meshIter.val());
        if(regionMesh.objectRegistry::foundObject<volScalarField>(variableName_))
        {
            type_ = "scalar";
            found = true;
            const volScalarField& tvf = regionMesh.objectRegistry::lookupObject<volScalarField>(variableName_);
            dimensions_.reset(tvf().dimensions());
            mesh_.reset(tvf().mesh());
        }
        else if(regionMesh.objectRegistry::foundObject<volVectorField>(variableName_))
        {
            found = true;
            type_ = "vector";
            const volVectorField& tvf = regionMesh.objectRegistry::lookupObject<volVectorField>(variableName_);
            dimensions_.reset(tvf().dimensions());
            mesh_.reset(tvf().mesh());
        }
        else if(regionMesh.objectRegistry::foundObject<volSymmTensorField>(variableName_))
        {
            found = true;
            type_ = "symmTensor";
            const volSymmTensorField& tvf = regionMesh.objectRegistry::lookupObject<volSymmTensorField>(variableName_);
            dimensions_.reset(tvf().dimensions());
            mesh_.reset(tvf().mesh());
        }
        else if(regionMesh.objectRegistry::foundObject<volTensorField>(variableName_))
        {
            found = true;
            type_ = "tensor";
            const volTensorField& tvf = regionMesh.objectRegistry::lookupObject<volTensorField>(variableName_);
            dimensions_.reset(tvf().dimensions());
            mesh_.reset(tvf().mesh());
        }
        else if(regionMesh.objectRegistry::foundObject<surfaceScalarField>(variableName_))
        {
            found = true;
            type_ = "surfaceScalar";
            const surfaceScalarField& tvf =
                regionMesh.objectRegistry::lookupObject<surfaceScalarField>(variableName_);
            dimensions_.reset(tvf().dimensions());
            mesh_.reset(tvf().mesh());
        }
        else if(regionMesh.objectRegistry::foundObject<surfaceVectorField>(variableName_))
        {
            found = true;
            type_ = "surfaceVector";
            const surfaceVectorField& tvf =
                regionMesh.objectRegistry::lookupObject<surfaceVectorField>(variableName_);
            dimensions_.reset(tvf().dimensions());
            mesh_.reset(tvf().mesh());
        }
        if (found)
        {
            db_.reset(regionMesh);
            Info << "Found " << variableName_ << " in " << db_().name() << " of type " << type_ << ", with dimensions " << dimensions_ << endl;
            break;
        }
    }
    if(!db_.valid())
    {
        const wordList candidates(candidateFieldNames(runTime()));

        FatalErrorIn(
        "iterationResidual::New(Time&, "
        "word&, "
        "ITstream&) ")
        << "Unknown iterationResidual type "
        << variableName_ << nl
        << "Searched fvMesh object registries:" << nl;

        forAllConstIters(meshes, meshIter)
        {
            const fvMesh& regionMesh = *(meshIter.val());

            FatalError
                << "  - " << regionMesh.name()
                << " (scalar="
                << regionMesh.objectRegistry::foundObject<volScalarField>(variableName_)
                << ", vector="
                << regionMesh.objectRegistry::foundObject<volVectorField>(variableName_)
                << ", symmTensor="
                << regionMesh.objectRegistry::foundObject<volSymmTensorField>(variableName_)
                << ", tensor="
                << regionMesh.objectRegistry::foundObject<volTensorField>(variableName_)
                << ", surfaceScalar="
                << regionMesh.objectRegistry::foundObject<surfaceScalarField>(variableName_)
                << ", surfaceVector="
                << regionMesh.objectRegistry::foundObject<surfaceVectorField>(variableName_)
                << ")"
                << nl;
        }

        FatalError
        << endl
        << "Trackable registered fields currently visible in the mesh registries:" << nl
        << candidates << nl << nl
        << "Only registered vol/surface scalar/vector/tensor fields can be used "
        << "for delta-based convergence checks." << nl
        << "They must also persist across the outer-iteration loop so prevIter() "
        << "state is meaningful."
        << exit(FatalError);
    }
}

void Foam::deltaVf::lookupAndMakeScalar()
{
    if (!db_.valid())
    {
        firstLookup();
    }

    if (isSurfaceType())
    {
        if (!sf_.valid())
        {
            makeSfScalar();
        }
    }
    else if (!vf_.valid())
    {
        makeVfScalar();
    }

    if(type_ == "scalar")
    {
        const volScalarField& src = db_().objectRegistry::lookupObject<volScalarField>(variableName_);

        if (vf().dimensions() != src.dimensions())
        {
            WarningInFunction
                << "Reinitializing residual field tracker for '" << variableName_
                << "' on mesh '" << db_().name() << "' due to dimension change "
                << vf().dimensions() << " -> " << src.dimensions() << endl;

            dimensions_.reset(src.dimensions());
            vf_.clear();
            makeVfScalar();
            prevIterStored_ = false;
        }

        vf() = src;
        dimensions_.reset(vf().dimensions());
    }
    else if (type_ == "vector")
    {
        const volVectorField& src = db_().objectRegistry::lookupObject<volVectorField>(variableName_);

        if (vf().dimensions() != src.dimensions())
        {
            WarningInFunction
                << "Reinitializing residual field tracker for '" << variableName_
                << "' on mesh '" << db_().name() << "' due to dimension change "
                << vf().dimensions() << " -> " << src.dimensions() << endl;

            dimensions_.reset(src.dimensions());
            vf_.clear();
            makeVfScalar();
            prevIterStored_ = false;
        }

        vf() = mag(src);
        dimensions_.reset(vf().dimensions());
    }
    else if (type_ == "symmTensor")
    {
        const volSymmTensorField& src = db_().objectRegistry::lookupObject<volSymmTensorField>(variableName_);

        if (vf().dimensions() != src.dimensions())
        {
            WarningInFunction
                << "Reinitializing residual field tracker for '" << variableName_
                << "' on mesh '" << db_().name() << "' due to dimension change "
                << vf().dimensions() << " -> " << src.dimensions() << endl;

            dimensions_.reset(src.dimensions());
            vf_.clear();
            makeVfScalar();
            prevIterStored_ = false;
        }

        vf() = sqrt32_*sqrt(magSqr(dev(src)));
        dimensions_.reset(vf().dimensions());
    }
    else if (type_ == "tensor")
    {
        const volTensorField& src = db_().objectRegistry::lookupObject<volTensorField>(variableName_);

        if (vf().dimensions() != src.dimensions())
        {
            WarningInFunction
                << "Reinitializing residual field tracker for '" << variableName_
                << "' on mesh '" << db_().name() << "' due to dimension change "
                << vf().dimensions() << " -> " << src.dimensions() << endl;

            dimensions_.reset(src.dimensions());
            vf_.clear();
            makeVfScalar();
            prevIterStored_ = false;
        }

        vf() = sqrt32_*sqrt(magSqr(dev(src)));
        dimensions_.reset(vf().dimensions());
    }
    else if (type_ == "surfaceScalar")
    {
        const surfaceScalarField& src =
            db_().objectRegistry::lookupObject<surfaceScalarField>(variableName_);

        if (sf().dimensions() != src.dimensions())
        {
            WarningInFunction
                << "Reinitializing residual field tracker for '" << variableName_
                << "' on mesh '" << db_().name() << "' due to dimension change "
                << sf().dimensions() << " -> " << src.dimensions() << endl;

            dimensions_.reset(src.dimensions());
            sf_.clear();
            makeSfScalar();
            prevIterStored_ = false;
        }

        sf() = src;
        dimensions_.reset(sf().dimensions());
    }
    else if (type_ == "surfaceVector")
    {
        const surfaceVectorField& src =
            db_().objectRegistry::lookupObject<surfaceVectorField>(variableName_);

        if (sf().dimensions() != src.dimensions())
        {
            WarningInFunction
                << "Reinitializing residual field tracker for '" << variableName_
                << "' on mesh '" << db_().name() << "' due to dimension change "
                << sf().dimensions() << " -> " << src.dimensions() << endl;

            dimensions_.reset(src.dimensions());
            sf_.clear();
            makeSfScalar();
            prevIterStored_ = false;
        }

        sf() = mag(src);
        dimensions_.reset(sf().dimensions());
    }
}


void Foam::deltaVf::makeVfScalar()
{
    word residualFieldName;
    if (type_ == "scalar")
    {
        residualFieldName = name();
    }
    else if (type_ == "vector")
    {
        residualFieldName = "mag("+name()+")";
    }
    else
    {
        residualFieldName = name()+"Eq";
    }
    vf_.reset(
        new volScalarField(
        IOobject
        (
         residualFieldName,
         runTime().timeName(),
         db_(),
        IOobject::NO_READ,
        IOobject::NO_WRITE),
        mesh_(),
        dimensionedScalar("",dimensions_(),0.0))
    );
}

void Foam::deltaVf::makeSfScalar()
{
    word residualFieldName;
    if (type_ == "surfaceScalar")
    {
        residualFieldName = name();
    }
    else if (type_ == "surfaceVector")
    {
        residualFieldName = "mag("+name()+")";
    }
    else
    {
        residualFieldName = name()+"Eq";
    }

    sf_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                residualFieldName,
                runTime().timeName(),
                db_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_(),
            dimensionedScalar("", dimensions_(), 0.0)
        )
    );
}

void Foam::deltaVf::makeDeltaVf()
{
    deltaVf_.reset(
        new volScalarField(
        IOobject
        (
         name()+"Residual",
         runTime().timeName(),
         vf_().db(),
        IOobject::NO_READ,
        writeField_
        ?IOobject::AUTO_WRITE
        :IOobject::NO_WRITE),
        vf_().mesh(),
        dimensionedScalar("",relative_?dimless:dimensions_(),0.0))
    );
}

void Foam::deltaVf::makeDeltaSf()
{
    deltaSf_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                name()+"Residual",
                runTime().timeName(),
                sf_().db(),
                IOobject::NO_READ,
                writeField_
              ? IOobject::AUTO_WRITE
              : IOobject::NO_WRITE
            ),
            sf_().mesh(),
            dimensionedScalar("", relative_ ? dimless : dimensions_(), 0.0)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deltaVf::deltaVf
(
    const Time& runTime,
    const word name,
    const ITstream stream,
    const bool writeField
)
:
    iterationResidual(runTime, "delta("+name+")", stream, writeField),
    writeField_(writeField),
    variableName_(name),
    sqrt32_(sqrt(3.0/2.0)),
    db_(),
    mesh_(),
    dimensions_(),
    type_("non"),
    vf_(),
    deltaVf_(),
    sf_(),
    deltaSf_(),
    prevIterStored_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::deltaVf::~deltaVf()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::deltaVf::calcResidual()
{

    lookupAndMakeScalar();

    if (isSurfaceType())
    {
        if(!deltaSf_.valid())
        {
            makeDeltaSf();
        }
        else
        {
            const dimensionSet expectedDims(relative_ ? dimless : sf().dimensions());

            if (deltaSf_().dimensions() != expectedDims)
            {
                deltaSf_.clear();
                makeDeltaSf();
                prevIterStored_ = false;
            }
        }

        if(!prevIterStored_)
        {
            const_cast<surfaceScalarField&>(sf()).storePrevIter();
            prevIterStored_ = true;
        }

        if(relative_)
        {
            dimensionedScalar dimensionedSmall("",dimensions_(),SMALL);
            deltaSf_.ref() =
                pos(mag(sf() - sf().prevIter()) - dimensionedSmall)
              * mag(sf() - sf().prevIter())
              / max(mag(sf()), dimensionedSmall);
        }
        else
        {
            deltaSf_.ref() = mag(sf() - sf().prevIter());
        }

        residual_ = operation(deltaSf_().primitiveField());
        const_cast<surfaceScalarField&>(sf()).storePrevIter();
        prevIterStored_ = true;

        if(!writeField_)
        {
            deltaSf_.clear();
        }
    }
    else
    {
        if(!deltaVf_.valid())
        {
            makeDeltaVf();
        }
        else
        {
            const dimensionSet expectedDims(relative_ ? dimless : vf().dimensions());

            if (deltaVf_().dimensions() != expectedDims)
            {
                deltaVf_.clear();
                makeDeltaVf();
                prevIterStored_ = false;
            }
        }

        if(!prevIterStored_)
        {
            const_cast<volScalarField&>(vf()).storePrevIter();
            prevIterStored_ = true;
        }

        if(relative_)
        {
            dimensionedScalar dimensionedSmall("",dimensions_(),SMALL);
            deltaVf_.ref() =
                pos(mag(vf() - vf().prevIter())-dimensionedSmall)
              * mag(vf() - vf().prevIter())
              / (max(mag(vf()),dimensionedSmall));
        }
        else
        {
            deltaVf_.ref() =  mag(vf() - vf().prevIter());
        }

        residual_ = operation(deltaVf_().primitiveField());
        const_cast<volScalarField&>(vf()).storePrevIter();
        prevIterStored_ = true;

        if(!writeField_)
        {
            deltaVf_.clear();
        }
    }

    return residual_;
}

void Foam::deltaVf::reset()
    {
        if
        (
            db_.valid()
         &&
            (
                (type_ == "scalar" && db_().objectRegistry::foundObject<volScalarField>(variableName_))
             || (type_ == "vector" && db_().objectRegistry::foundObject<volVectorField>(variableName_))
             || (type_ == "symmTensor" && db_().objectRegistry::foundObject<volSymmTensorField>(variableName_))
             || (type_ == "tensor" && db_().objectRegistry::foundObject<volTensorField>(variableName_))
             || (type_ == "surfaceScalar" && db_().objectRegistry::foundObject<surfaceScalarField>(variableName_))
             || (type_ == "surfaceVector" && db_().objectRegistry::foundObject<surfaceVectorField>(variableName_))
            )
        )
        {
            if (isSurfaceType())
            {
                const_cast<surfaceScalarField&>(sf()).storePrevIter();
            }
            else
            {
                const_cast<volScalarField&>(vf()).storePrevIter();
            }
            prevIterStored_ = true;
        }
        else
        {
            prevIterStored_ = false;
        }

        residual_ = GREAT;
    }
// * * * * * * * * * * * * * * Ostream operation  * * * * * * * * * * * * * * //

// ************************************************************************* //
