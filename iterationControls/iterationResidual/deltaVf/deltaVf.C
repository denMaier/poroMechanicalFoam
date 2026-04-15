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
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam{
     defineTypeNameAndDebug(deltaVf, 0);
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
            if (found)
            {
                db_.reset(regionMesh);
                Info << "Found " << variableName_ << " in " << db_().name() << " of type " << type_ << ", with dimensions " << dimensions_ << endl;
            }
        }
    if(!db_.valid())
    {
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
                << ")"
                << nl;
        }

        FatalError
        << endl
        << "Valid iterationResiduals are : " << endl
        << "any names of output fields"
        << exit(FatalError);
    }
}

void Foam::deltaVf::lookupAndMakeScalar()
{
    if (!db_.valid())
    {
        firstLookup();
    }

    if (!vf_.valid())
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
        deltaVf_.ref() = pos(mag(vf() - vf().prevIter())-dimensionedSmall) * mag(vf() - vf().prevIter())
                        /(max(mag(vf()),dimensionedSmall));
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
            )
        )
        {
            const_cast<volScalarField&>(vf()).storePrevIter();
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
