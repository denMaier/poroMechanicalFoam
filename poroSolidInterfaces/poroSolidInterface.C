/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "poroSolidInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvOption.H"
#include "fvMatrices.H"
#include "fvc.H"
#include "poroMechanicalLaw2.H"
#include "poroCouplingTerms.H"
#include "poroCouplingRegistry.H"

namespace Foam
{
        defineTypeNameAndDebug(poroSolidInterface, 0);
        defineRunTimeSelectionTable(poroSolidInterface, dictionary);
        addNamedToRunTimeSelectionTable(physicsModel, poroSolidInterface, physicsModel, poroSolid);
}

namespace
{
    bool looksLikeLegacyPoroMechanicalLaw(const Foam::mechanicalLaw& law)
    {
        return law.type().find("poroMechanicalLaw") != std::string::npos;
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::poroSolidInterface::makeIterCtrl()
    {
        iterCtrl_.reset(new iterationControl(runTime_,interactionProperties_,"Pressure-Displacement"));
    }

void Foam::poroSolidInterface::makePoroFluidCouplingSource()
    {
        dictionary poroSolidSourceDict = subOrEmptyDict("DtoPorePressureCouplingRegions"); //  from v1912 onwards: subDictOrAdd("kinematicPressureSource");
        poroSolidSourceDict.add("type","poroSolidToFluidCouplingSource");
        poroSolidSourceDict.lookupOrAddDefault<word>("selectionMode","all");//  from v1912 onwards: getOrAdd<word>("selectionMode","all");
        poroFluidRef().fvOptions().PtrList<fv::option>::append(fv::option::New("poroSolidToFluidCouplingSource",poroSolidSourceDict,poroFluidMesh()));
    }

Foam::tmp<Foam::volScalarField> Foam::poroSolidInterface::nDot
(
    const volScalarField& b,
    const volVectorField& U
) const
{
    return poroCouplingTerms::nDot(b, U);
}

Foam::tmp<Foam::volScalarField> Foam::poroSolidInterface::fixedStressStabil
(
    const volScalarField& b,
    const volScalarField& impK
) const
{
    return poroCouplingTerms::fixedStressStabil(b, impK);
}

void Foam::poroSolidInterface::updateCouplingTerms
(
    const volScalarField& couplingCoeff,
    const volScalarField& impK,
    const volVectorField& U,
    autoPtr<volScalarField>& nDotField,
    autoPtr<volScalarField>& fixedStressStabilField
) const
{
    poroCouplingTerms::updateCouplingFields
    (
        couplingCoeff,
        b(),
        impK,
        U,
        nDotField,
        fixedStressStabilField
    );
}

void Foam::poroSolidInterface::mapPressuresToSolidMesh
(
    autoPtr<volScalarField>& pSolidMeshField,
    autoPtr<volScalarField>& pRghSolidMeshField
)
{
    if(sharedMesh())
    {
        FatalErrorInFunction
            << "Pressure mapping requested although meshes are shared"
            << exit(FatalError);
    }

    if(!pSolidMeshField.valid() || !pRghSolidMeshField.valid())
    {
        FatalErrorInFunction
            << "Pressure fields on solid mesh are not initialized for non-shared mesh coupling"
            << exit(FatalError);
    }

    pRghSolidMeshField.ref() = solidToPoroFluid().mapSrcToTgt(poroFluid().p_rgh())();
    pSolidMeshField.ref() = solidToPoroFluid().mapSrcToTgt(poroFluid().p());
}

void Foam::poroSolidInterface::ensureSharedSolidFieldRegistered
(
    volScalarField& field,
    const word& ownerRegistryName
)
{
    if(!sharedMesh())
    {
        FatalErrorInFunction
            << "Shared solid-field registration requested although meshes are not shared"
            << exit(FatalError);
    }

    objectRegistry& solidRegistry =
        const_cast<objectRegistry&>(static_cast<const objectRegistry&>(solidMesh()));

    poroCouplingRegistry::ensureVolScalarFieldRegistered
    (
        solidRegistry,
        field,
        ownerRegistryName,
        sharedSolidRegistryFieldNames_
    );
}

void Foam::poroSolidInterface::updatePorosityFromSolidDisplacement()
{
    if(sharedMesh())
    {
        poroFluidRef().update_porosity(fvc::div(solid().D()), false);
    }
    else
    {
        tmp<volVectorField> DFluidMesh = solidToPoroFluid().mapTgtToSrc(solid().D());
        poroFluidRef().update_porosity(fvc::div(DFluidMesh), false);
        DFluidMesh.clear();
    }

    afterPorosityUpdate();
}

bool Foam::poroSolidInterface::checkMechanicalLawUpdateBiotCoeff(const mechanicalLaw& law, const label lawI, PtrList<volScalarField> &bs)
    {
        bool isCoupled = false;
        //Info << law.type() << " " << law.name()  << endl;

        const poroMechanicalLaw2* pmlPtr = dynamic_cast<const poroMechanicalLaw2*>(&law);

        if (pmlPtr)
        {
            bs.set
            (
                lawI,
                new volScalarField(pmlPtr->biotCoeff())
            );

            isCoupled = true;
        }
        else if (looksLikeLegacyPoroMechanicalLaw(law))
        {
            Warning() << " mechanicalLaw '" << law.type() << "' selected,"
                      << " Biot-Willis coefficient will only be regarded in effective stress formulation!" << nl
                      << " Please consider using a law with explicit Biot-coefficient support"
                      << " such as poroMechanicalLaw2 for consistent treatment!"
                      << endl;
            isCoupled = true;
            bs.set
            (
                lawI,
                new volScalarField(
                    IOobject(
                    "subMeshBiotCoeff",
                    runTime().timeName(),
                    law.rho()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                    ),
                    law.rho()().mesh(),
                    dimensionedScalar("",dimless,1.0)
                    )
            );
        }
        else
        {
            bs.set
            (
                lawI,
                new volScalarField(
                    IOobject(
                    "subMeshBiotCoeff",
                    runTime().timeName(),
                    law.rho()().db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                    ),
                    law.rho()().mesh(),
                    dimensionedScalar("",dimless,1.0)
                    )
            );
            Warning() << "mechanicalLaw " << law.type() << " for material " << law.name()
                            << " does not expose explicit Biot-coefficient support." << nl
                            << "This means there is no coupling from poroFluid to solid, "
                            << "and Biot Coefficient will be set to 1.0 in this region!" << nl
                            << "To account for other values, the storage term in poroHydraulicProperties can be modified."
                            << endl;
        }
        return isCoupled;
    }

void Foam::poroSolidInterface::makeBiotCoeff()
{
    tmp<volScalarField> tb(
        new volScalarField(
            IOobject(
                "biotCoeff",
                runTime().constant(),
                solidMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            solidMesh(),
            dimensionedScalar("",dimless,1.0)
        )
    );
    const PtrList<mechanicalLaw>* lawsPtr = dynamic_cast<const PtrList<mechanicalLaw>*>(&solid().mechanical());
    if (!lawsPtr)
    {
        FatalErrorInFunction
            << "Could not access mechanicalLaw list from solid mechanical model"
            << exit(FatalError);
    }
    const PtrList<mechanicalLaw>& laws = *lawsPtr;
    // Accumulated subMesh fields and then map to the base mesh
    PtrList<volScalarField> bs(laws.size());

    bool foundPMLaw = false;

    forAll(laws, lawI)
    {
        foundPMLaw = checkMechanicalLawUpdateBiotCoeff(laws[lawI],lawI,bs) || foundPMLaw;
    }
    if(!foundPMLaw)
    {
        WarningInFunction() << "No 'poroMechanicalLaw2' entry was found in mechanicalProperties." << nl
                            << "This means there is no explicit poroFluid-to-solid coupling law, "
                            << "and the Biot coefficient will default to 1.0 everywhere." << nl
                            << "If you need other values, either use poroMechanicalLaw2 or "
                            << "adjust the storage term in poroHydraulicProperties accordingly."
                            << endl;
    }

    // Map subMesh fields to the base mesh
    if (laws.size()==1)
    {
        tb.ref() = bs[0];
    }
    else
    {
        solid().mechanical().solSubMeshes().mapSubMeshVolFields<scalar>(bs, tb.ref());
    }
    // Clear subMesh fields
    bs.clear();
    if (sharedMesh_)
    {
       b_.reset(new volScalarField(tb()));
    }
    else
    {
        if(!solidToPoroFluid_.valid())
        {
            FatalErrorInFunction
                << "Cannot map the Biot coefficient because the solidToPoroFluid "
                   "mesh mapper has not been initialized." << nl
                << "This usually means non-shared mesh coupling was requested but "
                   "poroCouplingProperties did not create a valid mapper."
                << exit(FatalError);
        }
        const tmp<volScalarField> tbMapped(solidToPoroFluid().mapTgtToSrc(tb()));
        b_.reset(new volScalarField(tbMapped()));
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poroSolidInterface::poroSolidInterface(
    const word &type,
    Time &runTime,
    const word &region)
    : physicsModel(type, runTime),
      IOdictionary(
          IOobject(
              "poroCouplingProperties",
              runTime.constant(),
              runTime,
              IOobject::MUST_READ,
              IOobject::NO_WRITE)),
      runTime_(runTime),
      interactionProperties_(subDict(type + "Coeffs")),
      sharedMesh_(
          lookupOrAddDefault<Switch>(
              "sharedMesh", true)),
      porosityConstant_(lookupOrAddDefault<Switch>(
              "porosityConstant", true)),
      porosityConstantExplicit_(porosityConstant_
      ? Switch(true)
      : lookupOrAddDefault<Switch>(
              "porosityConstantExplicit", true)),
      solid_(),
      poroFluid_(),
      b_(),
      solidToPoroFluid_(),
      iterCtrl_(),
      intWork_()

{
    solid_ = solidModel::New(runTime, "solid");
    poroFluid_ = poroFluidModel::New(runTime, "poroFluid", sharedMesh_);

    if(!sharedMesh_)
    {
    const word mapMethodName(lookupOrAddDefault<word>("mapMethod","direct"));
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
    {
        FatalErrorInFunction
            << "Unknown mesh mapping method '" << mapMethodName << "' in "
            << "constant/poroCouplingProperties, sub-dictionary '" << type << "Coeffs'."
            << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_
            << exit(FatalError);
    }

    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    Switch consistent(
          lookupOrAddDefault<Switch>(
              "consistent", true));
        if(!consistent)
        {
            if(mapMethodName=="direct")
            {
                WarningInFunction
                    << "Direct mapping was selected for inconsistent region meshes. "
                    << "Switching to imMapNearest automatically."
                    << endl;
                mapMethod = meshToMesh::interpolationMethod::imMapNearest;
            }
        HashTable<word> patchMap;
        wordList cuttingPatches;

        readEntry("patchMap", patchMap);
        readEntry("cuttingPatches", cuttingPatches);
        solidToPoroFluid_.reset
                (
                    new meshToMesh
                    (
                        poroFluidMesh(),
                        solidMesh(),
                        mapMethod,
                        patchMap,
                        cuttingPatches
                    )
                );
        }
        else
        {
        solidToPoroFluid_.reset
                (
                    new meshToMesh
                    (
                        poroFluidMesh(),
                        solidMesh(),
                        mapMethod
                    )
                );
        }
    }

    intWork_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "internalWork",
                runTime.timeName(),
                solidMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh(),
            dimensionedScalar(dimPressure,0.0)
        )
    );

    makePoroFluidCouplingSource();
    Info << "Mesh is shared between Fields: " << sharedMesh_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::poroSolidInterface::~poroSolidInterface()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::poroSolidInterface> Foam::poroSolidInterface::New(
    Time &runTime,
    const word &region)
{
    word fsiTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the poroFluidModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary fsiProperties(
            IOobject(
                "poroCouplingProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE));

        fsiProperties.lookup("poroSolidInterface") >> fsiTypeName;
    }

    Info << "Selecting poroSolidInterface " << fsiTypeName << endl;
#if (OPENFOAM >= 2112)
        auto* ctorPtr = dictionaryConstructorTable(fsiTypeName);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                fsiTypeName,
                "poroSolidInterface",
                fsiTypeName,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }
#else
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fsiTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "poroSolidInterface::New(Time&, const word&)")
            << "Unknown poroSolidInterface type '" << fsiTypeName << "'."
            << endl
            << endl
            << "Valid poroSolidInterface types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

        auto* ctorPtr = cstrIter();
#endif

    return autoPtr<poroSolidInterface>(ctorPtr(runTime, region));
}

void Foam::poroSolidInterface::setDeltaT(Time &runTime)
{
    // For now, the poroFluid sets the time-step
    poroFluidRef().setDeltaT(runTime);
}


void Foam::poroSolidInterface::initializeFields()
{
}

void Foam::poroSolidInterface::prepareCouplingLoop()
{}

void Foam::poroSolidInterface::afterPorosityUpdate()
{}

void Foam::poroSolidInterface::afterFluidSolve()
{}

void Foam::poroSolidInterface::beforeSolidSolve()
{}

void Foam::poroSolidInterface::writeAdditionalFields(const Time&)
{}

void Foam::poroSolidInterface::movePoroFluidMesh()
{
    notImplemented("additionalPhysics mesh movement not yet implemented for fieldInteractions. Consider shared mesh for moving mesh interactions");
}

void Foam::poroSolidInterface::updateTotalFields()
{}

bool Foam::poroSolidInterface::evolveCouplingLoop()
{
    prepareCouplingLoop();
    initializeSolidHydraulicFields();

    SolverPerformance<vector>::debug = 0;

    couplingControl().reset();

    do
    {
        assembleCouplingTerms();

        if(!porosityConstantExplicit())
        {
            updatePorosityFromSolidDisplacement();
        }

        poroFluidRef().evolve();
        afterFluidSolve();

        if(!sharedMesh())
        {
            syncSolidHydraulicFields();
        }

        beforeSolidSolve();
        solidRef().evolve();
    }
    while(couplingControl().loop());

    clearCouplingTerms();
    couplingControl().write();

    Info << "Coupling Evolved" << endl;

    if(!porosityConstant())
    {
        updatePorosityFromSolidDisplacement();
    }

    return true;
}

bool Foam::poroSolidInterface::evolve()
{
    return evolveCouplingLoop();
}

void Foam::poroSolidInterface::writeFields(const Time &runTime)
{
    writeAdditionalFields(runTime);

    autoPtr<volSymmTensorField> DEpsilon;
    if(solid().incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            solidMesh().lookupObject<volTensorField>("grad(DD)");

        DEpsilon.reset(new volSymmTensorField(symm(gradDD)));
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            solidMesh().lookupObject<volTensorField>("grad(D)");

        DEpsilon.reset(new volSymmTensorField(symm(gradD)-symm(gradD.oldTime())));
    }

    intWork_.ref() =
        (solid().sigma()-solid().sigma().oldTime())
        && DEpsilon();
    intWork_().write();    

    poroFluidRef().writeFields(runTime);
    solidRef().writeFields(runTime);
}

void Foam::poroSolidInterface::end()
{
    this->IOobject::rename(this->IOobject::name()+".withDefaultValues");
    this->regIOobject::write();
    poroFluidRef().end();
}

// ************************************************************************* //
