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
#include "compatibilityFunctions.H"
#include <fstream>

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

Foam::tmp<Foam::fvMatrix<Foam::scalar>>
Foam::poroSolidInterface::fixedStressTimeMatrix
(
    const volScalarField& implicitCoupling,
    const volScalarField& p
) const
{
    // solid().U() is currently updated from fvc::ddt(solid().D()). Select
    // that exact scheme from the solid mesh, then instantiate its scalar
    // counterpart on the pressure mesh. This also works for separate meshes.
    ITstream& kinematicDdtScheme = solidMesh().ddtScheme
    (
        "ddt(" + solid().D().name() + ')'
    );

    return poroCouplingTerms::implicitCouplingMatrix
    (
        implicitCoupling,
        p,
        kinematicDdtScheme
    );
}

void Foam::poroSolidInterface::makeIterCtrl()
    {
        iterCtrl_.reset(new iterationControl(runTime_,interactionProperties_,"Pressure-Displacement"));
    }

void Foam::poroSolidInterface::resetSolidLinearSolverResiduals()
{
    firstSolidLinearInitialResidual_ = GREAT;
    lastSolidLinearInitialResidual_ = GREAT;
    solidLinearSolveCount_ = 0;
    solidLinearResidualCapability_ = false;

    dictionary& solverDict = const_cast<dictionary&>
    (
        solidMesh().data().solverPerformanceDict()
    );
    solverDict.remove("D");
    solverDict.remove("poroMechanicalFoamPreservesSolverPerformance");
}

void Foam::poroSolidInterface::captureSolidLinearSolverResiduals()
{
    const dictionary& solverDict =
        solidMesh().data().solverPerformanceDict();

    solidLinearResidualCapability_ =
        solverDict.lookupOrDefault<bool>
        (
            "poroMechanicalFoamPreservesSolverPerformance",
            false
        );

    if
    (
        !solidLinearResidualCapability_
     || !solverDict.found("D")
    )
    {
        return;
    }

    const List<SolverPerformance<vector>> history
    (
        solverDict.lookup("D")
    );

    if (history.empty())
    {
        return;
    }

    firstSolidLinearInitialResidual_ =
        cmptMax(mag(history.first().initialResidual()));
    lastSolidLinearInitialResidual_ =
        cmptMax(mag(history.last().initialResidual()));
    solidLinearSolveCount_ = history.size();
}

Foam::scalar Foam::poroSolidInterface::solidLinearSolverInitialResidual
(
    const bool first
) const
{
    if (!solidLinearResidualCapability_)
    {
        FatalErrorInFunction
            << "Solid linear-solver residual history is unavailable. "
            << "This solids4foam build does not advertise the "
            << "poroMechanicalFoam residual-preservation capability."
            << exit(FatalError);
    }

    if (!solidLinearSolveCount_)
    {
        FatalErrorInFunction
            << "The solid model recorded no D solve in its current evolve() call."
            << exit(FatalError);
    }

    return first
      ? firstSolidLinearInitialResidual_
      : lastSolidLinearInitialResidual_;
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
    tmp<volScalarField> tStabil(poroCouplingTerms::fixedStressStabil(b, impK));
    tmpRef(tStabil) *= fixedStressStabilScale_;
    return tStabil;
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

    fixedStressStabilField.ref() *= fixedStressStabilScale_;
}

void Foam::poroSolidInterface::printFixedStressDiagnostic
(
    const word& stage,
    const dimensionSet& matrixDimensions
) const
{
    const tmp<volScalarField> tImplicitCoupling(implicitCouplingDtoP());
    const tmp<volScalarField> tExplicitCoupling(explicitCouplingDtoP());

    const volScalarField& p = pField();

    volScalarField zeroExplicitCoupling
    (
        IOobject
        (
            "zeroExplicitCoupling",
            runTime().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p.mesh(),
        dimensionedScalar("zero", tExplicitCoupling().dimensions(), 0.0)
    );

    fvScalarMatrix fixedStressOption(p, matrixDimensions);
    const tmp<fvScalarMatrix> tFixedStressTime
    (
        fixedStressTimeMatrix(tImplicitCoupling(), p)
    );

    poroCouplingTerms::addCouplingSource
    (
        fixedStressOption,
        p,
        pCouplingRef(),
        tFixedStressTime(),
        zeroExplicitCoupling
    );

    const scalarField& V = p.mesh().V();
    const scalarField finalResidualByV
    (
        (
            -fixedStressOption.diag()*p.primitiveField()
          + fixedStressOption.source()
        )/V
    );

    if (stage == "after pressure solve" && Pstream::master())
    {
        const fileName csvPath(runTime().path()/"fixedStressDiagnostic.csv");

        std::ifstream existing(csvPath.c_str());
        const bool writeHeader = !existing.good();
        existing.close();

        std::ofstream csv(csvPath.c_str(), std::ios::app);
        if (writeHeader)
        {
            csv
                << "time,outer-iteration,inner-iteration,"
                << "fixed-stress-min-contribution,"
                << "fixed-stress-avg-contribution,"
                << "fixed-stress-max-contribution\n";
        }

        const label outerIteration = iterCtrl_.valid()
          ? iterCtrl_().index() + 1
          : -1;
        const label innerIteration =
            const_cast<poroSolidInterface&>(*this)
           .poroFluidRef().iterCtrl().index();

        csv
            << runTime().timeName() << ','
            << outerIteration << ','
            << innerIteration << ','
            << gMin(finalResidualByV) << ','
            << average(finalResidualByV) << ','
            << gMax(finalResidualByV) << '\n';
    }

    Info<< "Fixed-stress contribution (" << stage << ") at time "
        << runTime().timeName() << " min/avg/max: "
        << gMin(finalResidualByV) << " / "
        << average(finalResidualByV) << " / "
        << gMax(finalResidualByV) << endl;
}

void Foam::poroSolidInterface::printExplicitCouplingDiagnostic
(
    const word& stage,
    const dimensionSet& matrixDimensions
) const
{
    const tmp<volScalarField> tImplicitCoupling(implicitCouplingDtoP());
    const tmp<volScalarField> tExplicitCoupling(explicitCouplingDtoP());

    const volScalarField& p = pField();

    volScalarField zeroImplicitCoupling
    (
        IOobject
        (
            "zeroImplicitCoupling",
            runTime().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p.mesh(),
        dimensionedScalar("zero", tImplicitCoupling().dimensions(), 0.0)
    );

    fvScalarMatrix explicitOption(p, matrixDimensions);
    const tmp<fvScalarMatrix> tFixedStressTime
    (
        fixedStressTimeMatrix(zeroImplicitCoupling, p)
    );

    poroCouplingTerms::addCouplingSource
    (
        explicitOption,
        p,
        pCouplingRef(),
        tFixedStressTime(),
        tExplicitCoupling()
    );

    const scalarField& V = p.mesh().V();
    const scalarField optionDiagByV(explicitOption.diag()/V);
    const scalarField optionSourceByV(explicitOption.source()/V);
    const scalarField finalDiagByV(-explicitOption.diag()/V);
    const scalarField finalSourceByV(-explicitOption.source()/V);

    Info<< "Explicit coupling diagnostic (" << stage << ") at time "
        << runTime().timeName() << nl
        << "  explicit source nDot/with corrections min/avg/max: "
        << gMin(tExplicitCoupling().primitiveField()) << " / "
        << average(tExplicitCoupling().primitiveField()) << " / "
        << gMax(tExplicitCoupling().primitiveField()) << nl
        << "  option diag/V min/avg/max: "
        << gMin(optionDiagByV) << " / "
        << average(optionDiagByV) << " / "
        << gMax(optionDiagByV) << nl
        << "  option source/V min/avg/max: "
        << gMin(optionSourceByV) << " / "
        << average(optionSourceByV) << " / "
        << gMax(optionSourceByV) << nl
        << "  final diag/V min/avg/max after lhs==fvOptions: "
        << gMin(finalDiagByV) << " / "
        << average(finalDiagByV) << " / "
        << gMax(finalDiagByV) << nl
        << "  final source/V min/avg/max after lhs==fvOptions: "
        << gMin(finalSourceByV) << " / "
        << average(finalSourceByV) << " / "
        << gMax(finalSourceByV) << nl;

    if (p.mesh().objectRegistry::foundObject<volScalarField>("nDot"))
    {
        const volScalarField& nDotField =
            p.mesh().objectRegistry::lookupObject<volScalarField>("nDot");

        const volScalarField divUFromNDot
        (
            "divUFromNDot",
            nDotField
          / max
            (
                b(),
                dimensionedScalar("smallBiot", dimless, SMALL)
            )
        );

        Info<< "  nDot=b*div(U) min/avg/max: "
            << gMin(nDotField.primitiveField()) << " / "
            << average(nDotField.primitiveField()) << " / "
            << gMax(nDotField.primitiveField()) << nl
            << "  div(U) inferred from nDot/b min/avg/max: "
            << gMin(divUFromNDot.primitiveField()) << " / "
            << average(divUFromNDot.primitiveField()) << " / "
            << gMax(divUFromNDot.primitiveField()) << endl;
    }
    else
    {
        Info<< "  nDot field not registered on the fluid mesh" << endl;
    }
}

void Foam::poroSolidInterface::storeCouplingPressureReference()
{
    const volScalarField& p = pField();

    if (!pCouplingRef_.valid())
    {
        pCouplingRef_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pCouplingRef",
                    runTime().timeName(),
                    p.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                p
            )
        );
    }
    else
    {
        pCouplingRef_.ref() = p;
    }
}

const Foam::volScalarField& Foam::poroSolidInterface::pCouplingRef() const
{
    if (!pCouplingRef_.valid())
    {
        FatalErrorInFunction
            << "The coupling pressure reference has not been stored yet. "
            << "storeCouplingPressureReference() runs at the start of each "
            << "coupling iteration; coupling terms cannot be assembled "
            << "outside the coupling loop."
            << exit(FatalError);
    }

    return pCouplingRef_();
}

void Foam::poroSolidInterface::storeCouplingPressureReferenceAsPrevIter()
{
    volScalarField& p = poroFluidRef().pField();

    volScalarField pCurrent
    (
        IOobject
        (
            "pCouplingResidualCurrent",
            runTime().timeName(),
            p.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        p
    );

    p = pCouplingRef();
    p.storePrevIter();
    p = pCurrent;
}

void Foam::poroSolidInterface::addPressureEquationTerms
(
    fvMatrix<scalar>& eqn
) const
{
    const tmp<volScalarField> tImplicitCoupling(implicitCouplingDtoP());
    const tmp<volScalarField> tExplicitCoupling(explicitCouplingDtoP());

    if (fixedStressDiagnostic_)
    {
        printFixedStressDiagnostic("before pressure solve", eqn.dimensions());
        printExplicitCouplingDiagnostic("before pressure solve", eqn.dimensions());
    }

    const tmp<fvScalarMatrix> tFixedStressTime
    (
        fixedStressTimeMatrix(tImplicitCoupling(), pField())
    );

    poroCouplingTerms::addCouplingSource
    (
        eqn,
        pField(),
        pCouplingRef(),
        tFixedStressTime(),
        tExplicitCoupling()
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
      fixedStressStabilScale_(
          lookupOrAddDefault<scalar>(
              "fixedStressStabilScale", 1.5)),
      fixedStressDiagnostic_(
          lookupOrAddDefault<Switch>(
              "fixedStressDiagnostic", false)),
      pressureModeDiagnosticFields_(
          lookupOrAddDefault<wordList>
          (
              "pressureModeDiagnosticFields",
              wordList()
          )),
      solid_(),
      poroFluid_(),
      pressureModeDiagnostics_(),
      previousPressureModeAmplitudes_
      (
          pressureModeDiagnosticFields_.size(),
          0.0
      ),
      previousPressureDiagnostic_(),
      havePreviousPressureModeAmplitude_(false),
      pressureModeDiagnosticIteration_(0),
      b_(),
      solidToPoroFluid_(),
      iterCtrl_(),
      intWork_(),
      firstSolidLinearInitialResidual_(GREAT),
      lastSolidLinearInitialResidual_(GREAT),
      solidLinearSolveCount_(0),
      solidLinearResidualCapability_(false),
      pCouplingRef_()

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
    Info << "Fixed-stress stabilization scale: "
         << fixedStressStabilScale_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::poroSolidInterface::~poroSolidInterface()
{
    if(sharedMesh() && solid_.valid())
    {
        objectRegistry& solidRegistry =
            const_cast<objectRegistry&>
            (
                static_cast<const objectRegistry&>(solidMesh())
            );

        forAllConstIter(HashSet<word>, sharedSolidRegistryFieldNames_, iter)
        {
            solidRegistry.checkOut(iter.key());
        }
    }
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
{
    if (fixedStressDiagnostic_)
    {
        printFixedStressDiagnostic("after pressure solve", dimVolume/dimTime);
        printExplicitCouplingDiagnostic("after pressure solve", dimVolume/dimTime);
    }

    printPressureModeDiagnostic();
}

void Foam::poroSolidInterface::printPressureModeDiagnostic()
{
    if (pressureModeDiagnosticFields_.empty())
    {
        return;
    }

    if (pressureModeDiagnostics_.empty())
    {
        pressureModeDiagnostics_.setSize(pressureModeDiagnosticFields_.size());
        forAll(pressureModeDiagnosticFields_, modeI)
        {
            pressureModeDiagnostics_.set
            (
                modeI,
                new volScalarField
                (
                    IOobject
                    (
                        pressureModeDiagnosticFields_[modeI],
                        runTime_.constant(),
                        "poroFluid",
                        runTime_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    poroFluidMesh()
                )
            );
        }
    }

    const fvMesh& mesh = poroFluidMesh();
    const scalarField& volumes = mesh.V();
    const scalarField& pressure = pField().primitiveField();
    const scalar totalVolume = gSum(volumes);
    const scalar meanPressure = gSum(volumes*pressure)/totalVolume;
    scalarField amplitudes(pressureModeDiagnosticFields_.size(), 0.0);

    forAll(pressureModeDiagnosticFields_, modeI)
    {
        const scalarField& mode =
            pressureModeDiagnostics_[modeI].primitiveField();
        const scalar denominator = gSum(volumes*sqr(mode));

        if (denominator <= VSMALL)
        {
            FatalErrorInFunction
                << "Pressure mode diagnostic field '"
                << pressureModeDiagnosticFields_[modeI] << "' has zero norm"
                << exit(FatalError);
        }

        amplitudes[modeI] =
            gSum(volumes*mode*(pressure - meanPressure))/denominator;
    }

    ++pressureModeDiagnosticIteration_;
    if (havePreviousPressureModeAmplitude_)
    {
        const scalarField pressureIncrement =
            pressure - previousPressureDiagnostic_;
        scalarField checkerboardIncrement(pressure.size(), 0.0);
        forAll(pressureModeDiagnosticFields_, modeI)
        {
            checkerboardIncrement +=
                (amplitudes[modeI] - previousPressureModeAmplitudes_[modeI])
               *pressureModeDiagnostics_[modeI].primitiveField();
        }
        const scalarField complementIncrement =
            pressureIncrement - checkerboardIncrement;
        const scalar totalIncrementRms = Foam::sqrt
        (
            gSum(volumes*sqr(pressureIncrement))/totalVolume
        );
        const scalar checkerboardIncrementRms = Foam::sqrt
        (
            gSum(volumes*sqr(checkerboardIncrement))/totalVolume
        );
        const scalar complementIncrementRms = Foam::sqrt
        (
            gSum(volumes*sqr(complementIncrement))/totalVolume
        );

        const scalar incrementNorm =
            Foam::sqrt(sum(sqr(amplitudes - previousPressureModeAmplitudes_)));
        Info<< "Pressure mode subspace diagnostic: iteration "
            << pressureModeDiagnosticIteration_
            << ", incrementNorm " << incrementNorm << endl;
        Info<< "Pressure increment decomposition: iteration "
            << pressureModeDiagnosticIteration_
            << ", totalRms " << totalIncrementRms
            << ", checkerboardRms " << checkerboardIncrementRms
            << ", complementRms " << complementIncrementRms
            << ", checkerboardFraction "
            << checkerboardIncrementRms/max(totalIncrementRms, VSMALL)
            << endl;
    }
    else
    {
        Info<< "Pressure mode subspace diagnostic: iteration "
            << pressureModeDiagnosticIteration_
            << ", incrementNorm nan" << endl;
    }

    forAll(pressureModeDiagnosticFields_, modeI)
    {
        Info<< "Pressure mode diagnostic "
            << pressureModeDiagnosticFields_[modeI]
            << ": iteration " << pressureModeDiagnosticIteration_
            << ", amplitude " << amplitudes[modeI];
        if (havePreviousPressureModeAmplitude_)
        {
            Info<< ", increment "
                << amplitudes[modeI]
                 - previousPressureModeAmplitudes_[modeI];
        }
        else
        {
            Info<< ", increment nan";
        }
        Info<< endl;
    }

    previousPressureModeAmplitudes_ = amplitudes;
    previousPressureDiagnostic_ = pressure;
    havePreviousPressureModeAmplitude_ = true;
}

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

    previousPressureModeAmplitudes_ = 0.0;
    previousPressureDiagnostic_.clear();
    havePreviousPressureModeAmplitude_ = false;
    pressureModeDiagnosticIteration_ = 0;

    SolverPerformance<vector>::debug = 0;

    couplingControl().reset();

    do
    {
        // Freeze the fixed-stress pressure reference together with the
        // coupling terms for this coupling iteration
        storeCouplingPressureReference();

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
        resetSolidLinearSolverResiduals();
        solidRef().evolve();
        captureSolidLinearSolverResiduals();

        // The fluid sub-loop uses the pressure field's single prevIter slot
        // for relaxation/nonlinear sources. Restore it to the same reference
        // used by fixed-stress before the outer coupling convergence check.
        storeCouplingPressureReferenceAsPrevIter();
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
