#include "fvCFD.H"
#include "poroCouplingTerms.H"
#include "poroCouplingRegistry.H"
#include "poroPressureUnits.H"
#include "UniformDimensionedField.H"
#include "effectiveStressModel.H"
#include "varSatPoroMechanicalLawTerms.H"
#include "residualOperation.H"
#include "deltaVf.H"
#include "iterationControl.H"
#include "LinearSolverRes.H"
#include "SolverPerformance.H"
#include "MassBalanceTerms.H"

using namespace Foam;

namespace
{
static label failures = 0;

scalar relErr(const scalar actual, const scalar expected)
{
    return mag(actual - expected)/max(scalar(1), mag(expected));
}

void checkNear
(
    const word& name,
    const scalar actual,
    const scalar expected,
    const scalar tolerance = 1e-11
)
{
    if (relErr(actual, expected) > tolerance)
    {
        ++failures;
        Info<< "FAIL " << name
            << ": actual=" << actual
            << " expected=" << expected
            << " relErr=" << relErr(actual, expected)
            << " tolerance=" << tolerance << nl;
    }
    else
    {
        Info<< "PASS " << name << nl;
    }
}

void checkTrue(const word& name, const bool condition)
{
    if (!condition)
    {
        ++failures;
        Info<< "FAIL " << name << nl;
    }
    else
    {
        Info<< "PASS " << name << nl;
    }
}

bool containsWord(const wordList& words, const word& value)
{
    forAll(words, wordI)
    {
        if (words[wordI] == value)
        {
            return true;
        }
    }

    return false;
}

template<class GeoField>
scalar maxMag(const GeoField& field)
{
    scalar result = gMax(mag(field.primitiveField()));

    forAll(field.boundaryField(), patchI)
    {
        result = max(result, gMax(mag(field.boundaryField()[patchI])));
    }

    return result;
}

void testNDotRigidTranslation(const fvMesh& mesh)
{
    volScalarField couplingCoeff
    (
        IOobject("couplingCoeff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("couplingCoeff", dimless, 0.42),
        "zeroGradient"
    );

    volVectorField U
    (
        IOobject("U", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("U", dimVelocity, vector(0.13, -0.07, 0.02)),
        "zeroGradient"
    );

    const tmp<volScalarField> tnDot(poroCouplingTerms::nDot(couplingCoeff, U));

    checkNear("nDot field name", tnDot().name() == "nDot" ? 1.0 : 0.0, 1.0);
    checkNear("nDot is zero for rigid translation", maxMag(tnDot()), 0.0);
}

void testFixedStressStabil(const fvMesh& mesh)
{
    volScalarField biotCoeff
    (
        IOobject("biotCoeff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("biotCoeff", dimless, 0.3),
        "zeroGradient"
    );

    volScalarField impK
    (
        IOobject("impK", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("impK", dimPressure, 2.0e5),
        "zeroGradient"
    );

    const tmp<volScalarField> tStabil
    (
        poroCouplingTerms::fixedStressStabil(biotCoeff, impK)
    );

    checkNear("fixedStressStabil b^2/K golden", tStabil()[0], 4.5e-7);
    checkTrue("fixedStressStabil positive", tStabil()[0] > 0.0);
    checkTrue
    (
        "fixedStressStabil dimensions are inverse pressure",
        tStabil().dimensions() == dimless/dimPressure
    );
}

void testVarSatWeightingDecision(const fvMesh& mesh)
{
    const scalar S = 0.25;
    const scalar b = 0.8;

    volScalarField saturationWeightedBiot
    (
        IOobject("saturationWeightedBiot", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("saturationWeightedBiot", dimless, S*b),
        "zeroGradient"
    );

    volScalarField biotCoeff
    (
        IOobject("biotCoeffForVarSat", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("biotCoeffForVarSat", dimless, b),
        "zeroGradient"
    );

    volScalarField impK
    (
        IOobject("impKForVarSat", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("impKForVarSat", dimPressure, 2.0e5),
        "zeroGradient"
    );

    volVectorField U
    (
        IOobject("affineU", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("affineU", dimVelocity, Zero),
        "calculated"
    );

    forAll(U, cellI)
    {
        U[cellI] = mesh.C()[cellI];
    }
    forAll(U.boundaryField(), patchI)
    {
        U.boundaryFieldRef()[patchI] == mesh.Cf().boundaryField()[patchI];
    }

    const tmp<volScalarField> tDivU(fvc::div(U));
    const tmp<volScalarField> tnDot
    (
        poroCouplingTerms::nDot(saturationWeightedBiot, U)
    );
    const tmp<volScalarField> tStabil
    (
        poroCouplingTerms::fixedStressStabil(biotCoeff, impK)
    );

    checkNear("varSat affine velocity divergence golden", tDivU()[0], 3.0);
    checkNear("varSat nDot uses S*b weighting", tnDot()[0], 0.6);
    checkNear("varSat fixedStress uses b, not S*b", tStabil()[0], 3.2e-6);
    checkTrue
    (
        "varSat fixedStress would be smaller with S*b",
        tStabil()[0] > sqr(S*b)/2.0e5
    );
}

void testCouplingFieldLifecycle(fvMesh& mesh)
{
    volScalarField explicitCoeff
    (
        IOobject("lifecycleExplicitCoeff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("lifecycleExplicitCoeff", dimless, 0.2),
        "zeroGradient"
    );

    volScalarField stabilCoeff
    (
        IOobject("lifecycleStabilCoeff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("lifecycleStabilCoeff", dimless, 0.8),
        "zeroGradient"
    );

    volScalarField impK
    (
        IOobject("lifecycleImpK", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("lifecycleImpK", dimPressure, 2.0e5),
        "zeroGradient"
    );

    volVectorField U
    (
        IOobject("lifecycleAffineU", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("lifecycleAffineU", dimVelocity, Zero),
        "calculated"
    );

    forAll(U, cellI)
    {
        U[cellI] = mesh.C()[cellI];
    }
    forAll(U.boundaryField(), patchI)
    {
        U.boundaryFieldRef()[patchI] == mesh.Cf().boundaryField()[patchI];
    }

    autoPtr<volScalarField> nDotField;
    autoPtr<volScalarField> fixedStressField;

    poroCouplingTerms::updateCouplingFields
    (
        explicitCoeff,
        stabilCoeff,
        impK,
        U,
        nDotField,
        fixedStressField
    );

    objectRegistry& registry =
        const_cast<objectRegistry&>(static_cast<const objectRegistry&>(mesh));
    const volScalarField* firstNDotPtr = &nDotField();

    checkTrue("coupling lifecycle creates nDot", nDotField.valid());
    checkTrue("coupling lifecycle creates fixedStress", fixedStressField.valid());
    checkTrue("coupling lifecycle registers nDot", registry.foundObject<volScalarField>("nDot"));
    checkTrue
    (
        "coupling lifecycle registry owns nDot instance",
        &registry.lookupObject<volScalarField>("nDot") == firstNDotPtr
    );
    checkNear("coupling lifecycle first nDot uses explicit coeff", nDotField()[0], 0.6);
    checkNear("coupling lifecycle first fixedStress uses stabil coeff", fixedStressField()[0], 3.2e-6);

    explicitCoeff = dimensionedScalar("updatedExplicitCoeff", dimless, 0.5);
    stabilCoeff = dimensionedScalar("updatedStabilCoeff", dimless, 0.4);

    forAll(U, cellI)
    {
        U[cellI] = 2.0*mesh.C()[cellI];
    }
    forAll(U.boundaryField(), patchI)
    {
        U.boundaryFieldRef()[patchI] == 2.0*mesh.Cf().boundaryField()[patchI];
    }

    poroCouplingTerms::updateCouplingFields
    (
        explicitCoeff,
        stabilCoeff,
        impK,
        U,
        nDotField,
        fixedStressField
    );

    checkTrue("coupling lifecycle reuses nDot field", &nDotField() == firstNDotPtr);
    checkTrue
    (
        "coupling lifecycle registry still points at nDot",
        &registry.lookupObject<volScalarField>("nDot") == firstNDotPtr
    );
    checkNear("coupling lifecycle updates nDot in place", nDotField()[0], 3.0);
    checkNear("coupling lifecycle updates fixedStress from stabil coeff", fixedStressField()[0], 8.0e-7);
}

void testRelativeAccelerationFlux(const fvMesh& mesh)
{
    surfaceScalarField kf
    (
        IOobject("kf", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("kf", dimVelocity, 0.06)
    );

    volVectorField a
    (
        IOobject("a", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("a", dimAcceleration, vector(0, 0, -9.81)),
        "zeroGradient"
    );

    const dimensionedScalar magG("magG", dimAcceleration, 9.81);
    const tmp<surfaceVectorField> tq
    (
        poroCouplingTerms::relativeAccelerationFlux(kf, a, magG)
    );

    checkNear("relativeAccelerationFlux boundary z", tq().boundaryField()[0][0].z(), -0.06);
    checkTrue("relativeAccelerationFlux has velocity dimensions", tq().dimensions() == dimVelocity);
}

void testExplicitCouplingSourceSign(const fvMesh& mesh)
{
    volScalarField nDot
    (
        IOobject("nDotInput", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("nDotInput", dimless, 2.5),
        "zeroGradient"
    );

    surfaceVectorField qRelAcc
    (
        IOobject("qRelAccInput", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh.Cf()
    );

    const tmp<volScalarField> tSource
    (
        poroCouplingTerms::explicitCouplingSource(nDot, qRelAcc)
    );

    checkNear("explicitCouplingSource subtracts acceleration divergence", tSource()[0], -0.5);
    checkTrue("explicitCouplingSource keeps nDot dimensions", tSource().dimensions() == nDot.dimensions());
}

void testFvOptionCouplingRates(const fvMesh& mesh)
{
    volScalarField implicitCoeff
    (
        IOobject("implicitCoeff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("implicitCoeff", dimless/dimPressure, 8.0e-7),
        "zeroGradient"
    );

    volScalarField explicitSource
    (
        IOobject("explicitSource", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("explicitSource", dimless/dimTime, 1.75),
        "zeroGradient"
    );

    const dimensionedScalar deltaT("deltaT", dimTime, 0.25);
    const tmp<volScalarField> tImplicitRate
    (
        poroCouplingTerms::implicitCouplingRate(implicitCoeff, deltaT)
    );
    const tmp<volScalarField> tExplicitRate
    (
        poroCouplingTerms::explicitCouplingRate(explicitSource)
    );

    checkNear("fvOption implicit coupling is scaled by 1/deltaT", tImplicitRate()[0], 3.2e-6);
    checkTrue
    (
        "fvOption implicit rate dimensions include inverse time",
        tImplicitRate().dimensions() == dimless/(dimPressure*dimTime)
    );
    checkNear("fvOption explicit coupling enters with negative sign", tExplicitRate()[0], -1.75);
    checkTrue
    (
        "fvOption explicit rate keeps source dimensions",
        tExplicitRate().dimensions() == explicitSource.dimensions()
    );
}

void testFvOptionMatrixAssembly(const fvMesh& mesh)
{
    volScalarField p
    (
        IOobject("matrixAssemblyP", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("matrixAssemblyP", dimPressure, 10.0),
        "zeroGradient"
    );

    volScalarField implicitCoeff
    (
        IOobject("matrixImplicitCoeff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("matrixImplicitCoeff", dimless/dimPressure, 8.0e-7),
        "zeroGradient"
    );

    volScalarField explicitSource
    (
        IOobject("matrixExplicitSource", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("matrixExplicitSource", dimless/dimTime, 1.75),
        "zeroGradient"
    );

    fvScalarMatrix eqn(p, dimVol/dimTime);
    const dimensionedScalar deltaT("matrixDeltaT", dimTime, 0.25);

    poroCouplingTerms::addCouplingSource
    (
        eqn,
        p,
        implicitCoeff,
        explicitSource,
        deltaT
    );

    checkNear("fvOption matrix diagonal gets implicit rate", eqn.diag()[0], 3.2e-6);
    checkNear("fvOption matrix source gets negative explicit rate", eqn.source()[0], 1.75);
    checkTrue("fvOption matrix dimensions are volume per time", eqn.dimensions() == dimVol/dimTime);
}

void testResidualOperations()
{
    scalarField values(3);
    values[0] = 3.0;
    values[1] = 4.0;
    values[2] = -12.0;

    List<scalar> listValues(3);
    listValues[0] = 3.0;
    listValues[1] = 4.0;
    listValues[2] = -12.0;

    autoPtr<residualOperation> maxOp(residualOperation::New("max"));
    autoPtr<residualOperation> sumOp(residualOperation::New("sum"));
    autoPtr<residualOperation> l2Op(residualOperation::New("L2"));
    autoPtr<residualOperation> rmsOp(residualOperation::New("RMS"));

    checkNear("residualOperation max scalarField", maxOp->operation(values), 4.0);
    checkNear("residualOperation max List", maxOp->operation(listValues), 4.0);
    checkNear("residualOperation sum scalarField keeps signs", sumOp->operation(values), -5.0);
    checkNear("residualOperation sum List keeps signs", sumOp->operation(listValues), -5.0);
    checkNear("residualOperation L2 scalarField", l2Op->operation(values), 13.0);
    checkNear("residualOperation L2 List", l2Op->operation(listValues), 13.0);
    checkNear("residualOperation RMS scalarField", rmsOp->operation(values), Foam::sqrt(169.0/3.0));
    checkNear("residualOperation RMS List", rmsOp->operation(listValues), Foam::sqrt(169.0/3.0));

    scalarField emptyField(0);
    List<scalar> emptyList(0);
    checkNear("residualOperation RMS empty scalarField", rmsOp->operation(emptyField), 0.0);
    checkNear("residualOperation RMS empty List", rmsOp->operation(emptyList), 0.0);
}

void testDeltaVfResiduals(fvMesh& mesh)
{
    volScalarField trackedAbs
    (
        IOobject("trackedAbsDelta", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("trackedAbsDelta", dimLength, 2.0),
        "zeroGradient"
    );

    ITstream absStream("max 1e-12");
    deltaVf absResidual(mesh.time(), trackedAbs.name(), absStream, false);

    checkTrue("deltaVf finds registered scalar field", deltaVf::fieldExists(mesh.time(), trackedAbs.name()));
    checkTrue
    (
        "deltaVf candidates include scalar field",
        containsWord(deltaVf::candidateFieldNames(mesh.time()), trackedAbs.name())
    );
    checkTrue("deltaVf rejects missing field", !deltaVf::fieldExists(mesh.time(), "missingDeltaField"));

    checkNear("deltaVf first scalar residual stores baseline", absResidual.calcResidual(), 0.0);

    trackedAbs = dimensionedScalar("trackedAbsDelta", dimLength, 5.5);

    checkNear("deltaVf scalar absolute max delta", absResidual.calcResidual(), 3.5);

    volScalarField trackedRel
    (
        IOobject("trackedRelDelta", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("trackedRelDelta", dimPressure, 2.0),
        "zeroGradient"
    );

    ITstream relStream("max rel 1e-12");
    deltaVf relResidual(mesh.time(), trackedRel.name(), relStream, false);

    checkNear("deltaVf first relative residual stores baseline", relResidual.calcResidual(), 0.0);

    trackedRel = dimensionedScalar("trackedRelDelta", dimPressure, 4.0);

    checkNear
    (
        "deltaVf relative delta uses current magnitude denominator",
        relResidual.calcResidual(),
        0.5
    );

    volVectorField trackedVector
    (
        IOobject("trackedVectorDelta", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("trackedVectorDelta", dimVelocity, vector(3.0, 4.0, 0.0)),
        "zeroGradient"
    );

    ITstream vectorStream("max 1e-12");
    deltaVf vectorResidual(mesh.time(), trackedVector.name(), vectorStream, false);

    checkNear("deltaVf first vector residual stores baseline", vectorResidual.calcResidual(), 0.0);

    trackedVector = dimensionedVector("trackedVectorDelta", dimVelocity, vector(0.0, 0.0, 12.0));

    checkNear("deltaVf vector residual tracks magnitude delta", vectorResidual.calcResidual(), 7.0);

    surfaceScalarField trackedSurface
    (
        IOobject("trackedSurfaceDelta", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("trackedSurfaceDelta", dimVelocity, 0.25)
    );

    ITstream surfaceStream("max 1e-12");
    deltaVf surfaceResidual(mesh.time(), trackedSurface.name(), surfaceStream, false);

    checkNear("deltaVf first surface residual stores baseline", surfaceResidual.calcResidual(), 0.0);

    trackedSurface = dimensionedScalar("trackedSurfaceDelta", dimVelocity, 1.5);

    checkNear("deltaVf surface scalar residual tracks delta", surfaceResidual.calcResidual(), 1.25);
}

void testIterationControlDictionary(fvMesh& mesh)
{
    Time& runTime = const_cast<Time&>(mesh.time());

    dictionary futureDict
    (
        IStringStream
        (
            "enabled true;"
            "timeStart 100;"
            "iterations 5;"
        )()
    );

    iterationControl inactiveControl(runTime, futureDict, "futureLoop");

    checkNear("iterationControl disables future time window", inactiveControl.nCycles(), 0);
    checkTrue("iterationControl future window is inactive", !inactiveControl.status());

    volScalarField trackedControl
    (
        IOobject("trackedControlDelta", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("trackedControlDelta", dimLength, 1.0),
        "zeroGradient"
    );

    dictionary activeDict
    (
        IStringStream
        (
            "enabled true;"
            "iterations 4;"
            "interval 2;"
            "infoFrequency 1000;"
            "writeResidualField false;"
            "convergence"
            "{"
            "    trackedControlDelta max 1e-9;"
            "}"
        )()
    );

    iterationControl activeControl(runTime, activeDict, "activeLoop");

    checkNear("iterationControl reads iteration count", activeControl.nCycles(), 4);
    checkNear("iterationControl reads loop interval", activeControl.interval(), 2);
    checkTrue("iterationControl detects active status", activeControl.status());
    checkTrue
    (
        "iterationControl does not require mass balance for field residual",
        !activeControl.requiresMassBalanceResidual()
    );
    checkNear("iterationControl creates one residual", activeControl.residuals().size(), 1);
    checkNear
    (
        "iterationControl residual first call stores baseline",
        activeControl.residuals()[0].calcResidual(),
        0.0
    );

    trackedControl = dimensionedScalar("trackedControlDelta", dimLength, 3.25);

    checkNear
    (
        "iterationControl residual uses configured delta field",
        activeControl.residuals()[0].calcResidual(),
        2.25
    );
}

void testLinearSolverResiduals(fvMesh& mesh)
{
    mesh.data().solverPerformanceDict().clear();

    volScalarField scalarSolveField
    (
        IOobject("linearSolverScalar", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("linearSolverScalar", dimPressure, 1.0),
        "zeroGradient"
    );

    List<SolverPerformance<scalar>> scalarPerformance(1);
    scalarPerformance[0] =
        SolverPerformance<scalar>
        (
            "unitSolver",
            scalarSolveField.name(),
            0.125,
            1.0e-8,
            3,
            true,
            false
        );
    mesh.data().solverPerformanceDict().set(scalarSolveField.name(), scalarPerformance);

    ITstream scalarStream("linearSolver 1e-4");
    LinearSolverRes scalarResidual(mesh.time(), scalarSolveField.name(), scalarStream, false);

    checkNear("linearSolver scalar tolerance parse", scalarResidual.tolerance(), 1.0e-4);
    checkNear("linearSolver scalar initial residual", scalarResidual.calcResidual(), 0.125);

    volVectorField vectorSolveField
    (
        IOobject("linearSolverVector", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("linearSolverVector", dimVelocity, vector::zero),
        "zeroGradient"
    );

    List<SolverPerformance<vector>> vectorPerformance(1);
    vectorPerformance[0] =
        SolverPerformance<vector>
        (
            "unitSolver",
            vectorSolveField.name(),
            vector(0.05, 0.35, 0.2),
            vector(1.0e-9, 2.0e-9, 3.0e-9),
            Vector<label>(1, 2, 3),
            true,
            false
        );
    mesh.data().solverPerformanceDict().set(vectorSolveField.name(), vectorPerformance);

    ITstream vectorStream("linearSolver show");
    autoPtr<iterationResidual> vectorResidual
    (
        iterationResidual::New
        (
            mesh.time(),
            vectorSolveField.name(),
            vectorStream,
            false
        )
    );

    checkNear("linearSolver show tolerance disables convergence check", vectorResidual->tolerance(), -1.0);
    checkNear("linearSolver vector residual uses max component", vectorResidual->calcResidual(), 0.35);

    mesh.data().solverPerformanceDict().clear();
}

void testMassBalanceTerms(const fvMesh& mesh)
{
    volScalarField massBalance
    (
        IOobject("unitMassBalanceResidual", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("unitMassBalanceResidual", dimVolume, 4.0),
        "zeroGradient"
    );

    massBalance.oldTime();
    massBalance.storeOldTime();
    massBalance = dimensionedScalar("unitMassBalanceResidual", dimVolume, -1.0);

    const scalarField absoluteResidual
    (
        MassBalanceTerms::residualValues(massBalance, false)
    );
    const scalarField relativeResidual
    (
        MassBalanceTerms::residualValues(massBalance, true)
    );

    checkNear("MassBalance absolute residual uses magnitude", absoluteResidual[0], 1.0);
    checkNear("MassBalance relative residual uses stored old-time denominator", relativeResidual[0], 0.25);

    massBalance.storeOldTime();
    massBalance = dimensionedScalar("unitMassBalanceResidual", dimVolume, 0.0);

    const scalarField zeroRelativeResidual
    (
        MassBalanceTerms::residualValues(massBalance, true)
    );

    checkNear("MassBalance zero relative residual stays zero", zeroRelativeResidual[0], 0.0);
}

void testSharedRegistryRegistration(fvMesh& mesh)
{
    objectRegistry& registry =
        const_cast<objectRegistry&>(static_cast<const objectRegistry&>(mesh));
    HashSet<word> borrowedNames;

    volScalarField sharedP
    (
        IOobject
        (
            "sharedP",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("sharedP", dimPressure, 123.0),
        "zeroGradient"
    );

    checkTrue("shared registry initially lacks field", !registry.foundObject<volScalarField>("sharedP"));

    const bool inserted =
        poroCouplingRegistry::ensureVolScalarFieldRegistered
        (
            registry,
            sharedP,
            "fluidRegion",
            borrowedNames
        );

    checkTrue("shared registry inserts borrowed field", inserted);
    checkTrue("shared registry finds borrowed field", registry.foundObject<volScalarField>("sharedP"));
    checkTrue
    (
        "shared registry stores same field instance",
        &registry.lookupObject<volScalarField>("sharedP") == &sharedP
    );
    checkTrue("shared registry tracks borrowed field name", borrowedNames.found("sharedP"));

    const bool insertedAgain =
        poroCouplingRegistry::ensureVolScalarFieldRegistered
        (
            registry,
            sharedP,
            "fluidRegion",
            borrowedNames
        );

    checkTrue("shared registry registration is idempotent", !insertedAgain);
}

void runSharedRegistryConflict(fvMesh& mesh)
{
    objectRegistry& registry =
        const_cast<objectRegistry&>(static_cast<const objectRegistry&>(mesh));
    HashSet<word> borrowedNames;

    volScalarField existing
    (
        IOobject
        (
            "conflictingP",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("conflictingP", dimPressure, 1.0),
        "zeroGradient"
    );

    volScalarField requested
    (
        IOobject
        (
            "conflictingP",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("conflictingP", dimPressure, 2.0),
        "zeroGradient"
    );

    poroCouplingRegistry::ensureVolScalarFieldRegistered
    (
        registry,
        requested,
        "fluidRegion",
        borrowedNames
    );
}

void testPressureUnitScale(fvMesh& mesh)
{
    const dimensionedScalar pressureScale
    (
        poroPressureUnits::pressureScale(dimPressure, mesh)
    );

    checkNear("pressure unit scale is one for pressure fields", pressureScale.value(), 1.0);
    checkTrue("pressure unit scale is dimensionless for pressure fields", pressureScale.dimensions() == dimless);

    UniformDimensionedField<vector> gammaWater
    (
        IOobject
        (
            "gamma_water",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("gamma_water", dimPressure/dimLength, vector(0, 0, -9810.0))
    );

    const dimensionedScalar headScale
    (
        poroPressureUnits::pressureScale(dimLength, mesh)
    );

    checkNear("pressure unit scale uses |gamma_water| for head fields", headScale.value(), 9810.0);
    checkTrue
    (
        "pressure unit scale has gamma dimensions for head fields",
        headScale.dimensions() == dimPressure/dimLength
    );
}

void runHeadMissingGamma(fvMesh& mesh)
{
    poroPressureUnits::pressureScale(dimLength, mesh);
}

void testEffectiveStressModels(const fvMesh& mesh)
{
    volScalarField n
    (
        IOobject("nEff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("nEff", dimless, 0.35),
        "zeroGradient"
    );

    volScalarField S
    (
        IOobject("SEff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("SEff", dimless, 0.6),
        "zeroGradient"
    );

    volScalarField suction
    (
        IOobject("pEffSuction", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("pEffSuction", dimPressure, -100.0),
        "zeroGradient"
    );

    volScalarField pressure
    (
        IOobject("pEffPressure", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("pEffPressure", dimPressure, 100.0),
        "zeroGradient"
    );

    dictionary dict;

    autoPtr<effectiveStressModel> terzaghi =
        effectiveStressModel::New(dict, "terzaghi", mesh);
    autoPtr<effectiveStressModel> bishop =
        effectiveStressModel::New(dict, "bishop", mesh);
    autoPtr<effectiveStressModel> niemunis =
        effectiveStressModel::New(dict, "niemunis", mesh);
    autoPtr<effectiveStressModel> suctionCutOff =
        effectiveStressModel::New(dict, "suctionCutOff", mesh);

    const tmp<volScalarField> tTerzaghi(terzaghi->chi(n, S, suction));
    const tmp<volScalarField> tBishop(bishop->chi(n, S, suction));
    const tmp<volScalarField> tNiemunis(niemunis->chi(n, S, suction));
    const tmp<volScalarField> tSuctionOff(suctionCutOff->chi(n, S, suction));
    const scalar suctionOffChi = tSuctionOff()[0];
    const tmp<volScalarField> tPressureOn(suctionCutOff->chi(n, S, pressure));

    checkNear("effectiveStress terzaghi chi is one", tTerzaghi()[0], 1.0);
    checkNear("effectiveStress bishop chi equals saturation", tBishop()[0], 0.6);
    checkNear("effectiveStress niemunis chi golden", tNiemunis()[0], 0.86);
    checkNear("effectiveStress suctionCutOff ignores suction", suctionOffChi, 0.0);
    checkNear("effectiveStress suctionCutOff keeps positive pressure", tPressureOn()[0], 1.0);
    checkTrue("effectiveStress chi is dimensionless", tBishop().dimensions() == dimless);
    checkNear("effectiveStress bishop boundary chi equals saturation", tBishop().boundaryField()[0][0], 0.6);
}

void testVarSatMechanicalLawTerms(const fvMesh& mesh)
{
    const dimensionedScalar saturatedDensity("rhoSaturated", dimDensity, 2000.0);

    volScalarField S
    (
        IOobject("SForMechanicalTerms", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("SForMechanicalTerms", dimless, 0.6),
        "zeroGradient"
    );

    volScalarField n
    (
        IOobject("nForMechanicalTerms", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("nForMechanicalTerms", dimless, 0.25),
        "zeroGradient"
    );

    const dimensionedScalar rhoWater("rhoWater", dimDensity, 1000.0);
    const tmp<volScalarField> tRho
    (
        varSatPoroMechanicalLawTerms::mixtureDensity
        (
            saturatedDensity,
            S,
            n,
            rhoWater
        )
    );

    checkNear("varSatMechanical density golden", tRho()[0], 1900.0);
    checkTrue("varSatMechanical density dimensions", tRho().dimensions() == dimDensity);
    checkNear("varSatMechanical density boundary", tRho().boundaryField()[0][0], 1900.0);

    volSymmTensorField totalStress
    (
        IOobject("totalStressForMechanicalTerms", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor
        (
            "totalStressForMechanicalTerms",
            dimPressure,
            symmTensor(100.0, 10.0, 20.0, 200.0, 30.0, 300.0)
        ),
        "zeroGradient"
    );

    volScalarField chi
    (
        IOobject("chiForMechanicalTerms", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("chiForMechanicalTerms", dimless, 0.6),
        "zeroGradient"
    );

    volScalarField p
    (
        IOobject("pForMechanicalTerms", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("pForMechanicalTerms", dimPressure, 50.0),
        "zeroGradient"
    );

    const dimensionedScalar b("biotForMechanicalTerms", dimless, 0.8);
    const tmp<volSymmTensorField> tEffective
    (
        varSatPoroMechanicalLawTerms::initialEffectiveStress
        (
            totalStress,
            b,
            chi,
            p
        )
    );

    checkNear("varSatMechanical effective xx includes pore stress", tEffective()[0].xx(), 124.0);
    checkNear("varSatMechanical effective yy includes pore stress", tEffective()[0].yy(), 224.0);
    checkNear("varSatMechanical effective zz includes pore stress", tEffective()[0].zz(), 324.0);
    checkNear("varSatMechanical effective shear unchanged", tEffective()[0].xy(), 10.0);

    const tmp<volSymmTensorField> tTotal
    (
        varSatPoroMechanicalLawTerms::totalStress
        (
            tEffective(),
            b,
            chi,
            p
        )
    );

    checkNear("varSatMechanical total xx removes pore stress", tTotal()[0].xx(), totalStress[0].xx());
    checkNear("varSatMechanical total yy removes pore stress", tTotal()[0].yy(), totalStress[0].yy());
    checkNear("varSatMechanical total zz removes pore stress", tTotal()[0].zz(), totalStress[0].zz());
    checkNear("varSatMechanical total shear unchanged", tTotal()[0].xy(), totalStress[0].xy());
    checkTrue("varSatMechanical stress dimensions", tTotal().dimensions() == dimPressure);
}
}

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "registry-conflict",
        "Run the expected shared-registry duplicate-field failure path"
    );
    argList::addBoolOption
    (
        "head-missing-gamma",
        "Run the expected head-based pressure-unit failure path"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (args.found("head-missing-gamma"))
    {
        runHeadMissingGamma(mesh);
        FatalErrorInFunction
            << "head missing gamma test unexpectedly completed"
            << exit(FatalError);
    }

    if (args.found("registry-conflict"))
    {
        runSharedRegistryConflict(mesh);
        FatalErrorInFunction
            << "registry conflict test unexpectedly completed"
            << exit(FatalError);
    }

    testNDotRigidTranslation(mesh);
    testFixedStressStabil(mesh);
    testVarSatWeightingDecision(mesh);
    testCouplingFieldLifecycle(mesh);
    testRelativeAccelerationFlux(mesh);
    testExplicitCouplingSourceSign(mesh);
    testFvOptionCouplingRates(mesh);
    testFvOptionMatrixAssembly(mesh);
    testResidualOperations();
    testDeltaVfResiduals(mesh);
    testIterationControlDictionary(mesh);
    testLinearSolverResiduals(mesh);
    testMassBalanceTerms(mesh);
    testSharedRegistryRegistration(mesh);
    testPressureUnitScale(mesh);
    testEffectiveStressModels(mesh);
    testVarSatMechanicalLawTerms(mesh);

    if (failures)
    {
        FatalErrorInFunction
            << failures << " coupling term unit check(s) failed"
            << exit(FatalError);
    }

    Info<< "All coupling term unit checks passed" << nl;
    return 0;
}
