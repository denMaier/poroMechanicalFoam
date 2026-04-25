#include "fvCFD.H"
#include "poroCouplingTerms.H"
#include "poroCouplingRegistry.H"
#include "poroPressureUnits.H"
#include "UniformDimensionedField.H"
#include "effectiveStressModel.H"

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
    testSharedRegistryRegistration(mesh);
    testPressureUnitScale(mesh);
    testEffectiveStressModels(mesh);

    if (failures)
    {
        FatalErrorInFunction
            << failures << " coupling term unit check(s) failed"
            << exit(FatalError);
    }

    Info<< "All coupling term unit checks passed" << nl;
    return 0;
}
