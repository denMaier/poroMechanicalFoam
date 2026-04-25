#include "fvCFD.H"
#include "fvc.H"
#include "varSatPoroSolid.H"
#include "variablySaturatedPoroFluid.H"

using namespace Foam;
using namespace Foam::poroSolidInteractions;

namespace
{
static label failures = 0;

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

void checkNear
(
    const word& name,
    const scalar actual,
    const scalar expected,
    const scalar tolerance = 1e-11
)
{
    const scalar error = mag(actual - expected)/max(scalar(1), mag(expected));
    if (error > tolerance)
    {
        ++failures;
        Info<< "FAIL " << name
            << ": actual=" << actual
            << " expected=" << expected
            << " relErr=" << error
            << " tolerance=" << tolerance << nl;
    }
    else
    {
        Info<< "PASS " << name << nl;
    }
}

void checkNotNear
(
    const word& name,
    const scalar actual,
    const scalar unexpected,
    const scalar tolerance = 1e-11
)
{
    const scalar error = mag(actual - unexpected)/max(scalar(1), mag(unexpected));
    if (error <= tolerance)
    {
        ++failures;
        Info<< "FAIL " << name
            << ": actual=" << actual
            << " unexpected=" << unexpected
            << " relErr=" << error
            << " tolerance=" << tolerance << nl;
    }
    else
    {
        Info<< "PASS " << name << nl;
    }
}

class TestVarSatPoroSolid
:
    public varSatPoroSolid
{
public:
    TestVarSatPoroSolid(Time& runTime)
    :
        varSatPoroSolid(runTime, "solid")
    {}

    void prepare()
    {
        prepareCouplingLoop();
    }

    void initializeHydraulicFields()
    {
        initializeSolidHydraulicFields();
    }

    void assembleTerms()
    {
        assembleCouplingTerms();
    }

    void clearTerms()
    {
        clearCouplingTerms();
    }
};

void setUniform(volScalarField& field, const scalar value)
{
    field.internalFieldRef() = value;

    forAll(field.boundaryFieldRef(), patchI)
    {
        field.boundaryFieldRef()[patchI] = value;
    }
}

void imposePositiveDivergence(volVectorField& U)
{
    U.internalFieldRef() = vector::zero;

    forAll(U.boundaryFieldRef(), patchI)
    {
        vectorField& patchU = U.boundaryFieldRef()[patchI];
        const vectorField& patchSf = U.mesh().Sf().boundaryField()[patchI];

        forAll(patchU, faceI)
        {
            patchU[faceI] =
                patchSf[faceI].x() > 0
              ? vector(1, 0, 0)
              : vector::zero;
        }
    }
}

void testVarSatSharedMesh(Time& runTime)
{
    TestVarSatPoroSolid coupling(runTime);
    Info<< "poroFluid model has foundObject(p_rgh): "
        << coupling.poroFluid().foundObject<volScalarField>("p_rgh") << endl;
    Info<< "poroFluid mesh has S: "
        << coupling.poroFluidMesh().foundObject<volScalarField>("S") << endl;
    Info<< "poroFluid mesh has p: "
        << coupling.poroFluidMesh().foundObject<volScalarField>("p") << endl;
    Info<< "poroFluid mesh has p_rgh: "
        << coupling.poroFluidMesh().foundObject<volScalarField>("p_rgh") << endl;
    Info<< "poroFluid object names before explicit S() call: "
        << coupling.poroFluid().p().db().sortedNames() << endl;

    checkTrue("varSatPoroSolid fixture uses shared mesh", coupling.sharedMesh());
    checkTrue("varSatPoroSolid selected variably saturated fluid", coupling.poroFluid().name() == "varSatPoroFluid");
    checkTrue("varSatPoroSolid default cellZone owns fixture cell", coupling.poroFluidMesh().cellZones().whichZone(0) == 0);
    Info<< "CHECKPOINT: before varSat fluid cast" << endl;
    const poroFluidModels::variablySaturatedPoroFluid* varSatFluidPtr =
        dynamic_cast<const poroFluidModels::variablySaturatedPoroFluid*>
        (
            &coupling.poroFluid()
        );
    Info<< "CHECKPOINT: after varSat fluid cast" << endl;
    checkTrue("varSatPoroSolid exposes variably saturated fluid API", varSatFluidPtr);
    Info<< "poroFluid db name: " << coupling.poroFluid().p().db().name() << endl;
    Info<< "poroFluid db object names: "
        << coupling.poroFluid().p().db().sortedNames() << endl;
    Info<< "poroFluid DB direct foundObject<S> before API: "
        << coupling.poroFluid().foundObject<volScalarField>("S") << endl;
    Info<< "CHECKPOINT: after precondition checks" << endl;
    if (!varSatFluidPtr)
    {
        checkTrue("varSatPoroSolid coupling aborted due unsupported fluid type", false);
        return;
    }

    const volScalarField& coupledS = coupling.poroFluid().S();
    const volScalarField* saturationPointerBefore = &coupledS;
    checkTrue
    (
        "poroFluid API S() matches derived variablySaturatedS()",
        &coupledS == &varSatFluidPtr->S()
    );

    volScalarField& S = const_cast<volScalarField&>(coupledS);
    setUniform(S, 0.6);

    volVectorField& U = coupling.solidRef().U();
    imposePositiveDivergence(U);

    coupling.prepare();
    coupling.initializeHydraulicFields();
    checkTrue
    (
        "varSatPoroSolid keeps a single registered S field on fluid side",
        saturationPointerBefore == &coupling.poroFluid().S()
    );

    checkTrue
    (
        "varSatPoroSolid registers S on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("S")
    );
    const bool solidHasS = coupling.solidMesh().foundObject<volScalarField>("S");
    if (solidHasS)
    {
        checkTrue
        (
            "varSatPoroSolid shared S is fluid S instance",
            &coupling.solidMesh().lookupObject<volScalarField>("S") == &S
        );
    }
    checkTrue
    (
        "varSatPoroSolid registers n on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("n")
    );

    const tmp<volScalarField> tDivU(fvc::div(U));
    checkTrue("varSatPoroSolid fixture has nonzero velocity divergence", mag(tDivU()[0]) > SMALL);

    coupling.assembleTerms();

    const tmp<volScalarField> tExplicit(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicit(coupling.implicitCouplingDtoP());
    const volScalarField impK(coupling.solid().mechanical().impK());

    const scalar weightedExpected = S[0]*coupling.b()[0]*tDivU()[0];
    const scalar unweightedExpected = coupling.b()[0]*tDivU()[0];
    const scalar stabilExpected = sqr(coupling.b()[0])/impK[0];
    const scalar saturationScaledStabil = sqr(S[0]*coupling.b()[0])/impK[0];

    checkNear("varSat explicit coupling uses S*b weighting", tExplicit()[0], weightedExpected);
    checkNotNear("varSat explicit coupling is not unweighted b", tExplicit()[0], unweightedExpected);
    checkNear("varSat implicit coupling keeps b based stabilization", tImplicit()[0], stabilExpected);
    checkNotNear("varSat implicit coupling is not S*b stabilization", tImplicit()[0], saturationScaledStabil);

    coupling.clearTerms();
}
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    testVarSatSharedMesh(runTime);

    if (failures)
    {
        FatalErrorInFunction
            << failures << " variably saturated poroSolid interface unit check(s) failed"
            << exit(FatalError);
    }

    Info<< "All variably saturated poroSolid interface unit checks passed" << nl;

    return 0;
}
