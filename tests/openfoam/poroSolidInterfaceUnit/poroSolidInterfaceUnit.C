#include "fvCFD.H"
#include "poroSolid.H"
#include "poroMechanicalLaw2.H"
#include "poroTractionFvPatchVectorField.H"

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

class TestPoroSolid
:
    public poroSolid
{
public:
    TestPoroSolid(Time& runTime)
    :
        poroSolid(runTime, "solid")
    {}

    void initializeHydraulicFields()
    {
        initializeSolidHydraulicFields();
    }

    void assembleTerms()
    {
        storeCouplingPressureReference();
        assembleCouplingTerms();
    }

    void clearTerms()
    {
        clearCouplingTerms();
    }
};

void testSaturatedSharedMesh(Time& runTime)
{
    TestPoroSolid coupling(runTime);

    checkTrue("poroSolid fixture uses shared mesh", coupling.sharedMesh());
    checkTrue("poroFluid creates default cellZone on shared mesh", coupling.poroFluidMesh().cellZones().size() == 1);
    checkTrue("default cellZone owns fixture cell", coupling.poroFluidMesh().cellZones().whichZone(0) == 0);
    checkNear("poroSolid fixture reads Biot coefficient", coupling.b()[0], 0.75);

    coupling.initializeHydraulicFields();

    checkTrue
    (
        "poroSolid registers p on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("p")
    );
    checkTrue
    (
        "poroSolid registers p_rgh on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("p_rgh")
    );
    checkTrue
    (
        "poroSolid shared p is fluid p instance",
        &coupling.solidMesh().lookupObject<volScalarField>("p")
     == &coupling.poroFluid().p()
    );
    checkTrue
    (
        "poroSolid shared p_rgh is fluid p_rgh instance",
        &coupling.solidMesh().lookupObject<volScalarField>("p_rgh")
     == &coupling.poroFluid().p_rgh()
    );

    volScalarField& p = coupling.poroFluidRef().p();
    volScalarField& pRgh = coupling.poroFluidRef().p_rgh();

    p.primitiveFieldRef() = 2000.0;
    p.boundaryFieldRef() = 2000.0;
    pRgh.primitiveFieldRef() = 400.0;
    pRgh.boundaryFieldRef() = 400.0;

    const label topPatchID =
        coupling.solidMesh().boundaryMesh().findPatchID("top");
    checkTrue("poroTraction fixture has top patch", topPatchID >= 0);

    fvPatchVectorField& topPatch =
        coupling.solidRef().D().boundaryFieldRef()[topPatchID];
    poroTractionFvPatchVectorField* poroTractionPatch =
        dynamic_cast<poroTractionFvPatchVectorField*>(&topPatch);
    checkTrue("top patch uses poroTraction", poroTractionPatch);

    const PtrList<mechanicalLaw>& laws = coupling.solid().mechanical();
    const poroMechanicalLaw2* poroLaw =
        dynamic_cast<const poroMechanicalLaw2*>(&laws[0]);
    checkTrue("fixture exposes poroMechanicalLaw2", poroLaw);

    if (poroTractionPatch && poroLaw)
    {
        const word& selectedPressureName = poroLaw->pressureFieldName();
        checkTrue
        (
            "fixture law selects a supported pressure field",
            selectedPressureName == "p" || selectedPressureName == "p_rgh"
        );

        const scalar selectedPressure =
            selectedPressureName == "p" ? 2000.0 : 400.0;
        const scalar expectedTotalPressure =
            poroLaw->biotCoeffValue().value()*selectedPressure;

        const scalar oldTimeValue = runTime.value();
        const label oldTimeIndex = runTime.timeIndex();
        runTime.setTime(1.0, oldTimeIndex + 1);

        poroTractionPatch->setUpdated(false);
        poroTractionPatch->updateCoeffs();

        checkNear
        (
            "poroTraction applies b times the law-selected pressure",
            poroTractionPatch->pressure()[0],
            expectedTotalPressure
        );
        checkNear
        (
            "poroTraction keeps zero prescribed effective traction",
            mag(poroTractionPatch->traction()[0]),
            0.0
        );

        runTime.setTime(oldTimeValue, oldTimeIndex);
    }

    coupling.assembleTerms();

    const tmp<volScalarField> tExplicit(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicit(coupling.implicitCouplingDtoP());
    const volScalarField impK(coupling.solid().mechanical().impK());

    checkTrue("poroSolid explicit coupling has volume source dimensions", tExplicit().dimensions() == dimless/dimTime);
    checkTrue("poroSolid implicit coupling has inverse pressure dimensions", tImplicit().dimensions() == inv(dimPressure));
    checkNear("poroSolid explicit coupling is zero for uniform D", tExplicit()[0], 0.0);
    checkNear("poroSolid implicit coupling uses scaled b squared over impK", tImplicit()[0], 1.5*sqr(coupling.b()[0])/impK[0]);

    coupling.clearTerms();
}
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    testSaturatedSharedMesh(runTime);

    if (failures)
    {
        FatalErrorInFunction
            << failures << " poroSolid interface unit check(s) failed"
            << exit(FatalError);
    }

    Info<< "All poroSolid interface unit checks passed" << nl;

    return 0;
}
