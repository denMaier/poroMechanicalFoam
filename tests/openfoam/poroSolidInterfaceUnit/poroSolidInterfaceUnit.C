#include "fvCFD.H"
#include "poroSolid.H"

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

    coupling.assembleTerms();

    const tmp<volScalarField> tExplicit(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicit(coupling.implicitCouplingDtoP());
    const volScalarField impK(coupling.solid().mechanical().impK());

    checkTrue("poroSolid explicit coupling has volume source dimensions", tExplicit().dimensions() == dimless/dimTime);
    checkTrue("poroSolid implicit coupling has inverse pressure dimensions", tImplicit().dimensions() == inv(dimPressure));
    checkNear("poroSolid explicit coupling is zero for uniform D", tExplicit()[0], 0.0);
    checkNear("poroSolid implicit coupling uses b squared over impK", tImplicit()[0], sqr(coupling.b()[0])/impK[0]);

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
