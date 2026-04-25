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

template<class Type>
void setFieldToUniform
(
    GeometricField<Type, fvPatchField, volMesh>& fld,
    const Type& value
)
{
    fld.primitiveFieldRef() = value;
    fld.boundaryFieldRef() = value;

    fld.correctBoundaryConditions();
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

    void seedVelocityField(const scalar scale = 1.0)
    {
        volVectorField& U = solidRef().U();

        forAll(U, cellI)
        {
            U[cellI] = scale*U.mesh().C()[cellI];
        }

        forAll(U.boundaryFieldRef(), patchI)
        {
            forAll(U.boundaryFieldRef()[patchI], faceI)
            {
                U.boundaryFieldRef()[patchI][faceI] =
                    scale*U.mesh().Cf().boundaryField()[patchI][faceI];
            }
        }

        U.correctBoundaryConditions();
    }
};

void testVarSatSharedMeshRegistration(TestVarSatPoroSolid& coupling)
{
    checkTrue("varSatPoroSolid fixture uses shared mesh", coupling.sharedMesh());
    checkTrue
    (
        "varSatPoroSolid selected variably saturated fluid",
        coupling.poroFluid().name() == "varSatPoroFluid"
    );
    checkTrue("default cellZone exists", coupling.poroFluidMesh().cellZones().size() == 1);
    checkTrue
    (
        "default cellZone owns fixture cell",
        coupling.poroFluidMesh().cellZones().whichZone(0) == 0
    );
    checkNear("Biot coefficient from poroMechanicalLaw2", coupling.b()[0], 0.75);

    const auto* varSatFluidPtr =
        dynamic_cast<const poroFluidModels::variablySaturatedPoroFluid*>
        (
            &coupling.poroFluid()
        );
    checkTrue("poroFluid exposes variablySaturated API", varSatFluidPtr);

    if (!varSatFluidPtr)
    {
        return;
    }

    const volScalarField& coupledS = coupling.poroFluid().S();
    const volScalarField& varSatFluidS = varSatFluidPtr->S();

    checkTrue
    (
        "poroFluid API S() delegates to variablySaturated model",
        &coupledS == &varSatFluidS
    );

    coupling.prepare();
    coupling.initializeHydraulicFields();

    checkTrue
    (
        "shared p field is registered on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("p")
    );
    checkTrue
    (
        "shared p_rgh field is registered on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("p_rgh")
    );
    checkTrue
    (
        "shared S field is registered on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("S")
    );
    checkTrue
    (
        "shared n field is registered on solid mesh",
        coupling.solidMesh().foundObject<volScalarField>("n")
    );

    checkTrue
    (
        "shared p points to the fluid p instance",
        &coupling.solidMesh().lookupObject<volScalarField>("p")
      == &coupling.poroFluid().p()
    );
    checkTrue
    (
        "shared p_rgh points to the fluid p_rgh instance",
        &coupling.solidMesh().lookupObject<volScalarField>("p_rgh")
      == &coupling.poroFluid().p_rgh()
    );
    checkTrue
    (
        "shared S points to the fluid S instance",
        &coupling.solidMesh().lookupObject<volScalarField>("S")
      == &coupling.poroFluid().S()
    );

    const volScalarField& firstSolidS =
        coupling.solidMesh().lookupObject<volScalarField>("S");

    coupling.initializeHydraulicFields();

    const volScalarField& secondSolidS =
        coupling.solidMesh().lookupObject<volScalarField>("S");

    checkTrue
    (
        "re-registering shared hydraulic fields preserves object identity",
        &firstSolidS == &secondSolidS
    );

    coupling.clearTerms();
}

void testVarSatSharedMeshCoupling(TestVarSatPoroSolid& coupling)
{
    coupling.prepare();
    coupling.initializeHydraulicFields();

    volScalarField& fluidS = const_cast<volScalarField&>(coupling.poroFluid().S());
    setFieldToUniform(fluidS, 0.65);

    coupling.seedVelocityField(1.0);
    const volVectorField& U = coupling.solid().U();
    const tmp<volScalarField> tDivU(fvc::div(U));
    const scalar divU0 = tDivU()[0];
    checkTrue("seeded velocity is dimensioned as velocity", U.dimensions() == dimVelocity);
    checkTrue("seeded velocity has non-zero divergence", mag(divU0) > SMALL);
    checkTrue("seeded velocity is returned through solid().U()", &U == &coupling.solid().U());

    coupling.assembleTerms();

    const tmp<volScalarField> tExplicit0(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicit0(coupling.implicitCouplingDtoP());
    const tmp<volScalarField> tImpK(coupling.solid().mechanical().impK());

    checkNear
    (
        "explicit coupling uses S*b weighted divergence",
        tExplicit0()[0],
        fluidS[0]*coupling.b()[0]*divU0
    );
    checkNotNear
    (
        "explicit coupling is not pure b weighted divergence",
        tExplicit0()[0],
        coupling.b()[0]*divU0
    );
    checkNear
    (
        "implicit coupling is b^2 / impK",
        tImplicit0()[0],
        sqr(coupling.b()[0])/tImpK()[0]
    );
    checkNotNear
    (
        "implicit coupling is not S*b^2 / impK",
        tImplicit0()[0],
        sqr(fluidS[0]*coupling.b()[0])/tImpK()[0]
    );
    checkTrue
    (
        "explicit coupling has nDot dimensions",
        tExplicit0().dimensions() == dimless/dimTime
    );
    checkTrue
    (
        "implicit coupling has inverse pressure dimensions",
        tImplicit0().dimensions() == inv(dimPressure)
    );

    setFieldToUniform(fluidS, 0.20);
    coupling.assembleTerms();

    const tmp<volScalarField> tExplicit1(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicit1(coupling.implicitCouplingDtoP());

    checkNear
    (
        "re-assembly updates explicit coupling with saturation",
        tExplicit1()[0],
        fluidS[0]*coupling.b()[0]*divU0
    );
    checkNear
    (
        "implicit coupling does not depend on saturation in this model",
        tImplicit1()[0],
        tImplicit0()[0]
    );
    checkNotNear
    (
        "re-assembly changed explicit coupling value",
        tExplicit1()[0],
        tExplicit0()[0]
    );

    coupling.clearTerms();
    setFieldToUniform(fluidS, 0.65);
    coupling.seedVelocityField(1.0);
    coupling.assembleTerms();

    const tmp<volScalarField> tExplicit2(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicit2(coupling.implicitCouplingDtoP());

    checkNear
    (
        "clearTerms() allows explicit terms to reassemble",
        tExplicit2()[0],
        fluidS[0]*coupling.b()[0]*divU0
    );
    checkNear
    (
        "clearTerms() allows implicit terms to reassemble",
        tImplicit2()[0],
        tImplicit0()[0]
    );

    coupling.clearTerms();
}
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    TestVarSatPoroSolid coupling(runTime);
    testVarSatSharedMeshRegistration(coupling);
    testVarSatSharedMeshCoupling(coupling);

    if (failures)
    {
        FatalErrorInFunction
            << failures << " variably saturated poroSolid interface unit check(s) failed"
            << exit(FatalError);
    }

    Info<< "All variably saturated poroSolid interface unit checks passed" << nl;

    return 0;
}
