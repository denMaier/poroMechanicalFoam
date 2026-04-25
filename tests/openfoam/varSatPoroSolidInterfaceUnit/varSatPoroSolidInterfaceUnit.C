#include "fvCFD.H"
#include "fvc.H"
#include "fvOptionList.H"
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

    void afterFluid()
    {
        afterFluidSolve();
    }

    fvScalarMatrix fvOptionCouplingMatrix()
    {
        fvScalarMatrix eqn
        (
            poroFluidRef().pField(),
            dimVolume/dimTime
        );

        for
        (
            fv::option& source
          : static_cast<fv::optionList&>(poroFluidRef().fvOptions())
        )
        {
            source.addSup(eqn, 0);
        }

        return eqn;
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

    void seedAcceleratingDisplacementField(const scalar scale = 1.0)
    {
        volVectorField& D = solidRef().D();

        setFieldToUniform(D, vector::zero);
        D.oldTime();
        D.storeOldTime();

        forAll(D, cellI)
        {
            D[cellI] = scale*D.mesh().C()[cellI];
        }

        forAll(D.boundaryFieldRef(), patchI)
        {
            forAll(D.boundaryFieldRef()[patchI], faceI)
            {
                D.boundaryFieldRef()[patchI][faceI] =
                    scale*D.mesh().Cf().boundaryField()[patchI][faceI];
            }
        }

        D.correctBoundaryConditions();
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

void testVarSatFvOptionTransfersInterfaceTerms(TestVarSatPoroSolid& coupling)
{
    coupling.prepare();
    coupling.initializeHydraulicFields();

    volScalarField& fluidS = const_cast<volScalarField&>(coupling.poroFluid().S());
    setFieldToUniform(fluidS, 0.55);
    coupling.seedVelocityField(1.0);

    coupling.assembleTerms();

    const fvScalarMatrix optionMatrix(coupling.fvOptionCouplingMatrix());
    const tmp<volScalarField> tExplicit(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicit(coupling.implicitCouplingDtoP());
    const scalar deltaT = coupling.poroFluid().mesh().time().deltaTValue();
    const scalar expectedExplicitSource = 1.2375;
    const scalar expectedImplicitDiagonal = 4.6875e-8;

    checkNear
    (
        "fvOption diagonal is assembled implicit interface coupling over deltaT",
        optionMatrix.diag()[0],
        expectedImplicitDiagonal
    );
    checkNear
    (
        "fvOption source is assembled explicit interface coupling",
        optionMatrix.source()[0],
        expectedExplicitSource
    );
    checkNear
    (
        "interface implicit DTO matches independent fixture value",
        tImplicit()[0]/deltaT,
        expectedImplicitDiagonal
    );
    checkNear
    (
        "interface explicit DTO matches independent fixture value",
        tExplicit()[0],
        expectedExplicitSource
    );
    checkTrue
    (
        "fvOption transfer matrix keeps pressure equation dimensions",
        optionMatrix.dimensions() == dimVolume/dimTime
    );

    setFieldToUniform(fluidS, 0.25);
    coupling.assembleTerms();

    const fvScalarMatrix optionMatrixUpdated(coupling.fvOptionCouplingMatrix());
    const tmp<volScalarField> tExplicitUpdated(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicitUpdated(coupling.implicitCouplingDtoP());
    const scalar expectedUpdatedExplicitSource = 0.5625;

    checkNear
    (
        "fvOption transfer reflects reassembled explicit coupling",
        optionMatrixUpdated.source()[0],
        expectedUpdatedExplicitSource
    );
    checkNear
    (
        "reassembled explicit DTO matches independent fixture value",
        tExplicitUpdated()[0],
        expectedUpdatedExplicitSource
    );
    checkNotNear
    (
        "fvOption transfer source changes when saturation-weighted term changes",
        optionMatrixUpdated.source()[0],
        optionMatrix.source()[0]
    );
    checkNear
    (
        "fvOption transfer implicit diagonal remains saturation independent",
        optionMatrixUpdated.diag()[0],
        expectedImplicitDiagonal
    );

    coupling.clearTerms();
}

void testVarSatAccelerationTransferAndAfterFluidHook(TestVarSatPoroSolid& coupling)
{
    coupling.prepare();
    coupling.initializeHydraulicFields();

    volScalarField& fluidS = const_cast<volScalarField&>(coupling.poroFluid().S());
    setFieldToUniform(fluidS, 0.60);
    coupling.seedVelocityField(1.0);
    coupling.seedAcceleratingDisplacementField(1.0);

    const tmp<volScalarField> tDivU(fvc::div(coupling.solid().U()));
    const scalar nDotOnly =
        fluidS[0]*coupling.b()[0]*tDivU()[0];

    coupling.assembleTerms();

    const tmp<volScalarField> tExplicitWithAcceleration
    (
        coupling.explicitCouplingDtoP()
    );
    const tmp<volScalarField> tImplicitBeforeAfterFluid
    (
        coupling.implicitCouplingDtoP()
    );

    checkNotNear
    (
        "explicit DTO includes relative-acceleration source when D accelerates",
        tExplicitWithAcceleration()[0],
        nDotOnly
    );

    const fvScalarMatrix optionMatrix(coupling.fvOptionCouplingMatrix());
    const scalar deltaT = coupling.poroFluid().mesh().time().deltaTValue();
    const scalar expectedImplicitDiagonal = 4.6875e-8;

    checkNear
    (
        "fvOption source transfers acceleration-corrected explicit DTO",
        optionMatrix.source()[0],
        tExplicitWithAcceleration()[0]
    );
    checkNear
    (
        "fvOption diagonal transfers implicit DTO over deltaT with acceleration present",
        optionMatrix.diag()[0],
        expectedImplicitDiagonal
    );
    checkNear
    (
        "implicit DTO with acceleration present matches independent fixture value",
        tImplicitBeforeAfterFluid()[0]/deltaT,
        expectedImplicitDiagonal
    );

    coupling.afterFluid();

    const tmp<volScalarField> tExplicitAfterFluid(coupling.explicitCouplingDtoP());
    const tmp<volScalarField> tImplicitAfterFluid(coupling.implicitCouplingDtoP());

    checkNear
    (
        "afterFluidSolve clears transient acceleration and keeps nDot DTO",
        tExplicitAfterFluid()[0],
        nDotOnly
    );
    checkNear
    (
        "afterFluidSolve keeps implicit stabilization DTO available",
        tImplicitAfterFluid()[0],
        tImplicitBeforeAfterFluid()[0]
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
    testVarSatFvOptionTransfersInterfaceTerms(coupling);
    testVarSatAccelerationTransferAndAfterFluidHook(coupling);

    if (failures)
    {
        FatalErrorInFunction
            << failures << " variably saturated poroSolid interface unit check(s) failed"
            << exit(FatalError);
    }

    Info<< "All variably saturated poroSolid interface unit checks passed" << nl;

    return 0;
}
