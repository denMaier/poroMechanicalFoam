#include "fvCFD.H"
#include "saturationLaw.H"
#include "conductivityModel.H"
#include "storageLaw.H"
#include "poroHydraulicModel.H"

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

void checkBetween
(
    const word& name,
    const scalar actual,
    const scalar lower,
    const scalar upper,
    const scalar tolerance = 1e-12
)
{
    checkTrue
    (
        name,
        actual >= lower - tolerance && actual <= upper + tolerance
    );
}

dimensionedScalar dimlessScalar(const word& name, const scalar value)
{
    return dimensionedScalar(name, dimless, value);
}

dimensionedScalar pressureScalar
(
    const word& name,
    const dimensionSet& pressureDims,
    const scalar value
)
{
    return dimensionedScalar(name, pressureDims, value);
}

autoPtr<saturationLaw> makeVanGenuchten(volScalarField& p)
{
    dictionary properties;
    dictionary coeffs;
    coeffs.add("n", dimlessScalar("n", 1.5));
    coeffs.add("S_r", dimlessScalar("S_r", 0.1));
    coeffs.add("S_0", dimlessScalar("S_0", 0.95));
    coeffs.add("alpha", pressureScalar("alpha", dimless/p.dimensions(), 1e-4));
    properties.add("vanGenuchtenCoeffs", coeffs);
    return saturationLaw::New("vanGenuchten", properties, p);
}

autoPtr<saturationLaw> makeBrooksCorey(volScalarField& p)
{
    dictionary properties;
    dictionary coeffs;
    coeffs.add("n", dimlessScalar("n", 2.0));
    coeffs.add("S_r", dimlessScalar("S_r", 0.12));
    coeffs.add("S_pe", dimlessScalar("S_pe", 0.92));
    coeffs.add("p_e", pressureScalar("p_e", p.dimensions(), -100.0));
    properties.add("brooksCoreyCoeffs", coeffs);
    return saturationLaw::New("brooksCorey", properties, p);
}

autoPtr<conductivityModel> makeKozenyCarman(volScalarField& p)
{
    dictionary properties;
    dictionary coeffs;
    coeffs.add("D50", dimensionedScalar("D50", dimLength, 2e-4));
    coeffs.add("viscosity", dimensionedScalar("viscosity", dimless, 1.3e-6));
    properties.add("kozenyCarmanCoeffs", coeffs);
    return conductivityModel::New("kozenyCarman", properties, p);
}

autoPtr<storageLaw> makeStorageCoeff(volScalarField& p)
{
    dictionary properties;
    dictionary coeffs;
    coeffs.add("Ss", dimensionedScalar("Ss", dimless/p.dimensions(), 2.5e-8));
    properties.add("storageCoeffCoeffs", coeffs);
    return storageLaw::New("storageCoeff", properties, p);
}

autoPtr<storageLaw> makeKPrime(volScalarField& p, const bool pressureDependent)
{
    dictionary properties;
    properties.add("Cs", dimensionedScalar("Cs", dimless/p.dimensions(), 1e-8));

    dictionary coeffs;
    coeffs.add("Sw_0", dimlessScalar("Sw_0", 0.8));
    coeffs.add("Kw", pressureScalar("Kw", p.dimensions(), 2.2e9));
    coeffs.add("p_At", pressureScalar("p_At", p.dimensions(), 101325.0));
    coeffs.add("pDependent", pressureDependent);
    properties.add("KPrimeCoeffs", coeffs);
    return storageLaw::New("KPrime", properties, p);
}

autoPtr<storageLaw> makeSkemptonB(volScalarField& p)
{
    dictionary properties;
    properties.add("Cs", dimensionedScalar("Cs", dimless/p.dimensions(), 1e-8));

    dictionary coeffs;
    coeffs.add("BParameter", dimlessScalar("BParameter", 0.63));
    coeffs.add("Km", pressureScalar("Km", p.dimensions(), 7.5e7));
    properties.add("skemptonBCoeffs", coeffs);
    return storageLaw::New("skemptonB", properties, p);
}

autoPtr<storageLaw> makeMontenegroStorage(volScalarField& p)
{
    dictionary properties;
    properties.add("Cs", dimensionedScalar("Cs", dimless/p.dimensions(), 1e-8));

    dictionary coeffs;
    coeffs.add("S_e", dimlessScalar("S_e", 0.82));
    coeffs.add("p_At", pressureScalar("p_At", p.dimensions(), 101325.0));
    coeffs.add("p_e", pressureScalar("p_e", p.dimensions(), -10000.0));
    properties.add("montenegroCoeffs", coeffs);
    return storageLaw::New("montenegro", properties, p);
}

void testVanGenuchten(saturationLaw& law)
{
    const scalar Sr = 0.1;
    const scalar S0 = 0.95;

    checkNear("vanGenuchten.pStar golden", law.pStar(0), -4807.49856769136);
    checkNear("vanGenuchten.S(-2500) golden", law.S(-2500.0, 0), 0.917274756507531);
    checkNear("vanGenuchten.kr(-2500) golden", law.kr(-2500.0, 0), 0.264379532079508);
    checkNear("vanGenuchten.C(-2500) golden", law.C(-2500.0, 0), 1.81616612557229e-05);
    checkNear("vanGenuchten.S(-100) golden", law.S(-100.0, 0), 0.949716855408764);
    checkNear("vanGenuchten.kr(-100) golden", law.kr(-100.0, 0), 0.809925029846204);
    checkNear("vanGenuchten.C(-100) golden", law.C(-100.0, 0), 4.24433993710635e-06);

    checkNear("vanGenuchten.S(0) is S_0", law.S(0.0, 0), S0);
    checkNear("vanGenuchten.kr(0) saturated", law.kr(0.0, 0), 1.0);
    checkNear("vanGenuchten.C(0) saturated", law.C(0.0, 0), 0.0);
    checkNear("vanGenuchten.S(+100) is S_0", law.S(100.0, 0), S0);
    checkNear("vanGenuchten.kr(+100) saturated", law.kr(100.0, 0), 1.0);

    checkBetween("vanGenuchten.S(-2500) within physical bounds", law.S(-2500.0, 0), Sr, S0);
    checkBetween("vanGenuchten.S(-100) within physical bounds", law.S(-100.0, 0), Sr, S0);
    checkBetween("vanGenuchten.kr(-2500) within [0,1]", law.kr(-2500.0, 0), 0.0, 1.0);
    checkBetween("vanGenuchten.kr(-100) within [0,1]", law.kr(-100.0, 0), 0.0, 1.0);
    checkTrue("vanGenuchten.S increases with pressure", law.S(-2500.0, 0) < law.S(-100.0, 0));
    checkTrue("vanGenuchten.kr increases with pressure", law.kr(-2500.0, 0) < law.kr(-100.0, 0));
    checkTrue("vanGenuchten.C non-negative", law.C(-2500.0, 0) >= 0.0 && law.C(-100.0, 0) >= 0.0);
}

void testBrooksCorey(saturationLaw& law)
{
    const scalar Sr = 0.12;
    const scalar S0 = 0.92;
    const scalar pD = -100.0;

    checkNear("brooksCorey.pStar is air-entry pressure", law.pStar(0), pD);
    checkNear("brooksCorey.S(-250) golden", law.S(-250.0, 0), 0.248);
    checkNear("brooksCorey.kr(-250) golden", law.kr(-250.0, 0), 0.00065536);
    checkNear("brooksCorey.C(-250) golden", law.C(-250.0, 0), 0.001024);
    checkNear("brooksCorey.S(p_e) is S_pe", law.S(-100.0, 0), S0);
    checkNear("brooksCorey.kr(p_e) saturated", law.kr(-100.0, 0), 1.0);
    checkNear("brooksCorey.C(p_e) golden", law.C(-100.0, 0), 0.016);
    checkNear("brooksCorey.S above p_e is S_pe", law.S(-50.0, 0), S0);
    checkNear("brooksCorey.kr above p_e saturated", law.kr(-50.0, 0), 1.0);
    checkNear("brooksCorey.C above p_e zero", law.C(-50.0, 0), 0.0);

    checkBetween("brooksCorey.S(-250) within physical bounds", law.S(-250.0, 0), Sr, S0);
    checkBetween("brooksCorey.kr(-250) within [0,1]", law.kr(-250.0, 0), 0.0, 1.0);
    checkTrue("brooksCorey.S increases to air-entry pressure", law.S(-250.0, 0) < law.S(-100.0, 0));
    checkTrue("brooksCorey.kr increases to air-entry pressure", law.kr(-250.0, 0) < law.kr(-100.0, 0));
    checkTrue("brooksCorey.C non-negative", law.C(-250.0, 0) >= 0.0 && law.C(-100.0, 0) >= 0.0);
}

void testKozenyCarman(conductivityModel& model)
{
    checkNear("kozenyCarman.k golden", model.k(0), 2.18156525034832e-05);
    checkTrue("kozenyCarman.k positive", model.k(0) > 0.0);
}

void testStorageCoeff(storageLaw& law)
{
    checkNear("storageCoeff.Ss golden", law.Ss(0.37, -25000.0, 0), 2.5e-8);
    checkTrue("storageCoeff preserves runtime name", law.name() == "storageCoeff");
    checkTrue("storageCoeff is static", !law.updatesSs());
}

void testKPrime(storageLaw& staticLaw, storageLaw& pressureDependentLaw)
{
    checkNear
    (
        "KPrime pressure-independent Ss golden",
        staticLaw.Ss(0.37, -25000.0, 0),
        7.36757762824395e-07
    );
    checkNear
    (
        "KPrime pressure-dependent Ss golden",
        pressureDependentLaw.Ss(0.37, -25000.0, 0),
        9.75972704642230e-07
    );
    checkTrue
    (
        "KPrime pressure-dependent storage increases as p_At + p decreases",
        pressureDependentLaw.Ss(0.37, -25000.0, 0)
      > pressureDependentLaw.Ss(0.37, 0.0, 0)
    );
    checkTrue("KPrime updates storage", staticLaw.updatesSs());
}

void testSkemptonB(storageLaw& law)
{
    checkNear("skemptonB.Ss golden", law.Ss(0.37, -25000.0, 0), 1.41306878306878e-08);
    checkTrue("skemptonB is static", !law.updatesSs());
    checkTrue("skemptonB.Ss positive", law.Ss(0.37, -25000.0, 0) > 0.0);
}

void testMontenegroStorage(storageLaw& law)
{
    checkNear("montenegroStorage.Ss above p_e golden", law.Ss(0.37, -5000.0, 0), 7.33598691485139e-07);
    checkNear("montenegroStorage.Ss at p_e cut-off", law.Ss(0.37, -10000.0, 0), 0.0);
    checkNear("montenegroStorage.Ss below p_e cut-off", law.Ss(0.37, -15000.0, 0), 0.0);
    checkTrue("montenegroStorage updates storage", law.updatesSs());
}

void testPoroHydraulicModelZoneAssembly(fvMesh& mesh, volScalarField& p)
{
    checkTrue("poroHydraulicModel fixture has two hydraulic zones", mesh.cellZones().size() == 2);
    checkTrue("poroHydraulicModel fixture has zoneA", mesh.cellZones().findZoneID("zoneA") >= 0);
    checkTrue("poroHydraulicModel fixture has zoneB", mesh.cellZones().findZoneID("zoneB") >= 0);

    poroHydraulicModel model
    (
        p,
        dimensionedVector("g", dimAcceleration, vector(0, -9.81, 0))
    );

    const tmp<volScalarField> tn0(model.n0());
    checkNear("poroHydraulicModel zoneA n0", tn0()[0], 0.25);
    checkNear("poroHydraulicModel zoneB n0", tn0()[1], 0.45);

    const tmp<volScalarField> tSs(model.Ss(tn0(), p));
    checkNear("poroHydraulicModel zoneA Ss", tSs()[0], 1e-8);
    checkNear("poroHydraulicModel zoneB Ss", tSs()[1], 4e-8);

    const tmp<volScalarField> tk(model.k());
    checkNear("poroHydraulicModel zoneA k", tk()[0], 2e-6);
    checkNear("poroHydraulicModel zoneB k", tk()[1], 8e-6);
    checkTrue("poroHydraulicModel zone constants do not update Ss", !model.updatesSs());
    checkTrue("poroHydraulicModel zone constants do not update k", !model.updatesK());
}
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("p", dimPressure, -100.0)
    );

    volScalarField n
    (
        IOobject
        (
            "n",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("n", dimless, 0.37)
    );

    autoPtr<saturationLaw> vanGenuchten = makeVanGenuchten(p);
    testVanGenuchten(vanGenuchten());

    autoPtr<saturationLaw> brooksCorey = makeBrooksCorey(p);
    testBrooksCorey(brooksCorey());

    autoPtr<conductivityModel> kozenyCarman = makeKozenyCarman(p);
    testKozenyCarman(kozenyCarman());

    autoPtr<storageLaw> storageCoeff = makeStorageCoeff(p);
    testStorageCoeff(storageCoeff());

    autoPtr<storageLaw> kPrimeStatic = makeKPrime(p, false);
    autoPtr<storageLaw> kPrimePressureDependent = makeKPrime(p, true);
    testKPrime(kPrimeStatic(), kPrimePressureDependent());

    autoPtr<storageLaw> skemptonB = makeSkemptonB(p);
    testSkemptonB(skemptonB());

    autoPtr<storageLaw> montenegroStorage = makeMontenegroStorage(p);
    testMontenegroStorage(montenegroStorage());

    testPoroHydraulicModelZoneAssembly(mesh, p);

    if (failures)
    {
        FatalErrorInFunction
            << failures << " hydraulic model unit check(s) failed"
            << exit(FatalError);
    }

    Info<< "All hydraulic model unit checks passed" << nl;
    return 0;
}
