# Momentum-diagonal checkerboard stabilization coefficient

## Objective

Construct the checkerboard-stabilization coefficient from an approximation of
the solid momentum equation diagonal:

```text
kStab ~ stabFactor/A_D
```

where `A_D` is the cell-centred diagonal returned by `fvVectorMatrix::A()`.
The reciprocal is calculated on cells and interpolated linearly to faces.

This only changes how the stabilization coefficient is constructed. It is
independent of whether the pressure-increment term is assembled implicitly or
explicitly in the pressure equation.

The pressure stabilization continues to act on the physical pore-pressure
rate. Its consistent form is compact minus wide:

```cpp
fvc::laplacian(kStabf, fvc::ddt(p_rgh))
- fvc::div
  (
      kStabf*mesh.Sf()
    & fvc::interpolate(fvc::grad(fvc::ddt(p_rgh)))
  )
```

It does not adopt solids4Foam's smoothing operator for the solid hydrostatic
pressure/Lagrange multiplier.

## First approach: explicitly construct an approximate momentum matrix

The first approach to test is to build a lightweight momentum matrix in
`poroSolidStab::assembleCouplingTerms()` solely to obtain its diagonal. The
matrix is not solved.

For `linGeomTotalDispSolid`, the actual momentum equation contains the
implicit stiffness Laplacian, inertia, optional damping, equation relaxation,
and additional constraints. A useful approximation of its left-hand side is:

```cpp
volVectorField& D = solidRef().D();
const volScalarField& rho = solid().rho();

const volScalarField impK
(
    solid().mechanical().impK()
);

const surfaceScalarField impKf
(
    fvc::interpolate(impK)
);

fvVectorMatrix momentumApprox
(
    rho*fvm::d2dt2(D)
  - fvm::laplacian(impKf, D, "laplacian(DD,D)")
);

if (solid().dampingCoeff().value() > SMALL)
{
    momentumApprox +=
        solid().dampingCoeff()*rho*fvm::ddt(D);
}

momentumApprox.relax();
```

Only `momentumApprox.A()` is then retained long enough to construct the
coefficient:

```cpp
const tmp<volScalarField> tA(momentumApprox.A());
const volScalarField& A = tA();

const dimensionedScalar minA
(
    "minA",
    A.dimensions(),
    VSMALL
);

const volScalarField rA
(
    IOobject
    (
        "checkerboardMomentumRDiag",
        runTime().timeName(),
        solidMesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    stabFactor_/max(A, minA)
);
```

The precise temporary-field construction may need small adjustments for the
OpenFOAM API, but the intended operation is `stabFactor/A()` on solid cells.

### Shared mesh

For a shared solid/fluid mesh, interpolate the cell reciprocal directly:

```cpp
checkerboardStabilCoeff_ = fvc::interpolate(rA);
```

The default `linear` interpolation scheme provides the requested simple face
interpolation.

### Separate meshes

For different solid and fluid meshes, take the reciprocal before mapping:

```cpp
const tmp<volScalarField> trAFluid
(
    solidToPoroFluid().mapTgtToSrc(rA)
);

checkerboardStabilCoeff_ = fvc::interpolate(trAFluid());
```

The required order is therefore:

```text
construct A_D on solid cells
    -> calculate 1/A_D on solid cells
    -> map 1/A_D to fluid cells
    -> linearly interpolate 1/A_D to fluid faces
```

Mapping `A_D` first and taking its reciprocal afterward is not equivalent and
should be avoided.

## Pressure-equation use

Once the face coefficient has been constructed, either existing temporal
treatment can use it.

### Explicit finite-increment form

```cpp
const volScalarField pDotPredictor
(
    (pCouplingRef() - pField().oldTime())/runTime().deltaT()
);

eqn += fvc::laplacian(checkerboardStabilCoeff_(), pDotPredictor);
eqn -= fvc::div
(
    checkerboardStabilCoeff_()
   *(
        mesh.Sf()
      & fvc::interpolate(fvc::grad(pDotPredictor))
    )
);
```

Both explicit stencils are evaluated from the frozen temporal pressure
increment that produced the incoming solid predictor.

### Implicit finite-increment form

```cpp
const tmp<surfaceScalarField> tImplicitCoeff
(
    checkerboardStabilCoeff_()/runTime().deltaT()
);

eqn += fvm::laplacian(tImplicitCoeff(), pField());
eqn -= fvc::laplacian(tImplicitCoeff(), pField().oldTime());

const volScalarField pDotPredictor
(
    (pCouplingRef() - pField().oldTime())/runTime().deltaT()
);
eqn -= fvc::div
(
    checkerboardStabilCoeff_()
   *(
        mesh.Sf()
      & fvc::interpolate(fvc::grad(pDotPredictor))
    )
);
```

The finite pressure and displacement changes are incremental storage amounts,
so they do not use `fvSchemes/ddtSchemes`. Only true tangent-storage terms such
as `fvm::ddt(Ss, p)` need a selectable time scheme. The expressions above put
compact-minus-wide into the coupling `fvOption`. Since that option is on the
right-hand side of `lhs == fvOptions(...)`, OpenFOAM subtracts it during final
assembly and the pressure matrix receives the stabilizing wide-minus-compact
contribution.

## Expected scaling

For the total-displacement momentum equation, the approximate diagonal scales
roughly as:

```text
A_D ~ K/h^2 + rho/deltaT^2 + damping/deltaT.
```

Thus:

```text
1/A_D -> h^2/K
```

in the quasi-static stiffness-dominated limit. When inertia dominates,
`1/A_D` also introduces the desired time-step dependence naturally.

## Difference from retaining the actual solids4Foam diagonal

The explicitly constructed matrix avoids modifying solids4Foam and does not
depend on the lifetime of its local `DEqn` object. It can also be evaluated
before the fluid solve in every coupling iteration.

It is only an approximation of the actual solved momentum matrix. Depending
on the selected solid model, it may omit:

- solid-model-specific implicit terms;
- `fvOptions` contributions;
- imposed-cell constraints applied through `setCellDisps`;
- details of a nonlinear or block-coupled Jacobian.

Equation relaxation is included when `momentumApprox.relax()` is called.

A later alternative is to catch the temporary `DEqnA` field inside
`poroMechanicalLaw2::correct(volSymmTensorField&)` and retain it as a
persistent base-solid-mesh field. That alternative captures the actual
assembled diagonal, but it is only available after a solid solve. Since the
coupling loop solves the fluid first, the first fluid solve would still need
the explicitly constructed coefficient as a fallback.

## Suggested configuration

Coefficient construction and pressure-equation temporal treatment should be
configured independently:

```foam
stabilizationType          explicit;          // explicit or implicit pressure term
stabilizationCoefficient   momentumDiagonal; // stiffness or momentumDiagonal
stabFactor                 1.0;
```

Here, `momentumDiagonal` initially means the explicitly assembled approximate
momentum matrix described above. A later option such as
`retainedMomentumDiagonal` can select the captured solids4Foam diagonal if
needed.
