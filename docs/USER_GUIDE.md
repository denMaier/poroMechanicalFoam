# `poroMechanicalFoam` User Guide

## What This Repository Is For

`poroMechanicalFoam` is a porous-media extension for OpenFOAM and solids4Foam. It is intended for:

- groundwater flow in porous media
- coupled pore-pressure and deformation problems
- variably saturated problems using soil-water characteristic curves
- geotechnical workflows that need porous hydraulic laws and porous mechanical laws together

## Supported Versions

- OpenFOAM `v2412` or `v2512`
- solids4Foam `v2.3` (handled automatically by `Allwmake`)

## Installation

### Build with the bundled solids4Foam submodule

1. Source OpenFOAM.
2. Enter the repository root.
3. Run:

```bash
./Allwmake
```

On first build, `Allwmake` initializes `external/solids4foam` and sparse-checks out only the parts this repository needs.

### Build against an existing solids4Foam tree

```bash
export S4F_ROOT=/path/to/solids4foam
./Allwmake
```

### Optional Abaqus UMAT build

```bash
export S4F_USE_GFORTRAN=1
./Allwmake
```

## What Gets Built

The default build produces two executables placed in `$(FOAM_USER_APPBIN)`:

- `poroMechanicalFoam` — main solver
- `initPoroMechanicalFoam` — initialization utility

## Choosing a Workflow

`poroMechanicalFoam` selects its active physics at runtime from the case dictionaries, so the same executable handles all cases:

| Goal | `physicsModel` in `system/controlDict` |
|---|---|
| Pure solid mechanics | `solid` (solids4Foam built-in) |
| Saturated groundwater flow only | `poroFluid` |
| Saturated coupled poromechanics | `poroSolid` |
| Variably saturated coupled poromechanics | `poroSolid` (with `varSatPoroFluid` inside) |

Set the entry in `system/controlDict`:

```foam
application     poroMechanicalFoam;
```

## Required Case Dictionaries

### Saturated flow-only case

| File | Purpose |
|---|---|
| `system/controlDict` | application, time control |
| `system/fvSchemes` | discretization schemes |
| `system/fvSolution` | solvers and tolerances |
| `constant/g` | gravity vector |
| `constant/poroFluidProperties` | selects the fluid model type |
| `constant/poroHydraulicProperties` | hydraulic constitutive laws |

### Coupled poromechanics case (adds to the above)

| File | Purpose |
|---|---|
| `constant/poroCouplingProperties` | coupling interface and options |
| `constant/mechanicalProperties` | solid constitutive law |

---

## `constant/poroFluidProperties`

This file selects the pore-fluid solver. The `poroFluidModel` entry is the runtime selector:

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      poroFluidProperties;
}
// ----------------------------- //

poroFluidModel      poroFluid;          // saturated flow
// poroFluidModel   varSatPoroFluid;    // variably saturated, pressure form
// poroFluidModel   varSatPoroFluidHead; // variably saturated, head form

poroFluidCoeffs
{
    enabled         true;
    iterations      20;
    infoFrequency   1;

    convergence
    {
        p_rgh       max 1e-8;
    }
}
```

Available types and their corresponding coeffs sub-dict name:

| Type | Governing equation | Primary field | Coeffs sub-dict |
|---|---|---|---|
| `poroFluid` | saturated Darcy flow | `p_rgh` | `poroFluidCoeffs` |
| `varSatPoroFluid` | Richards equation, pressure form | `p_rgh` | `varSatPoroFluidCoeffs` |
| `varSatPoroFluidHead` | Richards equation, head form | `pHead` | `varSatPoroFluidHeadCoeffs` |

The `solutionAlgorithm` entry inside the coeffs sub-dict applies to variably saturated solvers and selects the Richards linearization: `Standard`, `Celia`, `LScheme`, or `Casulli`.

---

## `constant/poroHydraulicProperties`

This file provides the hydraulic constitutive laws. The model is zone-based: each cell zone in the mesh gets its own set of laws, read from a sub-dictionary named after the zone.

**Default-zone shortcut:** if the mesh has no cell zones, the solver automatically creates one named `defaultZone` covering all cells. When looking up the zone's sub-dict, the code falls back to the whole dictionary if no matching sub-dict is found — so a flat file with all entries at the top level works as the default for that zone. This is the recommended approach when you do not need zone-specific material properties.

For multi-zone meshes, add one named sub-dict per zone. Entries inside the zone sub-dict override the top-level defaults.

### Saturated case — single-material (flat/default-zone) layout

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      poroHydraulicProperties;
}
// ----------------------------- //

rho     rho [1 -3 0 0 0 0 0] 1000;   // fluid density
href    href [0 1 0 0 0 0 0] 0;      // reference water-table elevation

n               n [0 0 0 0 0 0 0] 0.35;
storageLaw      storageCoeff;
conductivityLaw constant;

storageCoeffCoeffs
{
    Ss      Ss [-1 1 2 0 0 0 0] 1e-5;  // specific storage [1/Pa]
}

k       k [0 1 -1 0 0 0 0] 1e-5;      // hydraulic conductivity K [m/s]
```

### Variably saturated case — van Genuchten, flat layout

```foam
rho     rho [1 -3 0 0 0 0 0] 1000;
href    href [0 1 0 0 0 0 0] 0;

n               n [0 0 0 0 0 0 0] 0.35;
storageLaw      storageCoeff;
conductivityLaw constant;
saturationLaw   vanGenuchten;

storageCoeffCoeffs
{
    Ss      Ss [-1 1 2 0 0 0 0] 1e-5;  // [1/Pa] for pressure form
}

k       k [0 1 -1 0 0 0 0] 1e-5;

vanGenuchtenCoeffs
{
    // alpha must satisfy [1/alpha] == [p]. For Pa-based pressure form:
    alpha   alpha [-1 1 2 0 0 0 0] 1e-4;  // [1/Pa]
    n       n [0 0 0 0 0 0 0]      1.5;
    S_r     S_r [0 0 0 0 0 0 0]    0.1;
    S_0     S_0 [0 0 0 0 0 0 0]    1.0;
}
```

### Multi-zone layout

When the mesh has named cell zones, add one sub-dict per zone. Each sub-dict provides the same entries as the flat layout above:

```foam
rho     rho [1 -3 0 0 0 0 0] 1000;
href    href [0 1 0 0 0 0 0] 0;

sand
{
    n               n [0 0 0 0 0 0 0] 0.38;
    storageLaw      storageCoeff;
    conductivityLaw constant;
    storageCoeffCoeffs { Ss Ss [-1 1 2 0 0 0 0] 2e-5; }
    k   k [0 1 -1 0 0 0 0] 1e-4;
}

clay
{
    n               n [0 0 0 0 0 0 0] 0.45;
    storageLaw      storageCoeff;
    conductivityLaw constant;
    storageCoeffCoeffs { Ss Ss [-1 1 2 0 0 0 0] 5e-5; }
    k   k [0 1 -1 0 0 0 0] 1e-8;
}
```

> **alpha and Ss dimensions by formulation:**
> - Pressure form (`p_rgh` in Pa): `[-1 1 2 0 0 0 0]` (= 1/Pa)
> - Head form (`pHead` in m): `[0 -1 0 0 0 0 0]` (= 1/m)
>
> The solver checks that `[1/alpha] == [p_rgh]` and aborts with a dimension error on mismatch.

### Available hydraulic laws

**Storage laws** (`storageLaw`):

| Type | Coeffs sub-dict | Key entries |
|---|---|---|
| `storageCoeff` | `storageCoeffCoeffs` | `Ss` |
| `skemptonB` | `skemptonBCoeffs` | `B`, `n0`, `Kf`, `Ks` |
| `KPrime` | `KPrimeCoeffs` | `KPrime` |
| `montana` | — | — |
| `pressureFunction` | — | tabulated |

**Conductivity laws** (`conductivityLaw`):

| Type | Key entries |
|---|---|
| `constant` | `k` — hydraulic conductivity K [m/s] |
| `anisotropicConstant` | `K` (tensor) |
| `kozenyCarman` | `k0`, `n0` |
| `porosityFunction` | tabulated |
| `gradientFunction` / `limitedGradient` | tabulated with depth |

**Saturation laws** (`saturationLaw`, variably saturated only):

| Type | Coeffs sub-dict |
|---|---|
| `vanGenuchten` | `vanGenuchtenCoeffs` |
| `brooksCorey` | `brooksCoreyCoeffs` |
| `saturated` | — (fully saturated placeholder) |
| `exponential` | — |
| `montenegroMaier` | — |
| `pressureFunction` | — |

---

## `constant/poroCouplingProperties`

Required for coupled cases. Selects the coupling interface and shared-mesh behavior.

### Saturated coupled case — shared mesh

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      poroCouplingProperties;
}
// ----------------------------- //

// Use 'poroSolid' for saturated coupling,
// 'varSatPoroSolid' for variably saturated coupling
poroSolidInterface  poroSolid;

sharedMesh          true;
porosityConstant    true;   // false = update porosity from deformation

poroSolidCoeffs
{
    enabled         true;
    iterations      20;
    infoFrequency   1;

    convergence
    {
        p_rgh       max 1e-7;
        D           max 1e-7;
    }
}
```

### Non-shared mesh (separate fluid and solid regions)

When fluid and solid use separate meshes, add mapping configuration:

```foam
sharedMesh      no;
mapMethod       direct;     // 'direct' requires conforming meshes
consistent      yes;        // yes = meshes share the same boundary layout
```

For non-conforming meshes, use `mapMethod mapNearest`, set `consistent no`, and provide `patchMap` and `cuttingPatches` entries.

For variably saturated coupling (`varSatPoroSolid`), also set `porosityConstant false` and `porosityConstantExplicit true` (explicit porosity update from deformation).

---

## `constant/mechanicalProperties`

This is the standard solids4Foam format. For poromechanical coupling, use `poroMechanicalLaw2` as the outer law. It wraps any solids4Foam effective-stress law and applies the Biot coupling.

### Minimal coupled example

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      mechanicalProperties;
}
// ----------------------------- //

planeStress     no;

mechanical
(
    soil
    {
        type                    poroMechanicalLaw2;
        biotCoeff               biotCoeff [0 0 0 0 0 0 0] 1;
        pressureFieldName       p;      // 'p' (total) or 'p_rgh' (excess)

        // Any solids4Foam mechanical law goes here as the effective-stress law:
        effectiveStressMechanicalLaw
        {
            type    linearElastic;
            rho     rho [1 -3 0 0 0 0 0] 2000;
            E       E   [1 -1 -2 0 0 0 0] 1e7;
            nu      nu  [0 0 0 0 0 0 0]   0.3;
        }
    }
);
```

Without `poroMechanicalLaw2`, the Biot coefficient defaults to 1.0 everywhere and a warning is printed. Use any solids4Foam linear-geometry law as the `effectiveStressMechanicalLaw`.

### Available repository-owned mechanical laws

| Type | Notes |
|---|---|
| `poroMechanicalLaw2` | Biot coupling wrapper for any effective-stress law |
| `varSatPoroMechanicalLaw` | Variably saturated variant; adds effective-stress model |
| `linearElasticMohrCoulombPlasticDilationCutoff` | Elasto-plastic with dilation cutoff |
| `ohdeElastic` | Compression-index elastic law |
| `SANISAND` | Advanced sand model |

### Effective-stress models for `varSatPoroMechanicalLaw`

| Type | Notes |
|---|---|
| `terzaghi` | χ = S (simple Bishop) |
| `bishop` | χ = S |
| `niemunis` | Niemunis effective-stress formulation |
| `saturationFunction` | χ from tabulated S |
| `pressureFunction` | χ from tabulated p |
| `suctionCutOff` | suction limit |

---

## `constant/g`

All fluid regions require a gravity field:

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedVectorField;
    location    "constant";
    object      g;
}
// ----------------------------- //

dimensions  [0 1 -2 0 0 0 0];
value       (0 -9.81 0);
```

For shared-mesh coupled cases, gravity is read from the solid region automatically — a single `constant/g` file is sufficient.

---

## Typical Workflow

1. Build: `./Allwmake`
2. Prepare a case directory with the required dictionaries above.
3. Generate the mesh and initial fields.
4. Run `initPoroMechanicalFoam` if a pre-consolidated stress state is needed.
5. Run `poroMechanicalFoam`.
6. Inspect time directories and function-object output.

---

## `initPoroMechanicalFoam`

This utility initializes a geostatic or pre-consolidated state before the main solve. It:

- constructs a `poroFluidModel` and `solidModel`
- maps hydraulic fields to the solid mesh
- computes saturation for variably saturated cases
- ramps density over time to obtain a consolidated geostatic state

Use it when the solid stress initial condition depends on an initialized hydraulic state.

---

## Setup Utilities

The following utilities are in `solvers/utilities/` but **not built by the default `Allwmake`**. Build them manually in their own directories if needed.

| Utility | What it does |
|---|---|
| `setInitialGeostatics` | Writes initial `sigma` from unit weight `gamma` and `K0` |
| `setTotalHeadToPressureHead` | Converts total head to `pHead` or `p` |
| `setPressureHeadToPotential` | Converts `pHead` to `Potential` or `p_rgh` |
| `set_P_s_ToSigma` | Builds `sigma` from mean stress `P` and deviatoric stress `s` |
| `makeNewField` | Creates an empty volume field for case setup |

---

## Boundary Conditions

### Hydraulic boundary conditions

Custom types for hydraulic patches (apply to `p_rgh` or `pHead`):

| Type | Use case |
|---|---|
| `fixedPressureHead` | Fixed pressure head at a boundary |
| `fixedPoroPotential` | Fixed hydraulic potential |
| `fixedPoroFlux` | Specified flux (Neumann) |
| `seepageOutlet` | Seepage face: switches between no-flow and fixed head |
| `limitedHeadInfiltration` | Head-limited infiltration |
| `standingWaveTheory` | Oscillating wave loading |
| `emptyingTank` | Transient tank-draining head condition |

### Solid boundary conditions

Custom types for mechanical patches:

| Type | Use case |
|---|---|
| `geoStaticTraction` | Lithostatic/geostatic initial traction |
| `hydrostaticTraction` | Water-pressure traction |
| `poroTraction` | Pore-pressure traction (saturated) |
| `varSatPoroTraction` | Pore-pressure traction (variably saturated) |
| `fixedWallNoTension` | Contact/no-tension wall constraint |
| `zDependendDisplacementOrTraction` | Depth-dependent loading |
| `borePileFilling` | Bore-pile filling load |
| `cyclicPoroTotalTraction` | Cyclic total traction |

---

## Function Objects

Add these to `system/controlDict` under `functions {}` for runtime monitoring:

| Type | Purpose |
|---|---|
| `poroIncrements` | Increment-style poromechanics diagnostics |
| `poroPatchFlux` | Flow through selected patches |
| `secondOrderWork` | Internal work / instability indicator |
| `solidInvariants` | Stress or strain invariants |
| `stressChange` | Stress changes over time |
| `setTimeStep` | Runtime time-step control |
| `tabulatedForcingSource` | Tabulated forcing applied to the domain |

Dictionary syntax should be taken from the source files in `functionObjects/` — no standalone manuals exist yet.

---

## Common Setup Pitfalls

**Gravity is zero or missing**
The solver will abort with a "gravity is zero" error. Ensure `constant/g` exists.

**Cell zone names must match `poroHydraulicProperties`**
Each named cell zone in the mesh needs a matching sub-dict in `poroHydraulicProperties`. If no matching sub-dict is found, the solver falls back to the top-level entries of the file — so a flat file with no zone sub-dicts works as the default for all zones. If the mesh has no cell zones at all, the solver auto-creates `defaultZone` and the same fallback applies.

**`alpha` dimensions in van Genuchten**
The solver checks that `[1/alpha] == [p]` and aborts with a dimension error on mismatch. Use `[-1 1 2 0 0 0 0]` for pressure form (Pa), `[0 -1 0 0 0 0 0]` for head form (m). The same rule applies to `Ss`.

**No `poroMechanicalLaw2` found**
Without this law, there is no explicit hydraulic-to-solid coupling and the Biot coefficient defaults to 1.0. A warning is printed. If you need non-unity Biot coefficients, use `poroMechanicalLaw2`.

**Shared mesh vs separate mesh**
`sharedMesh yes` means both physics use the same OpenFOAM mesh region. `sharedMesh no` requires separate mesh directories under `constant/poroFluid/` and `constant/solid/`. Mismatches cause mesh-not-found errors at startup.

---

## Parallel Runs

```bash
decomposePar
mpirun -np <N> poroMechanicalFoam -parallel
reconstructPar
```

---

## Case Templates

The `caseTemplates/` directory contains starter cases for each supported workflow. Copy the closest one and adapt it:

| Template | Workflow |
|---|---|
| `poroFluid-saturated` | Saturated flow-only, `p_rgh` |
| `poroFluid-varSat-p_rgh` | Variably saturated flow-only, pressure form |
| `poroFluid-varSat-pHead` | Variably saturated flow-only, head form |
| `poroMechanical-saturated-sharedMesh` | Coupled saturated poromechanics, shared mesh |
| `poroMechanical-varSat-sharedMesh` | Coupled variably saturated poromechanics, shared mesh |

Each template includes a `README.md` listing what still needs to be adapted before running.

The dictionary examples in this guide use the same structure as these templates. If something looks different, the templates are the authoritative source.

## Starting Points

- Copy the closest template from `caseTemplates/` as your starting point.
- For new coupled cases, start with `poroMechanical-saturated-sharedMesh` — it has the fewest moving parts.
- For variably saturated work, the `poroFluid-varSat-p_rgh` or `poroMechanical-varSat-sharedMesh` template includes the SWCC setup.
- Use `initPoroMechanicalFoam` for stress-equilibrated initial conditions before the main solve.
- Background and worked examples are in the BAW report linked from `README.md`.

## Where To Look Next

- `caseTemplates/` — starter cases for each workflow
- `README.md` — install and build basics
- `docs/solids4FoamModelsMinimal.md` — minimal solids4Foam dependency details
- `docs/REPOSITORY_DOCUMENTATION.md` — maintainer reference: runtime architecture, extension points, full model hierarchy
