# `poroMechanicalFoam` Maintainer Reference

## Purpose

`poroMechanicalFoam` is an OpenFOAM extension for coupled porous-media flow and deformation. It is built around solids4Foam's `physicsModel` abstraction and adds:

- pore-fluid solvers for saturated and variably saturated flow
- poro-solid coupling interfaces that connect pore pressure, saturation, porosity, and deformation
- hydraulic constitutive models for storage, conductivity, and saturation laws
- porous mechanical laws that consume hydraulic state on the solid side
- boundary conditions, function objects, and setup utilities tailored to geotechnical workflows

The repository is structured as a hybrid of:

- repository-owned code in this tree
- a pinned minimal solids4Foam dependency in `external/solids4foam`
- OpenFOAM/solids4Foam runtime selection, so most model choices are made in case dictionaries rather than in solver source

## Supported Environment

The current build scripts explicitly support:

- OpenFOAM `v2412`
- OpenFOAM `v2512`

The top-level build expects a sourced OpenFOAM environment and a solids4Foam source tree available either through:

- the included submodule at `external/solids4foam`
- or an external tree pointed to by `S4F_ROOT`

The repository pins solids4Foam to tag `v2.3` and uses sparse checkout for only the required upstream sources.

## Build Overview

The top-level entry point is [`Allwmake`](../Allwmake). Its build order is:

1. optionally build `abaqusUMATs/` if `S4F_USE_GFORTRAN=1`
2. build `materialModels/poroHydraulicModel` as `libporoHydraulicModels`
3. initialize the `external/solids4foam` submodule if needed
4. build `solids4FoamModelsMinimal`
5. build the main repository library `libporoModels`
6. build the default solver applications under `solvers/`

The main library source list is defined in [`Make/files`](../Make/files). It compiles most of the repository-owned model code into:

- `$(FOAM_USER_LIBBIN)/libporoModels`

The main library links against:

- `libsolids4FoamModelsMinimal`
- `libporoHydraulicModels`
- OpenFOAM core libraries required by the fluid, solid, and coupling code

## Built Targets

### Default build path

The default `./Allwmake` path builds:

- `libporoHydraulicModels`
- `libsolids4FoamModelsMinimal`
- `libporoModels`
- `poroMechanicalFoam`
- `initPoroMechanicalFoam`

This follows [`solvers/Allwmake`](../solvers/Allwmake), which currently builds only:

- [`solvers/poroMechanicalFoam`](../solvers/poroMechanicalFoam/Allwmake)
- [`solvers/initPoroMechanicalFoam`](../solvers/initPoroMechanicalFoam/Allwmake)

### Present but not in the default solver build

The repository also contains:

- `poroMechanicalStrengthReduction`
- several small setup utilities under `solvers/utilities/`

These are not compiled by the default `solvers/Allwmake`. They must be built manually in their own directories if needed.

The repository already contains a bring-up note for the strength-reduction solver:

- [`docs/poroMechanicalStrengthReduction-bringup.md`](poroMechanicalStrengthReduction-bringup.md)

## High-Level Runtime Architecture

The main solver [`solvers/poroMechanicalFoam/poroMechanicalFoam.C`](../solvers/poroMechanicalFoam/poroMechanicalFoam.C) is deliberately thin. It:

1. creates `runTime`
2. constructs a runtime-selected `physicsModel`
3. advances the time loop
4. calls:
   - `physics().setDeltaT(runTime)` when adaptive time stepping is enabled
   - `physics().evolve()`
   - `physics().updateTotalFields()`
   - `physics().writeFields(runTime)` on output steps
5. calls `physics().end()` at shutdown

This means the executable itself does not hard-code whether the case is:

- dry or drained solid mechanics
- saturated groundwater flow
- coupled poromechanics
- variably saturated poromechanics

That decision is delegated to runtime-selected model families.

## Core Model Families

### 1. Top-level physics selection

The top-level abstraction is solids4Foam's `physicsModel`.

In this repository, the active `physicsModel` is usually one of:

- a solids4Foam `solidModel`
- a repository `poroFluidModel`
- a repository `poroSolidInterface`

The practical meaning is:

- `solidModel`: pure solid mechanics
- `poroFluidModel`: flow-only porous-media problem
- `poroSolidInterface`: coupled porous flow and deformation

### 2. Pore-fluid solver family

The pore-fluid hierarchy is documented in:

- [`poroFluidModels/CLASS_SYSTEM.md`](../poroFluidModels/CLASS_SYSTEM.md)

The current hierarchy is:

```text
poroFluidModel
├── poroFluid
└── variablySaturatedPoroFluid
    ├── varSatPoroFluid
    └── varSatPoroFluidHead
```

Key roles:

- `poroFluid`
  - saturated groundwater flow
  - solves with `p_rgh`
- `varSatPoroFluid`
  - Richards equation in pressure form
  - solves with `p_rgh`
- `varSatPoroFluidHead`
  - Richards equation in hydraulic-head form
  - solves with `pHead`

Shared responsibilities in `poroFluidModel`:

- mesh ownership for the fluid region
- field initialization (`p`, `p_rgh`, `phi`, hydraulic gradient, porosity)
- access to `poroHydraulicModel`
- fvOptions handling
- iteration control
- adaptive time-step support

### 3. Richards linearization family

Variably saturated fluid solvers delegate nonlinear treatment of Richards equation to a second runtime-selected family:

```text
richardsLinearization
├── Standard
├── Celia
├── LScheme
├── Casulli
└── Newton
```

`Newton` exists in the source tree, but the main library source list currently comments it out in [`Make/files`](../Make/files). That means it is present as code but not active in the default build.

### 4. Poro-solid coupling family

The coupling hierarchy is documented in:

- [`poroSolidInterfaces/classSystem.md`](../poroSolidInterfaces/classSystem.md)

The current hierarchy is:

```text
poroSolidInterface
├── poroSolid
└── varSatPoroSolid
```

`poroSolidInterface` owns the shared coupling infrastructure:

- construction of `solidModel`
- construction of `poroFluidModel`
- reading `constant/poroCouplingProperties`
- optional `meshToMesh` mapping
- Biot coefficient handling
- outer coupling iteration control
- installation of the `poroSolidToFluidCouplingSource` fvOption

Derived classes provide formulation-specific field exchange and source terms:

- `poroSolid`: saturated poromechanics
- `varSatPoroSolid`: variably saturated poromechanics

### 5. Hydraulic constitutive family

The hydraulic constitutive backend is split into:

- [`poroHydraulicModel`](../materialModels/poroHydraulicModel/poroHydraulicModel.H)
- [`varSatPoroHydraulicModel`](../materialModels/poroHydraulicModel/varSatPoroHydraulicModel.H)

Responsibilities:

- reading `constant/poroHydraulicProperties`
- providing porosity reference data
- assembling storage coefficients
- assembling intrinsic conductivity
- reconstructing hydrostatic reference pressure
- for variably saturated cases: saturation, relative permeability, effective conductivity, and moisture capacity

This separates governing-equation logic in the fluid solver from constitutive-law evaluation in the material model.

## Repository Contents by Directory

### `abaqusUMATs/`

Optional Abaqus UMAT interoperability build. The top-level build skips this unless `S4F_USE_GFORTRAN` is set.

### `docs/`

Repository documentation and notes:

- minimal solids4Foam dependency note
- strength-reduction bring-up note
- Doxygen configuration
- this consolidated documentation

### `functionObjects/`

Custom OpenFOAM function objects for poromechanics diagnostics and control:

- `poroIncrements`
- `poroPatchFlux`
- `secondOrderWork`
- `setTimeStep`
- `solidInvariants`
- `stressChange`
- `tabulatedForcingSource`

These are intended for runtime monitoring, post-processing, or forcing/control workflows.

### `iterationControls/`

Outer and inner convergence control components used by the fluid and coupled solvers.

It contains:

- `iterationControl`
- `iterationResidual` types
- residual reduction operations

This is the convergence-control layer used by coupled fixed-point style solves and nonlinear Richards loops.

### `materialModels/`

Repository-owned constitutive code, split into:

- porous hydraulic laws
- additional mechanical laws that extend solids4Foam

#### `materialModels/poroHydraulicModel/`

Hydraulic constitutive backend. It includes runtime-selectable families for:

- saturation laws
- storage laws
- conductivity laws
- porosity submesh helpers

#### `materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/`

Repository mechanical laws and poromechanical wrappers:

- `linearElasticMohrCoulombPlasticDilationCutoff`
- `ohdeElastic`
- `SANISAND`
- `poroMechanicalLaw2`
- `varSatPoroMechanicalLaw`
- `MFront` prototype sources

`varSatPoroMechanicalLaw` also introduces a nested effective-stress model family.

### `poroFluidModels/`

Pore-fluid solver family, Richards linearization family, and custom hydraulic boundary conditions.

### `poroSolidInterfaces/`

Poro-solid coupling family, coupling fvOption, and custom solid-side boundary conditions for traction/displacement and geotechnical loading.

### `solidModels/`

Repository-owned solids4Foam-compatible solid models. At present this contains:

- `seismicLinGeomSolid`

### `solids4FoamModelsMinimal/`

Minimal solids4Foam build metadata that compiles a reduced upstream subset into `libsolids4FoamModelsMinimal`.

### `solvers/`

Top-level applications:

- `poroMechanicalFoam`
- `initPoroMechanicalFoam`
- `poroMechanicalStrengthReduction`
- small utilities under `solvers/utilities/`

## Runtime-Selectable Hydraulic Laws

The hydraulic constitutive family is assembled from zone-wise runtime-selected laws.

### Saturation laws

Available saturation law types in `materialModels/poroHydraulicModel/saturationLaws/`:

- `brooksCorey`
- `exponential`
- `montenegroMaier`
- `pressureFunction`
- `saturated`
- `vanGenuchten`

### Storage laws

Available storage law types in `materialModels/poroHydraulicModel/storageLaws/`:

- `KPrime`
- `montenegro`
- `pressureFunction`
- `skemptonB`
- `storageCoeff`

### Conductivity laws

Available conductivity model types in `materialModels/poroHydraulicModel/conductivityModels/`:

- `anisotropicConstant`
- `constant`
- `gradientFunction`
- `kozenyCarman`
- `limitedGradient`
- `porosityFunction`

## Runtime-Selectable Mechanical and Effective-Stress Laws

### Additional mechanical laws in this repository

Available repository-owned mechanical law types:

- `linearElasticMohrCoulombPlasticDilationCutoff`
- `ohdeElastic`
- `SANISAND`
- `poroMechanicalLaw2`
- `varSatPoroMechanicalLaw`

These extend the solids4Foam mechanical-law family and are used from `constant/mechanicalProperties`.

`MFront` code exists in the source tree as an initial prototype, but it is not wired into a supported or documented workflow and should not be treated as an available user-facing mechanical law.

### Effective-stress models for `varSatPoroMechanicalLaw`

Available effective-stress model types:

- `bishop`
- `niemunis`
- `pressureFunction`
- `saturationFunction`
- `suctionCutOff`
- `terzaghi`

## Boundary Conditions and Patch Fields

### Pore-fluid patch fields

Custom hydraulic patch field types under `poroFluidModels/poroHydraulicFvPatchFields/`:

- `codedFixedValue`
- `emptyingTank`
- `fixedPoroFlux`
- `fixedPoroPotential`
- `fixedPressureHead`
- `limitedHeadInfiltration`
- `seepageOutlet`
- `standingWaveTheory`
- `timeDependentOutletWall`

There is also an older `limitedHeadInfiltration` implementation in `oldLimitedHead/`.

### Poro-solid / solid patch fields

Custom solid-side patch field types under `poroSolidInterfaces/poroSolidFvPatchFields/`:

- `borePileFilling`
- `cyclicPoroTotalTraction`
- `fixedWallNoTension`
- `geoStaticTraction`
- `hydrostaticTraction`
- `poroTraction`
- `varSatPoroTraction`
- `zDependendDisplacementOrTraction`

These encode geotechnical loading and support conditions beyond stock solids4Foam behavior.

## Coupling Source Term

The deformation-to-flow coupling enters the fluid equations through:

- [`poroSolidInterfaces/poroSolidToFluidCouplingSource`](../poroSolidInterfaces/poroSolidToFluidCouplingSource/poroSolidToFluidCouplingSource.H)

This is implemented as an `fvOption`. It looks up the active `poroSolidInterface` and asks it for explicit and implicit deformation-to-pressure source contributions. This keeps the fluid equation assembly generic while allowing each poro-solid formulation to define its own coupling term.

## Function Objects

The repository adds the following function objects:

- `poroIncrements`
  - writes increment-style poromechanics diagnostics
- `poroPatchFlux`
  - reports fluxes on selected patches
- `secondOrderWork`
  - evaluates internal work / instability-related quantity
- `setTimeStep`
  - runtime time-step control helper
- `solidInvariants`
  - post-processes stress or strain invariants
- `stressChange`
  - reports stress changes
- `tabulatedForcingSource`
  - applies tabulated forcing to a case

Exact dictionary syntax should be taken from the source files in `functionObjects/` because no standalone usage manual for them is yet included in this repository.

## Iteration and Residual Infrastructure

The repository includes a reusable convergence-control layer:

- `iterationControl`
  - manages iteration stopping logic
- `iterationResidual`
  - abstract residual source
- residual operations:
  - `L2`
  - `RMS`
  - `max`
  - `sum`
- residual measures:
  - `delta`
  - `linearSolver`
  - `MassBalance`

This infrastructure is used by:

- the coupled poro-solid outer loop
- nonlinear variably saturated fluid solves
- adaptive stopping based on residuals configured in dictionaries

## Applications and Utilities

### `poroMechanicalFoam`

Main production solver. It creates a runtime-selected `physicsModel` and advances the time loop. It is the general-purpose entry point for:

- flow-only cases
- solid-only cases using the same framework
- coupled poromechanical cases

### `initPoroMechanicalFoam`

Initialization workflow for stress and hydraulic state. The source in [`solvers/initPoroMechanicalFoam/initPoroMechanicalFoam.C`](../solvers/initPoroMechanicalFoam/initPoroMechanicalFoam.C) shows that it:

- constructs a `poroFluidModel`
- constructs a `solidModel`
- maps hydraulic fields from fluid mesh to solid mesh
- computes saturation for variably saturated cases
- recalculates density on the solid side
- ramps density over time to obtain an initialized geostatic/consolidated state

### `poroMechanicalStrengthReduction`

Strength-reduction solver prototype in [`solvers/poroMechanicalStrengthReduction/poroMechanicalStrengthReduction.C`](../solvers/poroMechanicalStrengthReduction/poroMechanicalStrengthReduction.C). Current characteristics:

- marked in source as unfinished work
- updates `frictionAngle`, `cohesion`, and `dilationAngle` during progression
- assumes certain fields such as `TotalHead` and `Sw` exist
- not wired into the default solver build

The current repository note in `docs/poroMechanicalStrengthReduction-bringup.md` should be treated as the authoritative status summary.

### Setup utilities

Utilities present under `solvers/utilities/`:

- `makeNewField`
  - creates an empty volume field of a chosen type
- `setInitialGeostatics`
  - writes an initial geostatic `sigma` field from `gamma` and `K0`
- `setPressureHeadToPotential`
  - converts `pHead` into `Potential` or `p_rgh`
- `setTotalHeadToPressureHead`
  - converts total head into `pHead` or `p`
- `set_P_s_ToSigma`
  - builds `sigma` from mean stress `P` and deviatoric stress `s`

These utilities are useful for case preparation, but they are not part of the default top-level build.

## Build and Dependency Notes

### solids4Foam dependency strategy

The repository does not rely on a full solids4Foam build. Instead it builds a reduced library from a sparse-checked-out submodule:

- submodule path: `external/solids4foam`
- tag: `v2.3`
- sparse scope:
  - `src/solids4FoamModels`
  - `src/blockCoupledSolids4FoamTools`

See:

- [`docs/solids4FoamModelsMinimal.md`](solids4FoamModelsMinimal.md)

### Submodule behavior

The submodule is declared in [`.gitmodules`](../.gitmodules) with:

- `path = external/solids4foam`
- `url = https://github.com/solids4foam/solids4foam.git`
- `sparseCheckout = true`
- `ignore = dirty`

## Case-Level Expectations

Although this repository does not ship a full tutorial set, the source tree makes the main case-level expectations clear:

- `constant/poroHydraulicProperties`
  - selects hydraulic constitutive laws
- `constant/poroCouplingProperties`
  - selects the poro-solid interface and coupling options
- `constant/mechanicalProperties`
  - selects the solid mechanical law
- `system/controlDict`
  - chooses the application and function objects
- `system/fvSolution`
  - may contain iteration-control configuration

For coupled cases, the solid-side constitutive laws expect hydraulic fields with standard names such as:

- `p`
- `p_rgh`
- `S`
- `n`

The coupling layer either registers these directly on the solid mesh for shared-mesh cases or maps them there for non-shared meshes.

## Extension Points

### Add a new pore-fluid solver

1. Derive from `poroFluidModel` or `variablySaturatedPoroFluid`.
2. Register it with `addToRunTimeSelectionTable(poroFluidModel, ..., dictionary)`.
3. Implement `pField()`, `evolve()`, and time-step logic.
4. Add the source to [`Make/files`](../Make/files).

### Add a new Richards nonlinear strategy

1. Derive from `richardsLinearization`.
2. Register it in that runtime-selection table.
3. Implement the linearization/update interface.
4. Add it to the library source list.

### Add a new poro-solid coupling formulation

1. Derive from `poroSolidInterface`.
2. Implement the coupling hooks used by `evolveCouplingLoop()`.
3. Register it with `addToRunTimeSelectionTable(poroSolidInterface, ..., dictionary)`.
4. Add the new source to the build.

### Add a new hydraulic law

1. Derive from `saturationLaw`, `storageLaw`, or `conductivityModel`.
2. Register the law.
3. Ensure the relevant dictionary parsing is supported.
4. Add the source to `materialModels/poroHydraulicModel/Make/files`.

### Add a new mechanical law or effective-stress model

1. Derive from the relevant solids4Foam law base class.
2. Register it with the correct linear-geometry or effective-stress selection table.
3. Add the source to the root [`Make/files`](../Make/files).

## Current Documentation Gaps and Known Caveats

These points are clear from the repository state and should be treated as current limitations:

- no bundled tutorial cases are present in this repository
- function objects and patch fields do not yet have dedicated user manuals
- `poroMechanicalStrengthReduction` is present but not part of the default build and is explicitly marked unfinished
- `MFront` is present only as an initial prototype and is not wired into a supported workflow
- some source comments still reflect older foam-extend wording even though the active build target is OpenFOAM `v2412`/`v2512`
- some utilities and auxiliary targets must be built manually

## Recommended Reading Order for New Maintainers

1. [`README.md`](../README.md)
2. [`Allwmake`](../Allwmake)
3. [`poroFluidModels/CLASS_SYSTEM.md`](../poroFluidModels/CLASS_SYSTEM.md)
4. [`poroSolidInterfaces/classSystem.md`](../poroSolidInterfaces/classSystem.md)
5. [`docs/solids4FoamModelsMinimal.md`](solids4FoamModelsMinimal.md)
6. the solver entry points in `solvers/`

That path gives the fastest route from build setup to runtime architecture.
