# poroFluidModels class system

This folder implements the pore-fluid side of `poroMechanicalFoam`. The code is organized as a small runtime-selectable class hierarchy plus a second runtime-selectable hierarchy for Richards-equation linearization.

## 1. Top-level runtime selection

`poroFluidModel` is the abstract root class. It provides:

- common mesh and registry handling
- access to `poroFluidProperties`
- shared fields such as `p`, `p_rgh`, `phi`, `hydraulicGradient`, `n`
- time-step control and outer-iteration control
- the runtime-selection entry point `poroFluidModel::New(...)`

Concrete implementations are registered with:

```cpp
addToRunTimeSelectionTable(poroFluidModel, SomeClass, dictionary);
```

That means the active pore-fluid solver is chosen from dictionaries at run time rather than hard-coded in the solver executable.

## 2. Solver hierarchy

The currently relevant inheritance tree is:

```text
poroFluidModel
├── poroFluid
└── variablySaturatedPoroFluid
    ├── varSatPoroFluid
    └── varSatPoroFluidHead
```

### `poroFluid`

Use this class for saturated groundwater flow with `p_rgh` as the solved variable.

- primary unknown: `p_rgh`
- constitutive backend: `poroHydraulicModel`
- governing character: linear storage + Darcy diffusion
- typical role: simplest groundwater solver and coupling target when saturation changes are not needed

### `variablySaturatedPoroFluid`

This is an abstract intermediate base for Richards-equation solvers.

It adds the machinery that both variably saturated solvers need:

- `varSatPoroHydraulicModel` for saturation-dependent constitutive laws
- pore-fluid saturation field `S`
- an optional mass-balance residual field
- a runtime-selected `richardsLinearization` object

It does not choose the primary pressure variable itself. That decision is pushed down to its concrete children.

### `varSatPoroFluid`

This is the Richards solver written in excess-pressure form.

- primary unknown: `p_rgh`
- reconstructed total pressure: `p = p_rgh + p_Hyd`
- gravity handling: implicit through the hydrostatic decomposition
- face transport field: `kEffbyGammaf`

This version is convenient when the rest of the code base already works naturally with `p_rgh`.

### `varSatPoroFluidHead`

This is the Richards solver written in hydraulic-head form.

- primary unknown: `pHead`
- reconstructed total pressure: `p = gamma * pHead`
- gravity handling: explicit via `z` terms in the equation
- face transport field: `kEfff`

This version is convenient when hydraulic head is the more natural user-facing quantity or when boundary conditions are specified directly in head units.

## 3. Secondary runtime selection: Richards linearization

The variably saturated branch delegates the nonlinear treatment of Richards equation to a second class family:

```text
richardsLinearization
├── Standard
├── Celia
├── LScheme
├── Casulli
└── Newton   (present in tree, currently not enabled in Make/files)
```

`variablySaturatedPoroFluid` constructs this object from:

```cpp
poroFluidDict().get<word>("solutionAlgorithm")
```

The linearization object is responsible for:

- initializing iteration-specific coefficient fields
- contributing the linearized saturation term through `ddtS(...)`
- checking convergence of the inner nonlinear loop
- updating auxiliary coefficient `C`

So the variably saturated solver classes own the equation structure and field bookkeeping, while `richardsLinearization` owns the specific nonlinear update strategy.

## 4. Constitutive-model split

There is also a constitutive split between saturated and variably saturated solvers:

- `poroFluid` uses `poroHydraulicModel`
- `variablySaturatedPoroFluid` children use `varSatPoroHydraulicModel`

This keeps the governing-equation code separate from the material-law details:

- storage law
- conductivity law
- relative permeability / effective conductivity
- hydrostatic pressure reconstruction
- saturation law

## 5. Primary-field convention

`poroFluidModel` exposes a generic `pField()` interface so coupling code does not need to know whether the active solver is head-based or pressure-based.

- `poroFluid::pField()` -> `p_rgh`
- `varSatPoroFluid::pField()` -> `p_rgh`
- `varSatPoroFluidHead::pField()` -> `pHead`

Use `pField()` when writing generic code against the class hierarchy.
Use `p()`, `p_rgh()`, or solver-specific fields only when the formulation matters.

## 6. Practical extension points

To add a new pore-fluid solver:

1. Derive from `poroFluidModel` or `variablySaturatedPoroFluid`.
2. Register it with `addToRunTimeSelectionTable(...)`.
3. Decide which field is the primary unknown by implementing `pField()`.
4. Implement `evolve()`, `newDeltaT()`, and any additional write-out logic.

To add a new nonlinear strategy for Richards equation:

1. Derive from `richardsLinearization`.
2. Register it in the linearization runtime-selection table.
3. Implement `initalize(...)`, `checkConvergedAndUpdate(...)`, and `ddtS(...)`.

That separation is the key design choice in this folder: solver formulation and nonlinear strategy are selected independently.
