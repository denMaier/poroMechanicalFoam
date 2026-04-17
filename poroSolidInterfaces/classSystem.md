# poroSolidInterfaces class system

This folder contains the coupling layer between the poro-fluid solver family and the solids4Foam solid solver family. The classes here do not implement the hydraulic constitutive laws or the solid constitutive laws themselves. Instead, they orchestrate how both sides exchange fields and in which order they are solved.

## Runtime selection structure

`poroSolidInterface` is the abstract base class.

- It derives from `physicsModel` so it can be selected and evolved like the other top-level models in the solver stack.
- It also derives from `IOdictionary` and reads `constant/poroCouplingProperties`.
- The concrete implementation is selected at runtime through the `poroSolidInterface` entry in that dictionary.

Current derived classes in this folder:

- `poroSolid`
  - fully saturated formulation
- `varSatPoroSolid`
  - variably saturated formulation with additional saturation and porosity state mirrored onto the solid side

The base class also registers itself as a `physicsModel` named `poroSolid`, which is what the solver layer instantiates.

## Ownership and responsibilities

`poroSolidInterface` owns the cross-model infrastructure:

- construction of the `solidModel`
- construction of the `poroFluidModel`
- reading coupling options from `poroCouplingProperties`
- optional `meshToMesh` mapper for non-shared meshes
- extraction and caching of the Biot coefficient field `b`
- creation of the `iterationControl` object for the outer coupling loop
- installation of the `fv::poroSolidToFluidCouplingSource` into the fluid model

The derived classes own formulation-specific transient fields:

- explicit deformation-to-flow source `nDot`
- implicit fixed-stress stabilization coefficient
- mapped hydraulic fields required by the solid-side constitutive law
- any extra state needed by the chosen formulation

## Coupling path

The coupling is split into two directions.

### Solid to fluid

The deformation contribution enters the flow equation through the `fv::poroSolidToFluidCouplingSource` fvOption.

The source object does not calculate the terms itself. Instead, it looks up the active `poroSolidInterface` instance and asks it for:

- `explicitCouplingDtoP()`
- `implicitCouplingDtoP()`

This keeps the fvOption generic while allowing each derived class to assemble its own formulation-specific source terms.

### Fluid to solid

The fluid-to-solid coupling is consumed indirectly by the porous mechanical laws on the solid side.

Those material laws expect to find hydraulic fields with fixed names in the solid mesh object registry, primarily:

- `p`
- `p_rgh`
- `S` for variably saturated mechanics
- `n` where total density updates require porosity on the solid side

For shared meshes, the interface registers the fluid fields directly on the solid mesh registry.

For separate meshes, the interface creates and maintains mapped copies on the solid mesh, so the constitutive laws can remain agnostic to the mesh layout.

## Shared mesh vs mapped mesh

The `sharedMesh` switch in `poroCouplingProperties` controls the data path.

When `sharedMesh true`:

- no `meshToMesh` object is needed
- hydraulic fields can be checked directly into the solid mesh registry
- deformation fields can be used directly when assembling fluid coupling terms

When `sharedMesh false`:

- the base class constructs a `meshToMesh` interpolator
- solid kinematic fields are mapped to the fluid mesh when assembling flow sources
- fluid pressures, and in the variably saturated case also saturation and porosity, are mapped back to the solid mesh

This is the main reason the derived classes hold mapped field pointers such as `pSolidMesh_`, `pRghSolidMesh_`, `SSolidMesh_`, and `nSolidMesh_`.

## Outer iteration lifecycle

The actual time-step work is implemented in the derived `evolve()` methods, but both formulations follow the same high-level sequence:

1. reset the outer coupling residual controller
2. assemble deformation-to-flow coupling fields from the current solid state
3. optionally update fluid porosity before the fluid solve
4. evolve the poro-fluid model
5. map hydraulic state back to the solid side when the meshes differ
6. evolve the solid model
7. repeat until `iterationControl` declares convergence
8. update porosity for the finished time step
9. write coupling diagnostics if enabled

The main difference between the two derived classes is which fields are needed for steps 2 and 5.

## Formulation-specific differences

### `poroSolid`

The saturated implementation uses:

- Biot coefficient `b`
- solid velocity or mapped solid velocity
- solid tangent bulk stiffness `impK`
- solid displacement for porosity updates

It mirrors only pressure-related fields onto the solid side.

### `varSatPoroSolid`

The variably saturated implementation adds:

- saturation `S`, which weights the deformation-to-flow coupling
- porosity `n` mirrored onto the solid side for density updates
- a pressure conversion factor `magGammaW_` so the same class can work with pressure-based and head-based poro-fluid solvers

It also triggers `solidRef().recalculateRho()` during the outer loop so total density can react to updated hydraulic state.

## Extension points

To add another poro-solid formulation in this folder:

1. derive a new class from `poroSolidInterface`
2. register it with `addToRunTimeSelectionTable(poroSolidInterface, ..., dictionary)`
3. implement `evolve()`
4. implement `explicitCouplingDtoP()`
5. implement `implicitCouplingDtoP()`
6. decide which hydraulic fields must be exposed on the solid registry for the mechanical law
7. update `Make/files`

If the new formulation needs additional field exchange, the preferred pattern in the current design is:

- keep the base class responsible for shared infrastructure
- keep the derived class responsible for formulation-specific transient fields and mapped field ownership
