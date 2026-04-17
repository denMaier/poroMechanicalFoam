# thoughts

Potential robustness improvements noticed while documenting the `poroSolidInterfaces` classes.

## 1. Validate field availability earlier and more uniformly [done]

The code already contains useful runtime checks, but they are spread across constructors and solver methods.

It would be more robust to validate all required fields up front during construction, especially for:

- `S` in the variably saturated case
- `gamma_water` for head-based formulations
- mapped solid-side field allocation when `sharedMesh` is false

Implemented:

- lazy first-use validation for `S` in `varSatPoroSolid`
- lazy first-use validation for `gamma_water` when head units are actually used
- lazy initialization of mapped solid-side hydraulic fields before the first coupled solve

That now fails on first actual use instead of during construction, which is safer for lazy initialization while still improving diagnostics.

## 2. Consolidate mapped-field synchronization [done]

In the non-shared-mesh path, mapped copies of `p`, `p_rgh`, `S`, and `n` are refreshed in several places. The current approach works, but it spreads synchronization rules across constructors and `evolve()`.

Implemented:

- `poroSolid::syncSolidHydraulicFields()` for non-shared pressure refresh
- `varSatPoroSolid::syncSolidHydraulicFields()` for non-shared pressure and saturation refresh
- `varSatPoroSolid::syncSolidPorosityField()` for non-shared porosity refresh

That now keeps the non-shared mapped-field update paths localized instead of spreading them across the coupling loop.

## 3. Reduce duplication between `poroSolid` and `varSatPoroSolid`

Both derived classes duplicate most of the outer coupling loop:

- shared vs non-shared mesh branches
- `nDot_` assembly
- fixed-stress stabilization assembly
- porosity updates
- pressure mapping back to the solid mesh

This makes small behavior changes easy to apply in one class but miss in the other. A more robust design would move the common loop skeleton into `poroSolidInterface` and expose only formulation-specific hooks such as:

- assembleExplicitCoupling()
- assembleImplicitCoupling()
- syncHydraulicStateToSolid()
- postFluidSolve()

## 4. Strengthen type-based Biot-coefficient discovery

`poroSolidInterface::checkMechanicalLawUpdateBiotCoeff()` currently mixes:

- `dynamic_cast` for `poroMechanicalLaw2`
- string-based `law.type()` checks for `poroMechanicalLaw`

That is brittle if more porous mechanical laws are added.

A more robust design would expose a small virtual interface for mechanical laws that can provide:

- Biot coefficient
- whether hydraulic coupling is active

Then the interface layer would not need to know specific implementation class names.

## 5. Make object-registry ownership more explicit

The shared-mesh path relies on `objectRegistry::checkIn(...)` for fluid-owned fields so the solid-side constitutive laws can find them. That is convenient, but it is also fragile because ownership and lifetime are implicit.

It would be more robust to centralize that registration logic in helper functions and document:

- which registry owns the field
- whether the field may already be registered
- whether repeated construction/destruction in restart scenarios is safe

## 6. Clarify the saturation weighting used in stabilization

In `varSatPoroSolid`, the explicit coupling term uses a saturation-weighted Biot coefficient, while the fixed-stress stabilization term is still built from `b()` and `impK`.

That may be intentional, but the rationale is not obvious from the code. A comment or a small formulation note would make this safer to maintain. If the intent is actually to weight stabilization with saturation as well, this deserves a closer numerical review before any change.

## 7. Revisit the currently disabled `q_relAcc_` contribution

Both derived classes keep a `q_relAcc_` field and a helper that computes it, but the term is commented out in the active coupling assembly.

That leaves the code in an ambiguous state:

- it looks supported
- it is allocated as part of the design
- it is not actually used

More robust options would be:

- remove the dead path until it is validated
- or guard it with an explicit runtime switch and a documented experimental status
