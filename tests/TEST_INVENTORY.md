# Test Inventory

This inventory records the current test surface and flags checks that are
mostly structural or tautological. A tautological check here means an assertion
that mainly re-expresses the implementation body, constructor wiring, or a
literal repository policy, rather than independently validating externally
observable behavior.

## Python Repository Checks

### `tests/test_build_manifests.py`

- `test_active_sources_in_makefiles_exist`: verifies active source paths listed
  in OpenFOAM `Make/files` manifests exist.
- `test_makefiles_do_not_build_template_placeholders`: verifies active
  makefiles do not include template placeholder sources.
- `test_default_solver_build_matches_documented_targets`: checks current
  `solvers/Allwmake` policy: build `poroMechanicalFoam` and
  `initPoroMechanicalFoam`, exclude strength reduction from the default build.
- `test_top_level_build_rejects_unsupported_openfoam_versions`: checks
  `Allwmake` rejects unsupported OpenFOAM versions.

Low-value or tautological candidates:

- `test_default_solver_build_matches_documented_targets` and
  `test_top_level_build_rejects_unsupported_openfoam_versions` are policy
  string checks. They are useful as guardrails, but they do not exercise build
  behavior.

### `tests/test_case_templates.py`

- `test_every_case_template_has_minimum_openfoam_layout`: verifies every case
  template has the expected OpenFOAM directory and dictionary skeleton.
- `test_case_templates_use_the_main_solver`: verifies template `controlDict`
  files select `poroMechanicalFoam`.
- `test_template_runtime_model_names_are_registered_in_repo`: verifies model
  names referenced in template `*Properties` files have matching `TypeName`
  declarations somewhere in the repo.
- `test_flow_templates_have_the_field_selected_by_their_model`: verifies known
  flow templates include the selected primary pressure/head field.

Low-value or tautological candidates:

- `test_case_templates_use_the_main_solver` is a direct policy check and may
  need updating if solver-specific templates are intentionally added.
- `test_flow_templates_have_the_field_selected_by_their_model` is useful but
  hard-coded to three known templates, so it can miss new templates unless the
  expectation table is maintained.

### `tests/test_runtime_selection.py`

- `test_compiled_runtime_selected_classes_define_type_names`: verifies compiled
  sources that register with OpenFOAM runtime selection also call
  `defineTypeNameAndDebug`.
- `test_runtime_selected_classes_have_matching_header_typename`: verifies
  runtime-selected sources have a neighboring header declaring the registered
  class and containing `TypeName`.
- `test_richards_linearization_default_build_excludes_newton_consistently`:
  verifies the default library build includes `Standard` and excludes `Newton`.
- `test_registered_type_names_are_unique_within_model_family_directories`:
  verifies runtime `TypeName` values are unique inside major model families.

Low-value or tautological candidates:

- `test_richards_linearization_default_build_excludes_newton_consistently` is a
  literal build-manifest policy check, not behavioral coverage.

## OpenFOAM Unit Drivers

### `tests/openfoam/hydraulicModelUnit/hydraulicModelUnit.C`

Tested code:

- Saturation laws:
  - `vanGenuchten`: `pStar`, `S`, `kr`, `C`, saturated-pressure behavior,
    physical bounds, monotonicity, non-negative capacity.
  - `brooksCorey`: `pStar`, `S`, `kr`, `C`, air-entry behavior, above-entry
    behavior, physical bounds, monotonicity, non-negative capacity.
- Conductivity models:
  - `kozenyCarman`: permeability value and positivity.
- Storage laws:
  - `storageCoeff`: static specific storage and `updatesSs`.
  - `KPrime`: pressure-independent and pressure-dependent storage, pressure
    trend, `updatesSs`.
  - `skemptonB`: specific storage, static behavior, positivity.
  - `montenegro`: above-entry storage and cutoff behavior, `updatesSs`.
- `poroHydraulicModel` material-zone assembly:
  - real two-cell, two-zone OpenFOAM mesh created through `topoSet`.
  - zone-specific `n`, `Ss`, and `k` are assembled onto the expected cells.
  - constant zone laws report no state-dependent storage or conductivity update.

Low-value or tautological candidates:

- `storageCoeff preserves runtime name` only confirms the selected object name
  is retained after constructing the object with that same runtime name.
- Simple `updatesSs()` assertions are useful API contract checks, but they are
  nearly tautological for implementations that return a hard-coded boolean.
- Positivity/bounds assertions that follow immediately from stronger golden
  values are redundant but harmless. Keep them only if they improve diagnostic
  clarity.

### `tests/openfoam/poroSolidInterfaceUnit/poroSolidInterfaceUnit.C`

Tested code:

- `poroSolid` shared-mesh fixture construction.
- Default shared-mesh cell-zone creation.
- Biot coefficient loading through the solid mechanical law.
- Solid-side registration of hydraulic fields `p` and `p_rgh`.
- Object identity for shared `p` and `p_rgh` fields.
- Saturated coupling-term assembly:
  - explicit DTO dimensions and zero value for uniform displacement.
  - implicit DTO dimensions and `b^2 / impK` value.

Low-value or tautological candidates:

- `fixture uses shared mesh`, default cell-zone existence, and default
  cell-zone ownership are mostly fixture smoke tests. They are still useful if
  shared-mesh construction has historically been fragile.
- `Biot coefficient from fixture` mostly confirms the test fixture dictionary
  value is read. It becomes more useful when paired with coupling-term
  assertions that depend on `b`.

### `tests/openfoam/varSatPoroSolidInterfaceUnit/varSatPoroSolidInterfaceUnit.C`

Tested code:

- `varSatPoroSolid` shared-mesh fixture construction.
- Runtime selection of the variably saturated fluid model.
- Default shared-mesh cell-zone creation.
- Biot coefficient loading through `poroMechanicalLaw2`.
- `poroFluid::S()` delegation to `variablySaturatedPoroFluid::S()`.
- Solid-side registration and object identity for shared `p`, `p_rgh`, `S`,
  and `n` fields.
- Idempotence of hydraulic field re-registration.
- Variably saturated coupling assembly:
  - explicit DTO uses `S*b*div(U)`.
  - explicit DTO is not pure `b*div(U)`.
  - implicit DTO uses `b^2 / impK`.
  - implicit DTO is not saturation-weighted.
  - re-assembly updates explicit saturation-dependent terms in place.
  - `clearTerms()` permits re-assembly.
- Transfer of interface terms through the fluid `fvOptions` source:
  - implicit term goes to matrix diagonal over `deltaT`.
  - explicit term goes to matrix source.
  - source updates when saturation changes.
- Relative acceleration path:
  - accelerating displacement adds a transient source.
  - `afterFluidSolve()` clears the transient acceleration contribution while
    keeping the `nDot` and implicit DTO terms available.
- `Allrun` fatal-path coverage:
  - explicit and implicit DTO access before coupling-term assembly fails with
    the intended diagnostics.
  - `varSatPoroSolid` rejects a saturated `poroFluid` configuration before
    running with an invalid variably saturated coupling setup.

Low-value or tautological candidates:

- Shared-mesh, default cell-zone, and Biot coefficient checks are fixture smoke
  tests. They overlap with the saturated interface unit.
- The `fvOption` transfer assertions include independent fixture values for
  the saturation-dependent explicit source and saturation-independent implicit
  diagonal. They still also check DTO handoff wiring.
- `seeded velocity is returned through solid().U()` validates the test helper
  more than production behavior.

### `tests/openfoam/couplingTermsUnit/couplingTermsUnit.C`

Tested code:

- `poroCouplingTerms`:
  - `nDot` is named correctly and is zero for rigid translation.
  - `fixedStressStabil` computes `b^2/K`, stays positive, and has inverse
    pressure dimensions.
  - variably saturated weighting decision: `nDot` uses `S*b`, fixed-stress
    stabilization uses `b`.
  - `updateCouplingFields` creates, registers, reuses, and updates coupling
    fields.
  - `relativeAccelerationFlux` projects acceleration to flux dimensions.
  - `explicitCouplingSource` subtracts relative acceleration divergence.
  - `implicitCouplingRate`, `explicitCouplingRate`, and `addCouplingSource`
    apply the expected scaling/signs to matrix terms.
- `residualOperation`: `max`, `sum`, `L2`, `RMS`, including empty RMS.
- `deltaVf`: field discovery, missing-field rejection, absolute and relative
  residuals, vector and surface-field residuals.
- `iterationControl`: future-window disablement, active window parsing,
  interval parsing, residual creation, residual calculation.
- `LinearSolverRes`: scalar residual extraction, tolerance parsing, vector
  max-component residual, `show` mode.
- `MassBalanceTerms`: absolute and relative residual values, zero relative
  denominator behavior.
- `scalarDiskReader`: initial read, dimensions, write option, reread on time
  index change, missing-start-field fatal diagnostic.
- `poroCouplingRegistry`: borrowed-field registration, identity preservation,
  idempotence, duplicate-field fatal diagnostic.
- `poroPressureUnits`: pressure-dimensional scale, head-dimensional
  `gamma_water` scale, missing-`gamma_water` fatal diagnostic.
- Effective-stress models: `terzaghi`, `bishop`, `niemunis`, `suctionCutOff`,
  dimension and boundary behavior.
- `varSatPoroMechanicalLawTerms`: mixture density, initial effective stress,
  total stress reconstruction, dimensions and boundary behavior.

Low-value or tautological candidates:

- `nDot field name` only verifies the literal `IOobject` name set in
  `poroCouplingTerms::nDot`.
- `fixedStressStabil b^2/K golden`, `implicitCouplingRate`,
  `explicitCouplingRate`, and `addCouplingSource` are close to the production
  one-line implementations. They are still useful sign/dimension regression
  checks, but each independently has limited fault-finding power.
- `relativeAccelerationFlux boundary z` is largely a direct restatement of
  `kf * interpolate(a / magG)` for a uniform field.
- `pressure unit scale is one for pressure fields` and the head-scale check are
  very close to the implementation branches; the fatal diagnostic checks in
  `Allrun` add more meaningful guard coverage.
- `shared registry inserts borrowed field`, `finds borrowed field`, and
  `tracks borrowed field name` are direct observations of one small helper.
  The conflict fatal-path check is the more valuable behavior test.
- `varSatMechanical totalStress` checks now use an independent effective-stress
  input instead of only converting back from the effective stress produced
  moments earlier.

### `tests/openfoam/couplingTermsUnit/Allrun`

Tested behavior:

- Builds the `couplingTermsUnit` executable.
- Creates a minimal one-cell OpenFOAM mesh.
- Seeds disk scalar fields at times `0` and `1`.
- Runs the normal unit executable path.
- Verifies expected fatal diagnostics for:
  - missing scalar disk reader start field.
  - missing `gamma_water` for head-based pressure units.
  - duplicate shared-registry field registration.

Low-value or tautological candidates:

- The diagnostic substring checks are brittle to message wording, but they are
  not tautological: they confirm that fatal paths fail for the intended reason.

## Highest-Value Coverage Today

- Golden and invariant checks for saturation laws, storage laws, effective
  stress models, and variably saturated interface coupling.
- Integration checks that shared fluid fields are registered on the solid mesh
  by identity, especially for `S` and `n`.
- Re-assembly checks proving saturation changes update explicit coupling while
  implicit stabilization remains saturation-independent.
- Fatal-path checks in `couplingTermsUnit/Allrun`.
- Fatal-path check that `varSatPoroSolid` rejects the saturated `poroFluid`
  model.
- Fatal-path checks that unassembled explicit and implicit coupling DTOs cannot
  be consumed.
- Two-zone `poroHydraulicModel` assembly from disk-backed cellZones, including
  a regression check for named zone dictionary lookup.
- Disk reread-on-time-index-change and field-registration idempotence checks.

## Suggested Cleanup Priority

1. Keep fixture smoke tests, but do not add many more unless they guard a known
   fragile initialization path.
2. Prefer independent expected values over comparing one public method's output
   with another public method on the same object.
3. When testing inverse transformations, include an independent golden value for
   at least one direction to avoid round-trip tautologies.
4. Convert pure policy string checks into comments or build/test smoke checks if
   they start creating churn.
