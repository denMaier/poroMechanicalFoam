# PoroMechanicalFoam Debug Resume

Date: 2026-04-08

## Scope

We worked on safe fixes/refactors in `poroSolidInterfaces` and related coupling code, then chased build/runtime regressions.

## What We Changed

### Coupling correctness and safety

- Fixed mapped displacement use for porosity update on non-shared meshes in:
  - `poroSolidInterfaces/poroSolid/poroSolid.C`
- Fixed boundary relative permeability call in:
  - `materialModels/poroHydraulicModel/varSatPoroHydraulicModel.C`
- Fixed Biot-law detection accumulation in:
  - `poroSolidInterfaces/poroSolidInterface.C`
- Fixed include guard collision in:
  - `poroSolidInterfaces/varSatPoroSolid/varSatPoroSolid.H`
- Added multiple fail-fast validity checks for coupling fields/mappers in:
  - `poroSolidInterfaces/poroSolidInterface.H`
  - `poroSolidInterfaces/poroSolidInterface.C`
  - `poroSolidInterfaces/poroSolid/poroSolid.C`
  - `poroSolidInterfaces/varSatPoroSolid/varSatPoroSolid.C`

### Safe parent-class consolidation

- Moved duplicated helper logic to parent class:
  - `nDot(...)`
  - `fixedStressStabil(...)`
  - pressure mapping helper for non-shared meshes
- Updated daughters to call parent helpers.

### Build-system/linking fixes attempted and applied

- Fixed compile issue for helper return type qualification in:
  - `poroSolidInterfaces/poroSolidInterface.C`
- Fixed bad solver lib search path in:
  - `solvers/initPoroMechanicalFoam/Make/options`
- Added missing minimal-library sources to resolve link symbols in:
  - `solids4FoamModelsMinimal/Make/files`
  - `solids4FoamModelsMinimal/Make/files.openfoam`
  - included `robin`, `thermalRobin`, `sphericalCavityStressDisplacement`
- Improved `S4F_ROOT` fallback and fail-fast behavior in:
  - `Allwmake`
  - `solids4FoamModelsMinimal/Allwmake`

### Cell-zone runtime issue

- Fixed rerun crash in `poroFluidModel::addDefaultCellZone()` for meshes with existing point/face zones:
  - now copies existing zones, clears zone meshes, and re-adds with appended `defaultZone`.
  - file: `poroFluidModels/poroFluidModel.C`

### Iteration residual (`delta`) runtime issues

- Fixed unallocated `vf_`/construction path in:
  - `iterationControls/iterationResidual/deltaVf/deltaVf.C`
- Added prev-iteration tracking (`prevIterStored_`) to avoid `prevIter()` abort:
  - `iterationControls/iterationResidual/deltaVf/deltaVf.H`
  - `iterationControls/iterationResidual/deltaVf/deltaVf.C`

## Runtime Error Progression

1. `unallocated autoPtr ...` in `deltaVf::makeDeltaVf()` -> fixed.
2. `previous iteration field ... not stored` -> fixed.
3. `Unknown iterationResidual type nDot` in outer coupling residual setup/check.
4. After fixing `nDot` discoverability, next blocker appeared: dimension mismatch in `deltaVf::calcResidual()` at later timestep.
5. Both issues now resolved for benchmark rerun (see update below).

## What We Tried For `nDot` Lookup

- Kept `nDot_` alive until after outer loop convergence check (moved clear to after `do ... while` loop).
- Tried explicit `checkIn()` of `nDot_` after creation.
- Tried creating `nDot` with explicit `IOobject(... registerObject=true)`.
- Latest update: create `nDot` with explicit `IOobject(register=true)` in parent helper and transfer ownership directly from `tmp` into `nDot_` via `ptr()` (avoid copy losing registration behavior).

Touched files for latest `nDot` attempt:

- `poroSolidInterfaces/poroSolidInterface.C`
- `poroSolidInterfaces/poroSolid/poroSolid.C`
- `poroSolidInterfaces/varSatPoroSolid/varSatPoroSolid.C`

## Update (after this resume)

- Implemented persistent `nDot_` handling across outer-coupling iterations in both:
  - `poroSolidInterfaces/poroSolid/poroSolid.C`
  - `poroSolidInterfaces/varSatPoroSolid/varSatPoroSolid.C`
- `nDot_` is now allocated once per timestep/loop context, updated each iteration, and explicitly checked into the mesh registry if needed.
- Added failure diagnostics in `deltaVf::firstLookup()` to print which meshes were searched and field-type presence flags.
- Reproduced run in:
  - `/home/maier/OpenFOAM/maier-v2512/run/Benchmark/01_1DColumn`
- Result:
  - `nDot` is now discovered during convergence checks (`Found nDot in solid ...`).
  - Simulation progressed but then exposed a second issue: dimension mismatch in `deltaVf::calcResidual()`.
- Added robustness for residual field dimensions in:
  - `iterationControls/iterationResidual/deltaVf/deltaVf.C`
  - synchronize `dimensions_` with the actual tracked field each lookup
  - recreate `deltaVf_` when expected dimensions change
  - keep previous-iteration storage logic consistent after reinitialization
- Rebuilt `libporoModels` and reran the same benchmark case to completion (`End`, empty `calc_debug.err`).

## What Still Needs To Be Done

1. Validate the same fixes on at least one non-shared-mesh case and one `varSatPoroSolid` case.
2. If desired, reduce/guard the added `deltaVf::firstLookup()` diagnostics once broader validation is complete.

## Notes

- Build/test execution was mostly performed by user; fixes were made conservatively and iteratively from topmost errors.
- Current state: `nDot` discoverability issue is fixed in tested benchmark; benchmark also survives previously observed residual-dimension crash.
