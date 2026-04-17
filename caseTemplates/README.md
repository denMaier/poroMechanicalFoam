# Case Templates

This folder contains starter case folders for the supported workflows in this repository.

These templates are intended as:

- a checklist of required files
- a starting dictionary layout
- a place to begin adapting boundary conditions, materials, and geometry

They are not meant to be universal turnkey examples.

## Included Templates

- `poroFluid-saturated`
  - saturated flow-only case using `poroFluid`
  - solves for `p_rgh`
- `poroFluid-varSat-p_rgh`
  - variably saturated flow-only case using `varSatPoroFluid`
  - solves for `p_rgh`
- `poroFluid-varSat-pHead`
  - variably saturated flow-only case using `varSatPoroFluidHead`
  - solves for `pHead`
- `poroMechanical-saturated-sharedMesh`
  - coupled saturated poromechanics skeleton using `poroSolid`
  - shared-mesh setup
- `poroMechanical-varSat-sharedMesh`
  - coupled variably saturated poromechanics skeleton using `varSatPoroSolid`
  - shared-mesh setup

## Important Notes

- The fluid-only templates include a simple `blockMeshDict` and matching boundary names.
- The coupled templates are intentionally more skeletal.
  - They include the required poromechanical dictionaries.
  - You still need to adapt them to the solids4Foam solid-region setup you actually use.
- All templates use a flat `poroHydraulicProperties` with top-level entries. This works because the solver falls back to the top-level dict when no zone-specific sub-dict is found — and meshes without cell zones get an auto-created `defaultZone` that triggers this fallback.
- For multi-zone meshes, add one named sub-dict per zone (using the zone name as the key). Each sub-dict contains the same entries as the flat layout.
- Material values are placeholders. Replace them before running real cases.
- The coupled templates assume you will supply a valid `mechanicalProperties` structure for your chosen solids4Foam/solid mechanics workflow.

See [`docs/USER_GUIDE.md`](../docs/USER_GUIDE.md) for dictionary reference, available law types, and common setup pitfalls.

## Recommended Use

1. Copy the closest template to a new case directory.
2. Rename patches, fields, and regions to match your mesh.
3. Replace placeholder material values.
4. Adjust `fvSchemes`, `fvSolution`, and convergence tolerances.
5. Add any required initialization fields or preprocessing steps.

