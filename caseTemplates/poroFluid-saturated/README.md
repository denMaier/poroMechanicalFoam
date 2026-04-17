# Saturated Flow Template

Use this as a starting point for:

- flow-only porous-media cases
- saturated groundwater calculations
- `poroFluid` with `p_rgh` as the primary field

This template is single-region and uses a simple box mesh from `blockMesh`.

## Before Running

- adjust the mesh dimensions in `constant/polyMesh/blockMeshDict`
- replace the hydraulic material values in `constant/poroHydraulicProperties`
- replace the boundary conditions in `0/p_rgh`
- tune solver settings in `system/fvSolution`

