# Variably Saturated Flow Template (`p_rgh`)

Use this as a starting point for:

- flow-only variably saturated cases
- Richards-equation runs in pressure form
- `varSatPoroFluid`

This template is single-region and uses a simple box mesh from `blockMesh`.

## Before Running

- choose the Richards linearization in `constant/poroFluidProperties`
- replace the SWCC and conductivity values in `constant/poroHydraulicProperties`
- replace the boundary conditions in `0/p_rgh`
- tune the nonlinear and linear solver settings

