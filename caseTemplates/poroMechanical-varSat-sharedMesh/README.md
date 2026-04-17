# Coupled Variably Saturated Poromechanics Template

Use this as a starting point for:

- coupled variably saturated poromechanics
- `varSatPoroSolid` plus `varSatPoroFluid` or `varSatPoroFluidHead`
- shared-mesh workflows

This is a skeleton rather than a turnkey tutorial.

## What Is Included

- `system/` dictionaries for a coupled run
- `constant/poroFluidProperties`
- `constant/poroHydraulicProperties`
- `constant/poroCouplingProperties`
- a placeholder `constant/mechanicalProperties`
- placeholder notes for the initial fields

## What You Still Need To Do

- choose between `varSatPoroFluid` and `varSatPoroFluidHead`
- replace the placeholder `mechanicalProperties` with a valid solid setup
- replace the SWCC and effective-stress parameters
- add mesh and boundary-condition files that match your region names and patches
- tune coupling convergence and linear solvers

