# Coupled Saturated Poromechanics Template

Use this as a starting point for:

- coupled saturated poromechanics
- `poroSolid` plus `poroFluid`
- shared-mesh workflows

This is a skeleton rather than a turnkey tutorial.

## What Is Included

- `system/` dictionaries for a coupled run
- `constant/poroFluidProperties`
- `constant/poroHydraulicProperties`
- `constant/poroCouplingProperties`
- a placeholder `constant/mechanicalProperties`
- placeholder fields in `0/`

## What You Still Need To Do

- adapt the solid-region layout to the solids4Foam workflow you use
- replace the placeholder `mechanicalProperties` with a valid solid setup
- add mesh and boundary-condition files that match your region names and patches
- tune coupling convergence and linear solvers

