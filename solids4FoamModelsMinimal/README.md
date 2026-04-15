# solids4FoamModelsMinimal

This folder provides a minimal `wmake` target that builds a reduced solids4Foam
library focused on:

- `physicsModel`
- `solidModels` (excluding `coupled*` solvers)
- `materialModels`

The source list is derived from `src/solids4FoamModels/Make/files.openfoam`,
rewritten to point directly at `$(S4F_ROOT)/src/...`, and excludes fluid,
fluid-solid interface, and function object implementations except for the
small compatibility subset needed by retained thermal and analytical traction
patch-field code.

Additional minimal-target changes:

- OpenFOAM-only build target with direct upstream source paths (no local source-tree symlinks).
- PETSc helper sources are removed from the minimal source list.
- PETSc and Eigen integrations are explicitly disabled with
  `-DS4F_NO_USE_PETSC -DS4F_NO_USE_EIGEN`.
- `Make/options` is reduced to the core include paths/libraries required by the
  selected solids/material/numerics/high-order sources.
- Dependencies on `RBFMeshMotionSolver` and the `blockCoupledSolids4FoamTools`
  library were removed for this OpenFOAM-only target, although two helper
  sources are still compiled directly from `src/blockCoupledSolids4FoamTools`.

Build with:

```bash
cd src/solids4FoamModelsMinimal
./Allwmake
```

The resulting library name is:

- `libsolids4FoamModelsMinimal`
