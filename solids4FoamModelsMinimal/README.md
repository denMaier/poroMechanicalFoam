# solids4FoamModelsMinimal

This folder provides a minimal `wmake` target that builds a reduced solids4Foam
library focused on:

- `physicsModel`
- `solidModels` (excluding `coupled*` solvers)
- `materialModels`

The source list is derived from `src/solids4FoamModels/Make/files.openfoam` and
excludes fluid, fluid-solid interface, and function object implementations.

Additional minimal-target changes:

- OpenFOAM-only build target (`files.openfoam` only; no foam-extend target).
- PETSc helper sources are removed from the minimal source list.
- PETSc and Eigen integrations are explicitly disabled with
  `-DS4F_NO_USE_PETSC -DS4F_NO_USE_EIGEN`.
- `Make/options` is reduced to the core include paths/libraries required by the
  selected solids/material/numerics/high-order sources.
- Dependencies on `RBFMeshMotionSolver` and `blockCoupledSolids4FoamTools`
  were removed for this OpenFOAM-only target; the retained sources do not
  require these libraries in the OpenFOAM code path.

Build with:

```bash
cd src/solids4FoamModelsMinimal
./Allwmake
```

The resulting library name is:

- `libsolids4FoamModelsMinimal`
