# solids4FoamModelsMinimal integration (OpenFOAM v2412)

This repository uses a pinned solids4Foam dependency via submodule and builds a
dedicated minimal library (`libsolids4FoamModelsMinimal`).

## Supported platform

- OpenFOAM: **v2412 only**
- solids4Foam: **tag v2.3** (submodule at `external/solids4foam`)

## Setup

`Allwmake` handles submodule initialisation automatically on first build using
sparse checkout (cone mode), fetching `src/solids4FoamModels` and
`src/blockCoupledSolids4FoamTools` from tag v2.3:

```bash
./Allwmake
```

To use an existing solids4Foam source tree instead:

```bash
export S4F_ROOT=/path/to/solids4foam
./Allwmake
```

To set up the submodule manually (e.g. without running the build):

```bash
git submodule update --init external/solids4foam
git -C external/solids4foam sparse-checkout init --cone
git -C external/solids4foam sparse-checkout set src/solids4FoamModels src/blockCoupledSolids4FoamTools
git -C external/solids4foam checkout v2.3
```

If the shipped Git is too old for `git sparse-checkout`, `Allwmake` falls back
to configuring `core.sparseCheckout` and `.git/info/sparse-checkout` manually.

## Submodule scope

`.gitmodules` sets `sparseCheckout = true` for `external/solids4foam`.
The sparse-checkout cone (`src/solids4FoamModels` and
`src/blockCoupledSolids4FoamTools`) is configured by `Allwmake`
on first run and stored in `.git/modules/external/solids4foam/info/sparse-checkout`
(not tracked, but reproducibly recreated by the build script).

## What is built into `libsolids4FoamModelsMinimal`

The minimal target uses `solids4FoamModelsMinimal/Make/files` as the source
list and compiles a reduced solids4Foam models subset (physics/solid/material models
plus required numerics helpers) directly from `$(S4F_ROOT)/src/...`.

`solids4FoamModelsMinimal/Allwmake` keeps only build metadata in this repository,
prepares upstream `lnInclude` trees with `wmakeLnInclude`, and compiles selected
sources directly from:

- `external/solids4foam/src/solids4FoamModels/physicsModel`
- `external/solids4foam/src/solids4FoamModels/solidModels`
- `external/solids4foam/src/solids4FoamModels/materialModels`
- `external/solids4foam/src/solids4FoamModels/numerics`
- `external/solids4foam/src/solids4FoamModels/higherOrderHelpers`
- `external/solids4foam/src/solids4FoamModels/fluidModels/fvPatchFields`
- `external/solids4foam/src/solids4FoamModels/functionObjects/sphericalCavityAnalyticalSolution`
- `external/solids4foam/src/blockCoupledSolids4FoamTools`

## Include/link requirements

`solids4FoamModelsMinimal/Make/options` uses:

- Includes:
  - `external/solids4foam/src/solids4FoamModels/lnInclude`
  - OpenFOAM core include trees (`finiteVolume`, `meshTools`, `dynamicMesh`, `dynamicFvMesh`, `finiteArea`)
- Linked libraries:
  - `-lOpenFOAM -lfiniteVolume -lmeshTools -ldynamicFvMesh -ldynamicMesh -lfiniteArea -ltopoChangerFvMesh -lfvMotionSolvers`

`ThirdParty/eigen3` is no longer required by this repository's active build path.

Main poroMechanicalFoam library and solvers now link against `-lsolids4FoamModelsMinimal` instead of
`-lsolids4FoamModels`.

## Build flow

`./Allwmake` now does:

1. checks OpenFOAM is v2412
2. builds `poroHydraulicModels`
3. builds `solids4FoamModelsMinimal` (`solids4FoamModelsMinimal/Allwmake`)
4. builds poroMechanicalFoam libs/solvers
