# solids4FoamModelsMinimal integration (OpenFOAM v2412)

This repository uses a pinned solids4Foam dependency via submodule and builds a
dedicated minimal library (`libsolids4FoamModelsMinimal`).

## Supported platform

- OpenFOAM: **v2412 only**
- solids4Foam: **tag v2.3** (submodule at `external/solids4foam`)

## Setup

```bash
git submodule add https://github.com/solids4foam/solids4foam.git external/solids4foam
git -C external/solids4foam fetch --tags
git -C external/solids4foam checkout v2.3
```

Or point to an existing solids4Foam source tree:

```bash
export S4F_ROOT=/path/to/solids4foam
./Allwmake
```

## Submodule scope

Git submodules point to a repository commit, not a single subdirectory path.
So the submodule itself is still `external/solids4foam` (full repo identity).

If you want a smaller working tree, you can use sparse-checkout inside the
submodule after init/update and keep only the paths needed here:

- `src/solids4FoamModels`

## What is built into `libsolids4FoamModelsMinimal`

The minimal target uses `solids4FoamModelsMinimal/Make/files.openfoam` as the source
list and compiles a reduced solids4Foam models subset (physics/solid/material models
plus required numerics helpers).

`solids4FoamModelsMinimal/Allwmake` keeps only build metadata in this repository and
links the required source trees from:

- `external/solids4foam/src/solids4FoamModels/physicsModel`
- `external/solids4foam/src/solids4FoamModels/solidModels`
- `external/solids4foam/src/solids4FoamModels/materialModels`
- `external/solids4foam/src/solids4FoamModels/numerics`
- `external/solids4foam/src/solids4FoamModels/higherOrderHelpers`

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
