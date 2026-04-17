
# poroMechanicalFoam
poroMechanicalFoam is an extension to OpenFOAM and solids4Foam, designed to model porous flow and deformation coupled problems in porous media. It can handle entrapped gas bubbles and partial saturation through SWCC curves and advanced mechanical constitutive models for the porous skeleton deformation.
## Features
- Integration with solids4foam for solid mechanics calculations
- Handling of partial saturation with modern computational methods
- Flexibility for adjusting to specific needs
- Close alignment with OpenFOAM's solid mechanics developments
## Prerequisites
- OpenFOAM v2412 or v2512
- solids4Foam v2.3 (git submodule in `external/solids4foam`)
## Installation
1. Ensure you have the prerequisites installed.
2. Source OpenFOAM.
3. Clone this repository (the solids4Foam submodule is initialised automatically on first build).
   Or point to an existing solids4Foam source tree instead:
```
export S4F_ROOT=/path/to/solids4foam
```
4. Navigate to the poroMechanicalFoam installation folder.
5. Run the compilation script:
```
./Allwmake
```
   On first run, `Allwmake` fetches only the required subset of solids4Foam
   (`src/solids4FoamModels` and `src/blockCoupledSolids4FoamTools`) via a
   sparse checkout of tag v2.3. If the shipped Git is too old for
   `git sparse-checkout`, the script falls back to manual sparse-checkout
   configuration.

See `docs/solids4FoamModelsMinimal.md` for details on the minimal solids4Foam subset build.
## Usage
See the solver setup in the repository case files and the
[BAW report](https://www.baw.de/content/files/forschung_entwicklung/documents/B3952.00.04.70001.pdf)
for background and example workflows.

For user-facing usage instructions, see
[`docs/USER_GUIDE.md`](docs/USER_GUIDE.md).

For repository structure, build targets, runtime architecture, and extension
points, see the maintainer reference:
[`docs/REPOSITORY_DOCUMENTATION.md`](docs/REPOSITORY_DOCUMENTATION.md).
## Contributing
We welcome contributions to poroMechanicalFoam! If you'd like to contribute, please follow these steps:
1. Fork the repository
2. Create a new branch for your feature or bug fix
3. Make your changes and commit them with clear, descriptive messages
4. Push your changes to your fork
5. Submit a pull request to the main repository
Please ensure your code adheres to the existing style and include appropriate tests if applicable.
## License
This project is licensed under the GPL 3.0 License - see the [LICENSE](LICENSE) file for details.
## Acknowledgments
This project is developed and maintained by the Bundesanstalt für Wasserbau (BAW). We are grateful for their support and contributions to this work.
