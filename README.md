### Organizing the source code
Please place all your sources into the `src` folder.

Binary files must not be uploaded to the repository (including executables).

Mesh files should not be uploaded to the repository. If applicable, upload `gmsh` scripts with suitable instructions to generate the meshes (and ideally a Makefile that runs those instructions). If not applicable, consider uploading the meshes to a different file sharing service, and providing a download link as part of the building and running instructions.

### Compiling
To build the executable, make sure you have loaded the needed modules with
```bash
$ module load gcc-glibc dealii
```
Then run the following commands:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```
The executable will be created into `build`, and can be executed through
```bash
$ ./executable-name
```

====================================================================================

### Wave Equation

#### Problem statement

Implement a 2D finite element solver for the wave equation:

![Wave equation problem statement](img/equation.png)

TODO: discuss the time/space discretization, the numerical dissipation and dispersion behavior, and
computational/algorithmic aspects.

#### Simple approach

The code in `src/` uses:

- Space discretization: continuous `P1` finite elements (`FE_Q<2>(1)`).
- Time discretization: second-order centered finite differences (leapfrog/central difference).
- Time stepping: at least 500 time steps over `[0, T]` (possibly more if required by the CFL estimate).
- Current default run uses `T = 5` to visualize multiple oscillation cycles in the animation.
- Linear algebra per time step: solve a mass-matrix system with a direct sparse solve.
- Boundary conditions: imposed strongly each step through deal.II boundary-value elimination.
- Test case: standing wave on unit square,
  `u(x,y,t) = sin(πx) sin(πy) cos(√2 π t)`, with homogeneous Dirichlet boundary data and zero forcing.

#### Mesh handling

Mesh binaries are intentionally not kept in the repository.

- Gmsh geometry file: `mesh/unit_square.geo`
- Mesh generation script: `scripts/generate_mesh.sh`

Generate a mesh `mesh/unit_square.msh` locally (requires `gmsh`):

```bash
./scripts/generate_mesh.sh
```

#### Running

From `build`:

```bash
cmake ..
make -j
./wave-equation
```

Optional custom mesh path:

```bash
./wave-equation /custom/path/to/mesh.msh
```

Outputs are written as `solution-XXXX.vtu` files in a `solution/` folder.

- If run from `build/`, files are written to `../solution/` (repo root).
- If run from repo root, files are written to `solution/`.
