### Wave Equation

#### Problem statement

Implement a 2D finite element solver for the wave equation in 2D:

<p align="center">
    <img src="img/equation.png" alt="wave equation" width="400"/>
</p>

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
- Setup `standing`: standing wave on unit square,
  `u(x,y,t) = sin(πx) sin(πy) cos(√2 π t)`, with homogeneous Dirichlet boundary data and zero forcing.
- Setup `radial`: Gaussian radial pulse with reflective walls.
- Setup `radial_absorbing`: center-driven radial wave with a damping sponge layer near the boundary to reduce reflections.

#### Mesh handling

Mesh binaries are intentionally not kept in the repository.

- Gmsh geometry files: `mesh/unit_square.geo`, `mesh/radial_square.geo`
- Mesh generation script: `scripts/generate_mesh.sh`

Generate meshes locally (requires `gmsh`):

```bash
./scripts/generate_mesh.sh
```

This generates:
- `mesh/unit_square.msh`
- `mesh/radial_square.msh`

#### Running

From `build`:

```bash
cmake ..
make -j
./wave-equation --config ../configs/standing.cfg
```

Run radial pulse setup (PDF-like ring propagation):

```bash
./wave-equation --config ../configs/radial.cfg
```

Run radial pulse with absorbing wall:

```bash
./wave-equation --config ../configs/radial_absorbing.cfg
```

Config files included in the repository:

- `configs/standing.cfg`
- `configs/radial.cfg`
- `configs/radial_absorbing.cfg`

Config file format (`key = value`):

```bash
setup = standing|radial|radial_absorbing
refine = <int>
final_time = <double>
cfl = <double>
minimum_time_steps = <int>
output_every = <int>
source_amplitude = <double>
source_frequency = <double>
mesh_path = <path or empty>
```

`source_amplitude` and `source_frequency` are used by `radial_absorbing` to control
the emitted center wave strength and oscillation rate.

#### Makefile shortcuts

From repository root:

```bash
make mesh
make build
make run-standing
make run-radial
make run-radial-absorbing
make run CONFIG=configs/radial_absorbing.cfg
make clean-solution
```

To add more setups, add a new config file and corresponding setup logic in `src/`.

Show CLI help:

```bash
./wave-equation --help
```

Outputs are written as `solution-XXXX.vtu` files in a `solution/` folder.

- If run from `build/`, files are written to `../solution/` (repo root).
- If run from repo root, files are written to `solution/`.

---

## Original repository instructions

### Organizing the source code
Please place all your sources into the `src` folder.

Binary files must not be uploaded to the repository (including executables).

Mesh files should not be uploaded to the repository. If applicable, upload `gmsh` scripts with suitable instructions to generate the meshes (and ideally a Makefile that runs those instructions). If not applicable, consider uploading the meshes to a different file sharing service, and providing a download link as part of the building and running instructions.
