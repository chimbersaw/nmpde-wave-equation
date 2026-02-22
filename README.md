### Wave Equation

#### Problem statement

Implement a 2D finite element solver for the wave equation in 2D:

<p align="center">
    <img src="img/equation.png" alt="wave equation" width="400"/>
</p>

#### Mesh handling

Mesh binaries are intentionally not kept in the repository.

- Gmsh geometry files: `mesh/*.geo`
- Mesh generation script: `scripts/generate_mesh.sh`

Generate meshes locally (requires `gmsh`):

```bash
./scripts/generate_mesh.sh
```
---

## Original repository instructions

### Organizing the source code
Please place all your sources into the `src` folder.

Binary files must not be uploaded to the repository (including executables).

Mesh files should not be uploaded to the repository. If applicable, upload `gmsh` scripts with suitable instructions to generate the meshes (and ideally a Makefile that runs those instructions). If not applicable, consider uploading the meshes to a different file sharing service, and providing a download link as part of the building and running instructions.

## Solver usage

### Build

```bash
cmake -S . -B cmake-build-release -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-release -j
```

### Run with config

```bash
mpirun -n 4 ./cmake-build-release/wave-equation --config configs/standing_wave/theta_crank_nicolson.cfg
```

### Available presets

Standing wave (`configs/standing_wave/`):
- `theta_forward.cfg`
- `theta_crank_nicolson.cfg`
- `theta_backward.cfg`
- `newmark_avg_accel.cfg`
- `newmark_central_difference.cfg`

Gaussian pulse (`configs/gaussian_pulse/`):
- `theta_forward.cfg`
- `theta_crank_nicolson.cfg`
- `theta_backward.cfg`
- `newmark_avg_accel.cfg`
- `newmark_central_difference.cfg`
- `periodic_center_absorbing.cfg` (periodic center source with absorbing boundary)

Convergence (standing wave, `configs/convergence/`):
- `convergence_space.cfg`
- `convergence_time.cfg`
- `convergence_both.cfg`

### Config keys (key=value)

Core:
- `mode = solve | convergence_space | convergence_time | convergence_both`
- `method = theta | newmark`
- `mesh_file`, `fe_degree`, `wave_speed`, `dt`, `n_steps`
- `output_interval`, `output_dir`
- `scenario_u0`, `scenario_u1`, `scenario_f`, `scenario_bc`
- optional damping profile: `scenario_sigma` (default `zero`)

Method-specific:
- `theta` for `method=theta`
- `beta`, `gamma` for `method=newmark`

Convergence-specific:
- `convergence_mesh_files`
- `convergence_dt_values`
- `convergence_reference_case`
- `convergence_csv_space`
- `convergence_csv_time`

Absorbing boundary:
- set `scenario_bc = absorbing` to activate the first-order absorbing boundary contribution.

### Convergence outputs

Convergence runs write CSV files in `results/`:

- `results/convergence_space.csv`
- `results/convergence_time.csv`

Generate plots with:

```bash
python3 scripts/plot_convergence.py --space-csv results/convergence_space.csv --time-csv results/convergence_time.csv --output results/convergence.png
```

### Output files

- `solution/`: time-series visualization files (`.vtu` + `.pvtu`) for ParaView.
- `results/`: convergence CSV files and generated plots.
