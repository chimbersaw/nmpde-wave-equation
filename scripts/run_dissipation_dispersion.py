#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


METHODS = [
    {
        "key": "newmark_leapfrog",
        "title": r"Newmark ($\beta=0$)",
        "method": "newmark",
        "params": {"beta": 0.0, "gamma": 0.5},
    },
    {
        "key": "newmark_avg_accel",
        "title": r"Newmark ($\beta=0.25$)",
        "method": "newmark",
        "params": {"beta": 0.25, "gamma": 0.5},
    },
    {
        "key": "theta_forward",
        "title": r"$\theta$-FE ($\theta=0$)",
        "method": "theta",
        "params": {"theta": 0.0},
    },
    {
        "key": "theta_crank_nicolson",
        "title": r"$\theta$-CN ($\theta=0.5$)",
        "method": "theta",
        "params": {"theta": 0.5},
    },
    {
        "key": "theta_backward",
        "title": r"$\theta$-BE ($\theta=1$)",
        "method": "theta",
        "params": {"theta": 1.0},
    },
]

PLOT_ORDER = [
    "theta_backward",
    "theta_crank_nicolson",
    "newmark_avg_accel",
    "theta_forward",
    "newmark_leapfrog",
]


def parse_dt_values(text):
    return [float(value.strip()) for value in text.split(",") if value.strip()]


def format_dt(dt):
    return f"{dt:g}"


def safe_dt_name(dt):
    return format_dt(dt).replace(".", "p").replace("-", "m")


def write_config(path, root, method, dt, args):
    n_steps = max(1, int(round(args.final_time / dt)))
    diagnostics_csv = args.output_dir / "series" / f"{method['key']}_dt_{safe_dt_name(dt)}.csv"

    lines = [
        "mode = solve",
        f"method = {method['method']}",
        f"mesh_file = {args.mesh_file}",
        "fe_degree = 1",
        f"wave_speed = {args.wave_speed:g}",
        f"dt = {dt:.16g}",
        f"n_steps = {n_steps}",
        f"output_interval = {n_steps + 1}",
        f"output_dir = {args.output_dir / 'vtu'}",
        "write_solution = false",
        "u0 = standing_wave_5x5",
        "u1 = standing_wave_5x5_velocity",
        "f = standing_wave_5x5_forcing",
        "bc = zero_dirichlet",
        "convergence_reference_case = standing_wave_5x5_exact",
        f"diagnostics_csv = {diagnostics_csv}",
        "diagnostics_interval = 1",
        f"probe_x = {args.probe_x:g}",
        f"probe_y = {args.probe_y:g}",
        f"divergence_energy_ratio = {args.divergence_energy_ratio:g}",
    ]

    if method["method"] == "theta":
        lines.append(f"theta = {method['params']['theta']:g}")
    else:
        lines.append(f"beta = {method['params']['beta']:g}")
        lines.append(f"gamma = {method['params']['gamma']:g}")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return diagnostics_csv


def run_sweep(root, args):
    configs_dir = args.output_dir / "configs"
    configs_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "series").mkdir(parents=True, exist_ok=True)

    jobs = []
    for method in METHODS:
        for dt in args.dt_values:
            config_path = configs_dir / f"{method['key']}_dt_{safe_dt_name(dt)}.cfg"
            diagnostics_csv = write_config(config_path, root, method, dt, args)
            jobs.append((method, dt, config_path, diagnostics_csv))

    for index, (method, dt, config_path, _) in enumerate(jobs, start=1):
        print(f"[{index}/{len(jobs)}] {method['key']} dt={format_dt(dt)}", flush=True)
        subprocess.run(
            [str(args.executable), "--config", str(config_path)],
            cwd=root,
            check=True,
            timeout=args.timeout,
        )


def read_series(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def to_float(row, key):
    try:
        return float(row[key])
    except Exception:
        return math.nan


def cropped_xy(rows, x_key, y_key, cap=None, absolute_cap=False):
    xs = []
    ys = []
    for row in rows:
        x = to_float(row, x_key)
        y = to_float(row, y_key)
        if not math.isfinite(x) or not math.isfinite(y):
            break
        xs.append(x)
        ys.append(y)
        if cap is not None:
            test_value = abs(y) if absolute_cap else y
            if test_value > cap:
                ys[-1] = math.copysign(cap, y) if absolute_cap else cap
                break
    return xs, ys


def build_summary(output_dir, dt_values):
    rows = []
    for method in METHODS:
        for dt in dt_values:
            path = output_dir / "series" / f"{method['key']}_dt_{safe_dt_name(dt)}.csv"
            series = read_series(path)
            if not series:
                rows.append([method["key"], format_dt(dt), 0, "", "", "", "missing"])
                continue

            final = series[-1]
            errors = [abs(to_float(row, "probe_error")) for row in series if math.isfinite(to_float(row, "probe_error"))]
            energy_ratio = to_float(final, "energy_ratio")
            status = "ok"
            if not math.isfinite(energy_ratio):
                status = "nonfinite"
            elif energy_ratio > 1e4:
                status = "diverged"
            elif to_float(final, "time") + 1e-12 < max(to_float(row, "time") for row in series):
                status = "stopped"

            rows.append(
                [
                    method["key"],
                    format_dt(dt),
                    len(series),
                    f"{to_float(final, 'time'):.12g}",
                    f"{energy_ratio:.12g}",
                    f"{max(errors) if errors else math.nan:.12g}",
                    status,
                ]
            )

    path = output_dir / "dissipation_dispersion_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "dt", "n_rows", "final_time", "final_energy_ratio", "max_abs_probe_error", "status"])
        writer.writerows(rows)


def color_map(dt_values):
    ordered = sorted(dt_values)
    cmap = plt.get_cmap("viridis")
    if len(ordered) == 1:
        return {ordered[0]: cmap(0.5)}
    return {dt: cmap(0.08 + 0.84 * i / (len(ordered) - 1)) for i, dt in enumerate(ordered)}


def ordered_methods():
    by_key = {method["key"]: method for method in METHODS}
    return [by_key[key] for key in PLOT_ORDER]


def make_dashboard_axes(title, ylabel):
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), sharex=True)
    fig.suptitle(title, fontsize=15, fontweight="bold")
    flat_axes = axes.ravel()
    flat_axes[0].set_ylabel(ylabel)
    flat_axes[3].set_ylabel(ylabel)
    for ax in flat_axes[:5]:
        ax.set_xlabel("time")
        ax.grid(True, which="major", alpha=0.25, linewidth=0.8)
        ax.grid(True, which="minor", alpha=0.12, linewidth=0.5)
    return fig, flat_axes


def plot_energy(output_dir, dt_values):
    colors = color_map(dt_values)
    fig, axes = make_dashboard_axes("Discrete energy audit", r"$E(t)/E(0)$")

    for ax, method in zip(axes[:5], ordered_methods()):
        ax.axhline(1.0, color="#263238", linestyle="-", linewidth=0.9, alpha=0.7)
        for dt in sorted(dt_values, reverse=True):
            path = output_dir / "series" / f"{method['key']}_dt_{safe_dt_name(dt)}.csv"
            rows = read_series(path)
            xs, ys = cropped_xy(rows, "time", "energy_ratio", cap=1.25)
            if xs:
                marker_every = max(1, len(xs) // 10)
                ax.plot(
                    xs,
                    ys,
                    color=colors[dt],
                    linewidth=1.25,
                    marker="o",
                    markersize=2.2,
                    markevery=marker_every,
                    label=format_dt(dt),
                )
        ax.set_title(method["title"], fontsize=11, loc="left")
        ax.set_ylim(0.0, 1.25)

    axes[5].axis("off")
    handles, labels = axes[2].get_legend_handles_labels()
    axes[5].legend(
        handles,
        labels,
        title=r"time step $\Delta t$",
        loc="upper center",
        bbox_to_anchor=(0.5, 0.92),
        frameon=False,
        ncol=2,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_dir / "dissipation_energy.png", dpi=180)
    plt.close(fig)


def plot_probe(output_dir, dt_values):
    colors = color_map(dt_values)
    fig, axes = make_dashboard_axes("Point-probe trace against exact solution", r"$u_h(x_p,t)$")

    for ax, method in zip(axes[:5], ordered_methods()):
        exact_drawn = False
        ax.axhline(0.0, color="#263238", linewidth=0.7, alpha=0.45)
        for dt in sorted(dt_values, reverse=True):
            path = output_dir / "series" / f"{method['key']}_dt_{safe_dt_name(dt)}.csv"
            rows = read_series(path)
            if not rows:
                continue
            if not exact_drawn:
                xs_exact, ys_exact = cropped_xy(rows, "time", "exact_probe_u", cap=1.2, absolute_cap=True)
                if xs_exact:
                    ax.plot(
                        xs_exact,
                        ys_exact,
                        color="#1f2933",
                        linestyle="--",
                        linewidth=1.1,
                        alpha=0.75,
                        label="exact",
                    )
                    exact_drawn = True
            xs, ys = cropped_xy(rows, "time", "probe_u", cap=1.2, absolute_cap=True)
            if xs:
                marker_every = max(1, len(xs) // 12)
                ax.plot(
                    xs,
                    ys,
                    color=colors[dt],
                    linewidth=1.2,
                    marker=".",
                    markersize=2.5,
                    markevery=marker_every,
                    label=format_dt(dt),
                )
        ax.set_title(method["title"], fontsize=11, loc="left")
        ax.set_ylim(-1.15, 1.15)

    axes[5].axis("off")
    handles, labels = axes[2].get_legend_handles_labels()
    axes[5].legend(
        handles,
        labels,
        title=r"time step $\Delta t$",
        loc="upper center",
        bbox_to_anchor=(0.5, 0.92),
        frameon=False,
        ncol=2,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_dir / "dispersion_probe.png", dpi=180)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[1]

    parser = argparse.ArgumentParser(description="Run and plot dissipation/dispersion diagnostics.")
    parser.add_argument("--executable", type=Path, default=root / "cmake-build-release" / "wave-equation")
    parser.add_argument("--mesh-file", type=Path, default=root / "mesh" / "square_5x5_h0p25.msh")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "dissipation_dispersion")
    parser.add_argument("--dt-values", default="0.15,0.1,0.05,0.02,0.01,0.005,0.001")
    parser.add_argument("--final-time", type=float, default=5.0)
    parser.add_argument("--wave-speed", type=float, default=4.0)
    parser.add_argument("--probe-x", type=float, default=2.5)
    parser.add_argument("--probe-y", type=float, default=2.5)
    parser.add_argument("--divergence-energy-ratio", type=float, default=1e8)
    parser.add_argument("--timeout", type=float, default=70.0)
    parser.add_argument("--skip-run", action="store_true", help="Only replot existing diagnostics CSV files.")
    args = parser.parse_args()

    args.output_dir = args.output_dir.resolve()
    args.mesh_file = args.mesh_file.resolve()
    args.executable = args.executable.resolve()
    args.dt_values = parse_dt_values(args.dt_values)

    if not args.skip_run:
        run_sweep(root, args)

    build_summary(args.output_dir, args.dt_values)
    plot_energy(args.output_dir, args.dt_values)
    plot_probe(args.output_dir, args.dt_values)
    print(f"Wrote {args.output_dir / 'dissipation_energy.png'}")
    print(f"Wrote {args.output_dir / 'dispersion_probe.png'}")


if __name__ == "__main__":
    main()
