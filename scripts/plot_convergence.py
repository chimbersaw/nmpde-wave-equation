#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter


def load_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def to_float(row, key):
    try:
        return float(row[key])
    except Exception:
        return float("nan")


def to_int(row, key, default):
    try:
        return int(float(row[key]))
    except Exception:
        return default


def expected_time_order(row):
    method = row.get("method", "")
    if method == "theta":
        theta = to_float(row, "theta")
        return 2 if abs(theta - 0.5) < 1e-12 else 1
    if method == "newmark":
        gamma = to_float(row, "gamma")
        return 2 if abs(gamma - 0.5) < 1e-12 else 1
    return None


def method_label(row):
    method = row.get("method", "")
    if method == "theta":
        theta = to_float(row, "theta")
        if abs(theta) < 1e-12:
            return r"$\theta$-method forward Euler ($\theta=0$)"
        if abs(theta - 0.5) < 1e-12:
            return r"$\theta$-method Crank-Nicolson ($\theta=0.5$)"
        if abs(theta - 1.0) < 1e-12:
            return r"$\theta$-method backward Euler ($\theta=1$)"
        return rf"$\theta$-method ($\theta={theta:g}$)"
    if method == "newmark":
        beta = to_float(row, "beta")
        gamma = to_float(row, "gamma")
        if abs(beta) < 1e-12 and abs(gamma - 0.5) < 1e-12:
            return r"Newmark leapfrog ($\beta=0,\ \gamma=0.5$)"
        if abs(beta - 0.25) < 1e-12 and abs(gamma - 0.5) < 1e-12:
            return r"Newmark average acceleration ($\beta=0.25,\ \gamma=0.5$)"
        return rf"Newmark ($\beta={beta:g},\ \gamma={gamma:g}$)"
    return ""


def add_reference_slope(ax, xs, ys, order, label, anchor="min", scale=1.0):
    points = sorted(
        [(x, y) for x, y in zip(xs, ys) if x > 0 and y > 0],
        key=lambda t: t[0],
    )
    if len(points) < 2 or order is None:
        return

    x0, y0 = points[0] if anchor == "min" else points[-1]
    ref_x = [p[0] for p in points]
    ref_y = [scale * y0 * (x / x0) ** order for x in ref_x]
    ax.loglog(ref_x, ref_y, "k--", linewidth=1.2, alpha=0.65, label=label)


def format_tick(value, _pos):
    return f"{value:g}"


def plot_scaling(csv_path, output_path):
    rows = load_csv(csv_path)
    points = sorted(
        [(to_int(r, "n_processes", 0), to_float(r, "wall_time_s")) for r in rows],
        key=lambda t: t[0],
    )
    points = [(n, t) for n, t in points if n > 0 and t > 0]
    if not points:
        raise RuntimeError(f"No positive scaling data found in {csv_path}")

    ranks = [p[0] for p in points]
    times = [p[1] for p in points]
    speedups = [times[0] / t for t in times]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    title = method_label(rows[0]) if rows else "Strong scaling"
    fig.suptitle(title, fontsize=14)

    axes[0].plot(ranks, times, "o-", label="measured")
    for n, t in zip(ranks, times):
        axes[0].annotate(f"{t:g}s", (n, t), textcoords="offset points", xytext=(0, 7), ha="center")
    axes[0].set_title("Wall-clock time")
    axes[0].set_xlabel("MPI processes")
    axes[0].set_ylabel("time, s")
    axes[0].set_xticks(ranks)
    axes[0].grid(True, ls=":")

    axes[1].plot(ranks, speedups, "o-", label="measured")
    axes[1].plot(ranks, ranks, "k--", linewidth=1.2, alpha=0.65, label="ideal")
    axes[1].set_title("Strong scaling speedup")
    axes[1].set_xlabel("MPI processes")
    axes[1].set_ylabel("speedup vs n=1")
    axes[1].set_xticks(ranks)
    axes[1].grid(True, ls=":")
    axes[1].legend()

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output_path, dpi=160)


def main():
    parser = argparse.ArgumentParser(description="Plot convergence CSV files")
    parser.add_argument("--space-csv", default="results/convergence_newmark_avg_accel_space.csv")
    parser.add_argument("--time-csv", default="results/convergence_newmark_avg_accel_time.csv")
    parser.add_argument("--output", default="results/convergence_newmark_avg_accel.png")
    parser.add_argument("--scaling-csv", default=None)
    parser.add_argument("--scaling-output", default="results/strong_scaling_research.png")
    parser.add_argument("--scaling-only", action="store_true")
    args = parser.parse_args()

    if args.scaling_csv and args.scaling_only:
        plot_scaling(args.scaling_csv, args.scaling_output)
        return

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    title = ""

    if Path(args.space_csv).exists():
      rows = load_csv(args.space_csv)
      if rows and not title:
        title = method_label(rows[0])
      xs = [to_float(r, "h") for r in rows]
      l2 = [to_float(r, "l2_error") for r in rows]
      h1 = [to_float(r, "h1_error") for r in rows]
      fe_degree = to_int(rows[0], "fe_degree", 1) if rows else 1
      paired = sorted(zip(xs, l2, h1), key=lambda t: t[0], reverse=True)
      xs = [p[0] for p in paired]
      l2 = [p[1] for p in paired]
      h1 = [p[2] for p in paired]
      axes[0].loglog(xs, l2, "o-", label="L2")
      axes[0].loglog(xs, h1, "s-", label="H1")
      add_reference_slope(axes[0], xs, l2, fe_degree + 1, rf"$O(h^{fe_degree + 1})$", scale=0.8)
      add_reference_slope(axes[0], xs, h1, fe_degree, rf"$O(h^{fe_degree})$", scale=0.8)
      axes[0].set_title("Spatial convergence")
      axes[0].set_xlabel("h")
      axes[0].set_ylabel("error")
      axes[0].grid(True, which="both", ls=":")
      axes[0].legend()

    if Path(args.time_csv).exists():
      rows = load_csv(args.time_csv)
      if rows and not title:
        title = method_label(rows[0])
      xs = [to_float(r, "dt") for r in rows]
      l2 = [to_float(r, "l2_error") for r in rows]
      h1 = [to_float(r, "h1_error") for r in rows]
      time_order = expected_time_order(rows[0]) if rows else None
      paired = sorted(zip(xs, l2, h1), key=lambda t: t[0], reverse=True)
      xs = [p[0] for p in paired]
      l2 = [p[1] for p in paired]
      h1 = [p[2] for p in paired]
      axes[1].loglog(xs, l2, "o-", label="L2")
      axes[1].loglog(xs, h1, "s-", label="H1")
      if time_order is not None:
        add_reference_slope(axes[1], xs, l2, time_order, rf"$O(\Delta t^{time_order})$", anchor="max", scale=0.8)
      axes[1].xaxis.set_major_locator(FixedLocator(xs))
      axes[1].xaxis.set_major_formatter(FuncFormatter(format_tick))
      axes[1].set_title("Temporal convergence")
      axes[1].set_xlabel("dt")
      axes[1].set_ylabel("error")
      axes[1].grid(True, which="both", ls=":")
      axes[1].legend()

    if title:
      fig.suptitle(title, fontsize=14)
      fig.tight_layout(rect=[0, 0, 1, 0.94])
    else:
      fig.tight_layout()
    fig.savefig(output_path, dpi=160)

    if args.scaling_csv:
        plot_scaling(args.scaling_csv, args.scaling_output)


if __name__ == "__main__":
    main()
