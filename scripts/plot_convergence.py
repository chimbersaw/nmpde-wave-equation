#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


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


def make_plot(rows, x_key, title, output_path):
    xs = [to_float(r, x_key) for r in rows]
    l2 = [to_float(r, "l2_error") for r in rows]
    h1 = [to_float(r, "h1_error") for r in rows]

    paired = sorted(zip(xs, l2, h1), key=lambda t: t[0], reverse=True)
    xs = [p[0] for p in paired]
    l2 = [p[1] for p in paired]
    h1 = [p[2] for p in paired]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(xs, l2, "o-", label="L2 error")
    ax.loglog(xs, h1, "s-", label="H1 error")

    ax.set_xlabel(x_key)
    ax.set_ylabel("error")
    ax.set_title(title)
    ax.grid(True, which="both", ls=":")
    ax.legend()

    for i in range(1, len(xs)):
        if xs[i] > 0 and xs[i - 1] > 0 and l2[i] > 0 and l2[i - 1] > 0:
            order = math.log(l2[i - 1] / l2[i]) / math.log(xs[i - 1] / xs[i])
            ax.annotate(f"p~{order:.2f}", (xs[i], l2[i]))

    fig.tight_layout()
    fig.savefig(output_path, dpi=160)


def main():
    parser = argparse.ArgumentParser(description="Plot convergence CSV files")
    parser.add_argument("--space-csv", default="solution/convergence_space.csv")
    parser.add_argument("--time-csv", default="solution/convergence_time.csv")
    parser.add_argument("--output", default="solution/convergence.png")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    if Path(args.space_csv).exists():
      rows = load_csv(args.space_csv)
      xs = [to_float(r, "h") for r in rows]
      l2 = [to_float(r, "l2_error") for r in rows]
      h1 = [to_float(r, "h1_error") for r in rows]
      paired = sorted(zip(xs, l2, h1), key=lambda t: t[0], reverse=True)
      xs = [p[0] for p in paired]
      l2 = [p[1] for p in paired]
      h1 = [p[2] for p in paired]
      axes[0].loglog(xs, l2, "o-", label="L2")
      axes[0].loglog(xs, h1, "s-", label="H1")
      axes[0].set_title("Spatial convergence")
      axes[0].set_xlabel("h")
      axes[0].set_ylabel("error")
      axes[0].grid(True, which="both", ls=":")
      axes[0].legend()

    if Path(args.time_csv).exists():
      rows = load_csv(args.time_csv)
      xs = [to_float(r, "dt") for r in rows]
      l2 = [to_float(r, "l2_error") for r in rows]
      h1 = [to_float(r, "h1_error") for r in rows]
      paired = sorted(zip(xs, l2, h1), key=lambda t: t[0], reverse=True)
      xs = [p[0] for p in paired]
      l2 = [p[1] for p in paired]
      h1 = [p[2] for p in paired]
      axes[1].loglog(xs, l2, "o-", label="L2")
      axes[1].loglog(xs, h1, "s-", label="H1")
      axes[1].set_title("Temporal convergence")
      axes[1].set_xlabel("dt")
      axes[1].set_ylabel("error")
      axes[1].grid(True, which="both", ls=":")
      axes[1].legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=160)


if __name__ == "__main__":
    main()
