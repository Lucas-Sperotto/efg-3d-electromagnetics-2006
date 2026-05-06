#!/usr/bin/env python3
"""Plot the x = 5.33 cube-plane CSV exported by reproduce_cube_sparse."""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize


DEFAULT_INPUT = Path("data/output/cube_plane_x_5_33_refine15.csv")
DEFAULT_OUTPUT_DIR = Path("data/output/figures")
DEFAULT_L = 10.0
DEFAULT_LEVELS = 40
DEFAULT_DPI = 220

REQUIRED_COLUMNS = ["x", "y", "z", "V_num", "V_exact", "abs_error"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate contour plots for the cube potential plane x = 5.33."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input plane CSV. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for generated figures and metrics. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--L",
        type=float,
        default=DEFAULT_L,
        help=f"Cube side length used for interior-plane metrics. Default: {DEFAULT_L}",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=DEFAULT_LEVELS,
        help=f"Number of filled contour levels. Default: {DEFAULT_LEVELS}",
    )
    return parser.parse_args()


def read_plane_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"input CSV does not exist: {path}")

    data = pd.read_csv(path)
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in data.columns]

    if missing_columns:
        missing = ", ".join(missing_columns)
        raise ValueError(f"input CSV is missing required columns: {missing}")

    for column in REQUIRED_COLUMNS:
        data[column] = pd.to_numeric(data[column], errors="raise")

    if not np.isfinite(data[REQUIRED_COLUMNS].to_numpy()).all():
        raise ValueError("input CSV contains non-finite values in required columns")

    if data.empty:
        raise ValueError(f"input CSV has no data rows: {path}")

    duplicate_count = int(data.duplicated(subset=["y", "z"]).sum())
    if duplicate_count > 0:
        raise ValueError(f"input CSV contains duplicate y,z samples: {duplicate_count}")

    return data


def make_grid(data, value_column):
    grid = data.pivot(index="z", columns="y", values=value_column)
    grid = grid.sort_index(axis=0).sort_index(axis=1)

    missing_count = int(grid.isna().sum().sum())
    if missing_count > 0:
        raise ValueError(f"input CSV does not form a complete y-z grid: {missing_count} gaps")

    y_values = grid.columns.to_numpy(dtype=float)
    z_values = grid.index.to_numpy(dtype=float)
    y_grid, z_grid = np.meshgrid(y_values, z_values)
    values = grid.to_numpy(dtype=float)

    return y_grid, z_grid, values


def relative_discrete_error(v_num, v_exact):
    error = v_num - v_exact
    numerator = float(np.sum(error * error))
    denominator = float(np.sum(v_exact * v_exact))

    if denominator == 0.0:
        return np.nan

    return float(np.sqrt(numerator / denominator))


def compute_metrics(data, L):
    if L <= 0.0:
        raise ValueError("L must be positive")

    y = data["y"].to_numpy(dtype=float)
    z = data["z"].to_numpy(dtype=float)
    v_num = data["V_num"].to_numpy(dtype=float)
    v_exact = data["V_exact"].to_numpy(dtype=float)
    abs_error = data["abs_error"].to_numpy(dtype=float)

    tolerance = max(1.0e-12, abs(L) * 1.0e-12)
    interior = (
        (y > tolerance)
        & (y < L - tolerance)
        & (z > tolerance)
        & (z < L - tolerance)
    )

    if not interior.any():
        raise ValueError("interior-plane mask is empty")

    return {
        "points": int(len(data)),
        "x_plane_min": float(data["x"].min()),
        "x_plane_max": float(data["x"].max()),
        "y_min": float(data["y"].min()),
        "y_max": float(data["y"].max()),
        "z_min": float(data["z"].min()),
        "z_max": float(data["z"].max()),
        "V_num_min": float(np.min(v_num)),
        "V_num_max": float(np.max(v_num)),
        "V_exact_min": float(np.min(v_exact)),
        "V_exact_max": float(np.max(v_exact)),
        "plane_max_abs_error": float(np.max(abs_error)),
        "plane_mean_abs_error": float(np.mean(abs_error)),
        "plane_relative_error": relative_discrete_error(v_num, v_exact),
        "interior_points": int(np.count_nonzero(interior)),
        "interior_max_abs_error": float(np.max(abs_error[interior])),
        "interior_mean_abs_error": float(np.mean(abs_error[interior])),
        "interior_relative_error": relative_discrete_error(
            v_num[interior], v_exact[interior]
        ),
    }


def format_float(value):
    if np.isnan(value):
        return "nan"

    return f"{value:.17g}"


def metrics_text(input_path, output_dir, grid_shape, metrics):
    lines = [
        "cube_plane_x_5_33_metrics",
        f"input_csv: {input_path}",
        f"output_dir: {output_dir}",
        f"points: {metrics['points']}",
        f"grid_y_count: {grid_shape[1]}",
        f"grid_z_count: {grid_shape[0]}",
        f"x_plane_min: {format_float(metrics['x_plane_min'])}",
        f"x_plane_max: {format_float(metrics['x_plane_max'])}",
        f"y_min: {format_float(metrics['y_min'])}",
        f"y_max: {format_float(metrics['y_max'])}",
        f"z_min: {format_float(metrics['z_min'])}",
        f"z_max: {format_float(metrics['z_max'])}",
        f"V_num_min: {format_float(metrics['V_num_min'])}",
        f"V_num_max: {format_float(metrics['V_num_max'])}",
        f"V_exact_min: {format_float(metrics['V_exact_min'])}",
        f"V_exact_max: {format_float(metrics['V_exact_max'])}",
        f"plane_max_abs_error: {format_float(metrics['plane_max_abs_error'])}",
        f"plane_mean_abs_error: {format_float(metrics['plane_mean_abs_error'])}",
        f"plane_relative_error: {format_float(metrics['plane_relative_error'])}",
        f"interior_points: {metrics['interior_points']}",
        f"interior_max_abs_error: {format_float(metrics['interior_max_abs_error'])}",
        f"interior_mean_abs_error: {format_float(metrics['interior_mean_abs_error'])}",
        f"interior_relative_error: {format_float(metrics['interior_relative_error'])}",
    ]

    return "\n".join(lines) + "\n"


def save_metrics(path, text):
    path.write_text(text, encoding="utf-8")


def save_metrics_csv(path, metrics):
    rows = []

    for key, value in metrics.items():
        if isinstance(value, (int, np.integer)):
            formatted_value = str(int(value))
        else:
            formatted_value = format_float(float(value))

        rows.append({"metric": key, "value": formatted_value})

    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(path, index=False)


def plot_single_contour(
    y_grid,
    z_grid,
    values,
    title,
    colorbar_label,
    output_path,
    levels,
    cmap="viridis",
):
    fig, axis = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
    contour = axis.contourf(y_grid, z_grid, values, levels=levels, cmap=cmap)

    axis.set_xlabel("y")
    axis.set_ylabel("z")
    axis.set_title(title)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(float(np.min(y_grid)), float(np.max(y_grid)))
    axis.set_ylim(float(np.min(z_grid)), float(np.max(z_grid)))
    fig.colorbar(contour, ax=axis, label=colorbar_label)

    fig.savefig(output_path, dpi=DEFAULT_DPI)
    plt.close(fig)


def plot_comparison(
    y_grid, z_grid, v_num, v_exact, abs_error, output_path, levels, plane_x
):
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7), constrained_layout=True)

    potential_norm = Normalize(
        vmin=float(min(np.min(v_num), np.min(v_exact))),
        vmax=float(max(np.max(v_num), np.max(v_exact))),
    )
    potential_levels = np.linspace(potential_norm.vmin, potential_norm.vmax, levels)

    error_norm = Normalize(
        vmin=float(np.min(abs_error)),
        vmax=float(np.max(abs_error)),
    )
    error_levels = np.linspace(error_norm.vmin, error_norm.vmax, levels)

    panels = [
        (axes[0], v_num, "V_num", "viridis", potential_norm, potential_levels, "V"),
        (axes[1], v_exact, "V_exact", "viridis", potential_norm, potential_levels, "V"),
        (
            axes[2],
            abs_error,
            "abs_error",
            "magma",
            error_norm,
            error_levels,
            "|V_num - V_exact|",
        ),
    ]

    for axis, values, title, cmap, norm, panel_levels, colorbar_label in panels:
        contour = axis.contourf(
            y_grid,
            z_grid,
            values,
            levels=panel_levels,
            cmap=cmap,
            norm=norm,
        )
        axis.set_xlabel("y")
        axis.set_ylabel("z")
        axis.set_title(title)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(float(np.min(y_grid)), float(np.max(y_grid)))
        axis.set_ylim(float(np.min(z_grid)), float(np.max(z_grid)))
        fig.colorbar(contour, ax=axis, label=colorbar_label)

    fig.suptitle(f"Cube potential contours on x = {plane_x:.2f}")
    fig.savefig(output_path, dpi=DEFAULT_DPI)
    plt.close(fig)


def main():
    args = parse_args()

    if args.levels < 2:
        print("levels must be at least 2.", file=sys.stderr)
        return 1

    try:
        data = read_plane_csv(args.input)
        y_grid, z_grid, v_num = make_grid(data, "V_num")
        _, _, v_exact = make_grid(data, "V_exact")
        _, _, abs_error = make_grid(data, "abs_error")
        metrics = compute_metrics(data, args.L)

        args.output_dir.mkdir(parents=True, exist_ok=True)

        v_num_path = args.output_dir / "cube_plane_x_5_33_V_num_contour.png"
        v_exact_path = args.output_dir / "cube_plane_x_5_33_V_exact_contour.png"
        abs_error_path = args.output_dir / "cube_plane_x_5_33_abs_error_contour.png"
        comparison_path = args.output_dir / "cube_plane_x_5_33_comparison.png"
        metrics_path = args.output_dir / "cube_plane_x_5_33_metrics.txt"
        metrics_csv_path = args.output_dir / "cube_plane_x_5_33_metrics.csv"

        plane_x = 0.5 * (metrics["x_plane_min"] + metrics["x_plane_max"])

        plot_single_contour(
            y_grid,
            z_grid,
            v_num,
            f"Numerical potential V_num on x = {plane_x:.2f}",
            "V_num",
            v_num_path,
            args.levels,
        )
        plot_single_contour(
            y_grid,
            z_grid,
            v_exact,
            f"Analytical potential V_exact on x = {plane_x:.2f}",
            "V_exact",
            v_exact_path,
            args.levels,
        )
        plot_single_contour(
            y_grid,
            z_grid,
            abs_error,
            f"Absolute error on x = {plane_x:.2f}",
            "abs_error",
            abs_error_path,
            args.levels,
            cmap="magma",
        )
        plot_comparison(
            y_grid,
            z_grid,
            v_num,
            v_exact,
            abs_error,
            comparison_path,
            args.levels,
            plane_x,
        )

        text = metrics_text(args.input, args.output_dir, v_num.shape, metrics)
        save_metrics(metrics_path, text)
        save_metrics_csv(metrics_csv_path, metrics)
    except (OSError, ValueError, FileNotFoundError) as exc:
        print(f"Could not plot cube plane: {exc}", file=sys.stderr)
        return 1

    print(f"input CSV: {args.input}")
    print(f"points: {metrics['points']}")
    print(f"grid shape: {v_num.shape[1]} y samples x {v_num.shape[0]} z samples")
    print(f"V_num contour: {v_num_path}")
    print(f"V_exact contour: {v_exact_path}")
    print(f"abs_error contour: {abs_error_path}")
    print(f"comparison: {comparison_path}")
    print(f"metrics: {metrics_path}")
    print(f"metrics CSV: {metrics_csv_path}")
    print(f"plane max abs error: {format_float(metrics['plane_max_abs_error'])}")
    print(f"plane mean abs error: {format_float(metrics['plane_mean_abs_error'])}")
    print(f"plane relative error: {format_float(metrics['plane_relative_error'])}")
    print(f"interior max abs error: {format_float(metrics['interior_max_abs_error'])}")
    print(f"interior mean abs error: {format_float(metrics['interior_mean_abs_error'])}")
    print(f"interior relative error: {format_float(metrics['interior_relative_error'])}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
