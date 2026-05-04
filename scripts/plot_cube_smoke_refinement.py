#!/usr/bin/env python3
"""Plot convergence diagnostics for cube_smoke_refinement_summary.csv."""

import argparse
import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_INPUT = "data/output/cube_smoke_refinement_summary.csv"
DEFAULT_INTERNAL_REL_PNG = "data/output/cube_smoke_refinement_internal_rel.png"
DEFAULT_INTERNAL_MAX_PNG = "data/output/cube_smoke_refinement_internal_max.png"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot internal-error convergence from cube smoke refinement data."
    )
    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT,
        help=f"Input refinement summary CSV. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--internal-rel-output",
        default=DEFAULT_INTERNAL_REL_PNG,
        help=f"Output PNG for internal relative error. Default: {DEFAULT_INTERNAL_REL_PNG}",
    )
    parser.add_argument(
        "--internal-max-output",
        default=DEFAULT_INTERNAL_MAX_PNG,
        help=f"Output PNG for internal max absolute error. Default: {DEFAULT_INTERNAL_MAX_PNG}",
    )
    return parser.parse_args()


def parse_float(row, column, line_number):
    try:
        value = float(row[column])
    except KeyError as exc:
        raise ValueError(f"missing required column: {column}") from exc
    except ValueError as exc:
        raise ValueError(f"invalid {column} on line {line_number}: {row[column]}") from exc

    if not math.isfinite(value):
        raise ValueError(f"non-finite {column} on line {line_number}: {row[column]}")

    return value


def parse_int(row, column, line_number):
    try:
        value = int(row[column])
    except KeyError as exc:
        raise ValueError(f"missing required column: {column}") from exc
    except ValueError as exc:
        raise ValueError(f"invalid {column} on line {line_number}: {row[column]}") from exc

    return value


def read_rows(path):
    rows = []

    with open(path, newline="") as csv_file:
        reader = csv.DictReader(csv_file)

        for line_number, row in enumerate(reader, start=2):
            label = row.get("label", "")
            if label == "":
                raise ValueError(f"missing label on line {line_number}")

            rows.append(
                {
                    "label": label,
                    "nx": parse_int(row, "nx", line_number),
                    "internal_relative_error": parse_float(
                        row, "internal_relative_error", line_number
                    ),
                    "internal_max_abs_error": parse_float(
                        row, "internal_max_abs_error", line_number
                    ),
                    "is_int15": label.endswith("_int15"),
                }
            )

    if not rows:
        raise ValueError(f"input CSV has no data rows: {path}")

    return rows


def split_cases(rows):
    normal = sorted((row for row in rows if not row["is_int15"]), key=lambda row: row["nx"])
    int15 = sorted((row for row in rows if row["is_int15"]), key=lambda row: row["nx"])

    return normal, int15


def ensure_output_directory(path):
    directory = os.path.dirname(path)

    if directory:
        os.makedirs(directory, exist_ok=True)


def plot_metric(rows, y_column, ylabel, title, output_path):
    normal, int15 = split_cases(rows)

    plt.figure(figsize=(8, 5))

    if normal:
        plt.plot(
            [row["nx"] for row in normal],
            [row[y_column] for row in normal],
            marker="o",
            linewidth=2.0,
            label="células = intervalos nodais",
        )

    if int15:
        plt.plot(
            [row["nx"] for row in int15],
            [row[y_column] for row in int15],
            marker="s",
            linestyle="--",
            linewidth=2.0,
            label="integração 15x15x15",
        )

    plt.xlabel("nx")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, linestyle=":", linewidth=0.8)
    plt.legend()
    plt.tight_layout()

    ensure_output_directory(output_path)
    plt.savefig(output_path, dpi=160)
    plt.close()


def main():
    args = parse_args()

    try:
        rows = read_rows(args.input)
        plot_metric(
            rows,
            "internal_relative_error",
            "erro relativo discreto interno",
            "Convergência cube smoke: erro relativo interno",
            args.internal_rel_output,
        )
        plot_metric(
            rows,
            "internal_max_abs_error",
            "erro máximo absoluto interno",
            "Convergência cube smoke: erro máximo interno",
            args.internal_max_output,
        )
    except OSError as exc:
        print(f"Could not read/write plot data: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"Could not parse refinement summary: {exc}", file=sys.stderr)
        return 1

    print(f"input CSV: {args.input}")
    print(f"internal relative error plot: {args.internal_rel_output}")
    print(f"internal max absolute error plot: {args.internal_max_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
