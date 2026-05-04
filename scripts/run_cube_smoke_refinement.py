#!/usr/bin/env python3
"""Run reproduce_cube_smoke for several refinement configurations.

This script is intentionally standard-library only. It orchestrates the
existing C executable and the existing CSV analyzer; it does not change any C
code or CMake configuration.
"""

import csv
import os
import subprocess
import sys
from pathlib import Path


L = 10.0
SUMMARY_CSV = Path("data/output/cube_smoke_refinement_summary.csv")


CASES = [
    {
        "label": "3x3x3",
        "nx": 3,
        "ny": 3,
        "nz": 3,
        "nx_cells": 2,
        "ny_cells": 2,
        "nz_cells": 2,
        "sample_nx": 5,
        "sample_ny": 5,
        "sample_nz": 5,
        "output": "data/output/cube_smoke_3x3x3.csv",
    },
    {
        "label": "5x5x5",
        "nx": 5,
        "ny": 5,
        "nz": 5,
        "nx_cells": 4,
        "ny_cells": 4,
        "nz_cells": 4,
        "sample_nx": 9,
        "sample_ny": 9,
        "sample_nz": 9,
        "output": "data/output/cube_smoke_5x5x5.csv",
    },
    {
        "label": "7x7x7",
        "nx": 7,
        "ny": 7,
        "nz": 7,
        "nx_cells": 6,
        "ny_cells": 6,
        "nz_cells": 6,
        "sample_nx": 11,
        "sample_ny": 11,
        "sample_nz": 11,
        "output": "data/output/cube_smoke_7x7x7.csv",
    },
    {
        "label": "9x9x9",
        "nx": 9,
        "ny": 9,
        "nz": 9,
        "nx_cells": 8,
        "ny_cells": 8,
        "nz_cells": 8,
        "sample_nx": 13,
        "sample_ny": 13,
        "sample_nz": 13,
        "output": "data/output/cube_smoke_9x9x9.csv",
    },
    {
        "label": "7x7x7_int15",
        "nx": 7,
        "ny": 7,
        "nz": 7,
        "nx_cells": 15,
        "ny_cells": 15,
        "nz_cells": 15,
        "sample_nx": 11,
        "sample_ny": 11,
        "sample_nz": 11,
        "output": "data/output/cube_smoke_7x7x7_int15.csv",
    },
    {
        "label": "9x9x9_int15",
        "nx": 9,
        "ny": 9,
        "nz": 9,
        "nx_cells": 15,
        "ny_cells": 15,
        "nz_cells": 15,
        "sample_nx": 13,
        "sample_ny": 13,
        "sample_nz": 13,
        "output": "data/output/cube_smoke_9x9x9_int15.csv",
    },
]

# Heavy optional cases for later manual experiments:
# {
#     "label": "11x11x11",
#     "nx": 11, "ny": 11, "nz": 11,
#     "nx_cells": 10, "ny_cells": 10, "nz_cells": 10,
#     "sample_nx": 15, "sample_ny": 15, "sample_nz": 15,
#     "output": "data/output/cube_smoke_11x11x11.csv",
# }
# {
#     "label": "15x15x15",
#     "nx": 15, "ny": 15, "nz": 15,
#     "nx_cells": 15, "ny_cells": 15, "nz_cells": 15,
#     "sample_nx": 21, "sample_ny": 21, "sample_nz": 21,
#     "output": "data/output/cube_smoke_15x15x15.csv",
# }


def run_command(command, repo_root):
    completed = subprocess.run(
        command,
        cwd=repo_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    if completed.returncode != 0:
        print(f"Command failed: {' '.join(str(part) for part in command)}")
        if completed.stdout:
            print(completed.stdout)
        if completed.stderr:
            print(completed.stderr, file=sys.stderr)
        raise RuntimeError(f"command returned {completed.returncode}")

    return completed


def read_analyzer_summary(path):
    rows_by_scope = {}

    with path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)

        for row in reader:
            rows_by_scope[row["scope"]] = row

    for required_scope in ("global", "internal"):
        if required_scope not in rows_by_scope:
            raise ValueError(f"summary {path} is missing scope {required_scope}")

    return rows_by_scope


def case_summary_path(case):
    return Path("data/output") / f"{case['label']}_error_summary.csv"


def run_case(case, repo_root, executable, analyzer):
    output_path = Path(case["output"])
    summary_path = case_summary_path(case)

    print(f"Running {case['label']} -> {output_path}")

    app_command = [
        str(executable),
        str(case["nx"]),
        str(case["ny"]),
        str(case["nz"]),
        str(case["nx_cells"]),
        str(case["ny_cells"]),
        str(case["nz_cells"]),
        str(case["sample_nx"]),
        str(case["sample_ny"]),
        str(case["sample_nz"]),
        str(output_path),
    ]
    app_result = run_command(app_command, repo_root)

    if app_result.stdout:
        print(app_result.stdout.strip())

    analyzer_command = [
        sys.executable,
        str(analyzer),
        str(L),
        "--input",
        str(output_path),
        "--output",
        str(summary_path),
    ]
    analyzer_result = run_command(analyzer_command, repo_root)

    if analyzer_result.stdout:
        print(analyzer_result.stdout.strip())

    summary = read_analyzer_summary(repo_root / summary_path)
    global_row = summary["global"]
    internal_row = summary["internal"]

    return {
        "label": case["label"],
        "nx": case["nx"],
        "ny": case["ny"],
        "nz": case["nz"],
        "nodes": case["nx"] * case["ny"] * case["nz"],
        "nx_cells": case["nx_cells"],
        "ny_cells": case["ny_cells"],
        "nz_cells": case["nz_cells"],
        "integration_cells": case["nx_cells"] * case["ny_cells"] * case["nz_cells"],
        "sample_points": case["sample_nx"] * case["sample_ny"] * case["sample_nz"],
        "global_max_abs_error": global_row["max_abs_error"],
        "global_relative_error": global_row["relative_discrete_error"],
        "internal_max_abs_error": internal_row["max_abs_error"],
        "internal_relative_error": internal_row["relative_discrete_error"],
    }


def write_refinement_summary(repo_root, rows):
    output_path = repo_root / SUMMARY_CSV
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "label",
        "nx",
        "ny",
        "nz",
        "nodes",
        "nx_cells",
        "ny_cells",
        "nz_cells",
        "integration_cells",
        "sample_points",
        "global_max_abs_error",
        "global_relative_error",
        "internal_max_abs_error",
        "internal_relative_error",
    ]

    with output_path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerows(rows)

    return output_path


def print_table(rows):
    headers = [
        "label",
        "nodes",
        "cells",
        "samples",
        "global_rel",
        "internal_rel",
        "internal_max",
    ]
    widths = {
        "label": 16,
        "nodes": 7,
        "cells": 8,
        "samples": 8,
        "global_rel": 14,
        "internal_rel": 14,
        "internal_max": 14,
    }

    print("\nRefinement summary")
    print(
        " ".join(
            header.ljust(widths[header]) if header == "label" else header.rjust(widths[header])
            for header in headers
        )
    )
    print(
        " ".join(
            ("-" * widths[header]) for header in headers
        )
    )

    for row in rows:
        values = {
            "label": row["label"],
            "nodes": str(row["nodes"]),
            "cells": str(row["integration_cells"]),
            "samples": str(row["sample_points"]),
            "global_rel": row["global_relative_error"],
            "internal_rel": row["internal_relative_error"],
            "internal_max": row["internal_max_abs_error"],
        }
        print(
            " ".join(
                values[header].ljust(widths[header])
                if header == "label"
                else values[header].rjust(widths[header])
                for header in headers
            )
        )


def main():
    repo_root = Path(__file__).resolve().parents[1]
    executable = repo_root / "build" / "reproduce_cube_smoke"
    analyzer = repo_root / "scripts" / "analyze_cube_smoke.py"

    if not executable.exists() or not os.access(executable, os.X_OK):
        print(f"Executable not found or not executable: {executable}", file=sys.stderr)
        print("Run cmake --build build first.", file=sys.stderr)
        return 1

    if not analyzer.exists():
        print(f"Analyzer script not found: {analyzer}", file=sys.stderr)
        return 1

    rows = []

    try:
        for case in CASES:
            rows.append(run_case(case, repo_root, executable, analyzer))
    except (RuntimeError, OSError, ValueError) as exc:
        print(f"Refinement run failed: {exc}", file=sys.stderr)
        return 1

    output_path = write_refinement_summary(repo_root, rows)
    print_table(rows)
    print(f"\nsummary CSV: {output_path.relative_to(repo_root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
