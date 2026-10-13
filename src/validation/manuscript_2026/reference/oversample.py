#!/usr/bin/env python3
"""Exact, fixed expansion of the external two-cell AND reference; no jitter."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

LABEL = "oversampled two-cell AND reference—exploratory"
TARGET = "CCL8+CCL13+ Myeloid DC"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def require_external_output(output):
    repository = Path(__file__).resolve().parents[4]
    require(not output.resolve().is_relative_to(repository),
            "Reference outputs can contain cell/sample identifiers; choose a directory outside the repository")


def expansion(mapping, expected_cells=3759):
    require(len(mapping) == expected_cells, "Unexpected selected-cell count")
    require(all(set(("source_cell", "phenotype", "source_group")) <= set(m) for m in mapping),
            "Selection map lacks required fields")
    require(all(all(m[k] for k in ("source_cell", "phenotype", "source_group")) for m in mapping),
            "Selection map contains blank required fields")
    require(len({m["source_cell"] for m in mapping}) == len(mapping), "Source cells must be unique")
    target = [i for i, row in enumerate(mapping) if row["phenotype"] == TARGET]
    require(len(target) == 2, "Exactly two original AND cells are required")
    require(len({mapping[i]["source_group"] for i in target}) == 1,
            "Selected target source-sample support differs from the approved reference")
    return [i for i in range(len(mapping)) for _ in range(10 if i in target else 1)]


def build(reference, selected_cells, output):
    require_external_output(output)
    require(not output.exists(), "Output must be a new directory")
    with selected_cells.open(newline="") as handle:
        mapping = list(csv.DictReader(handle))
    indices = expansion(mapping)
    with reference.open() as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    require(header[1:] == [m["phenotype"] for m in mapping], "Header/map mismatch")
    require(len(indices) == 3777, "Unexpected expanded column count")
    output.mkdir(parents=True)
    destination = output / "AND_2cell_oversampled10x_exploratory.txt"
    genes = set()
    with reference.open() as source, destination.open("w", newline="") as target:
        source.readline()
        target.write("\t".join([header[0]] + [header[i + 1] for i in indices]) + "\n")
        for number, line in enumerate(source, 2):
            row = line.rstrip("\r\n").split("\t")
            require(len(row) == len(header), f"Ragged reference row {number}")
            require(row[0] not in genes and bool(row[0]), f"Duplicate/blank feature at row {number}")
            values = np.asarray(row[1:], dtype=float)
            require(np.isfinite(values).all() and (values >= 0).all(),
                    f"Invalid expression at row {number}")
            genes.add(row[0])
            target.write("\t".join([row[0]] + [row[i + 1] for i in indices]) + "\n")
    require(len(genes) == 20674 and "CCL13" not in genes, "Unexpected SCT feature set")
    counts = Counter()
    expanded = []
    for out_index, original in enumerate(indices, 1):
        counts[original] += 1
        expanded.append(dict(output_column=out_index, original_column=original + 1,
                             **mapping[original], copy_index=counts[original],
                             synthetic_duplicate=counts[original] > 1))
    with (output / "column_mapping_LOCAL_ONLY.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(expanded[0]))
        writer.writeheader()
        writer.writerows(expanded)
    report = dict(analysis_label=LABEL, passed=True, original_columns=3759,
                  exported_columns=3777, background_columns=3757, target_columns=20,
                  unique_target_cells=2, copies_per_target=10, source_groups_for_target=1,
                  genes=len(genes), CCL13_in_export=False, finite_nonnegative=True,
                  frozen_reference_sha256=sha256(reference),
                  exported_reference_sha256=sha256(destination),
                  warning="Duplication adds no independent cells, donors or biological replication.")
    (output / "reference_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=LABEL)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--selected-cells", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.reference, args.selected_cells, args.output), indent=2))
