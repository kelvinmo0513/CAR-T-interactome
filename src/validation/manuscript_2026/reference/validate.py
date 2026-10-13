#!/usr/bin/env python3
"""Independent readback and downloaded-service checks for the selected reference."""
import argparse
import csv
import json
from collections import Counter
from itertools import zip_longest
from pathlib import Path

import numpy as np
import pandas as pd

from oversample import LABEL, TARGET, require, require_external_output, sha256


def readback(original, derived, selected_cells, column_mapping):
    with selected_cells.open(newline="") as handle:
        original_map = list(csv.DictReader(handle))
    with column_mapping.open(newline="") as handle:
        mapping = list(csv.DictReader(handle))
    require(len(original_map) == 3759 and len(mapping) == 3777, "Unexpected mapping dimensions")
    require(len({m["source_cell"] for m in original_map}) == 3759, "Original cell IDs are not unique")
    require([int(m["output_column"]) for m in mapping] == list(range(1, 3778)),
            "Output columns are out of order")
    lookup = [int(m["original_column"]) for m in mapping]
    require(all(1 <= k <= 3759 for k in lookup), "Invalid original column index")
    counts = Counter()
    for index, row in zip(lookup, mapping):
        counts[index] += 1
        require(all(row[key] == original_map[index - 1][key]
                    for key in ("source_cell", "source_group", "phenotype")), "Mapping identity mismatch")
        require(int(row["copy_index"]) == counts[index], "Invalid copy index")
    targets = [i + 1 for i, m in enumerate(original_map) if m["phenotype"] == TARGET]
    require(len(targets) == 2 and all(counts[i] == 10 for i in targets), "Target copies differ")
    require(len({original_map[i - 1]["source_group"] for i in targets}) == 1, "Source-group support differs")
    expected_lookup = []
    for i in range(1, 3760):
        expected_lookup.extend([i] * (10 if i in targets else 1))
    require(lookup == expected_lookup, "Columns were reordered or background columns changed")
    genes = set()
    with original.open() as source, derived.open() as target:
        for number, pair in enumerate(zip_longest(source, target)):
            src, dst = pair
            require(src is not None and dst is not None, "Feature row counts differ")
            before = src.rstrip("\r\n").split("\t")
            after = dst.rstrip("\r\n").split("\t")
            require(len(before) == 3760 and len(after) == 3778, "Matrix dimensions differ")
            require(before[0] == after[0], f"Feature order mismatch at row {number}")
            require(after[1:] == [before[i] for i in lookup],
                    f"Expression/header text is not an exact copy at row {number}")
            if number == 0:
                require(before[1:] == [m["phenotype"] for m in original_map], "Source phenotype mismatch")
            else:
                require(before[0] not in genes and bool(before[0]), "Invalid feature identifier")
                genes.add(before[0])
                values = np.asarray(after[1:], dtype=float)
                require(np.isfinite(values).all() and (values >= 0).all(), "Nonfinite/negative expression")
    require(len(genes) == 20674 and "CCL13" not in genes, "Unexpected feature set")
    return dict(analysis_label=LABEL, passed=True, rows_verified=len(genes),
                all_headers_and_expression_strings_exact=True, background_columns_unchanged=3757,
                target_columns=20, unique_target_cells=2, copies_each=10, source_sample_labels=1,
                CCL13_in_SCT=False, finite_nonnegative=True,
                original_sha256=sha256(original), derived_sha256=sha256(derived))


def numeric_matrix(frame, name):
    require(frame.index.is_unique and frame.columns.is_unique, f"{name}: duplicate identifiers")
    values = frame.to_numpy(dtype=float)
    require(np.isfinite(values).all() and (values >= 0).all(), f"{name}: invalid values")
    return values


def service(signature, classes, profiles, mixture, scores, output):
    require_external_output(output)
    require(not output.exists(), "Output must be a new directory")
    sig = pd.read_csv(signature, sep="\t", index_col=0)
    cls = pd.read_csv(classes, sep="\t", header=None, index_col=0)
    pro = pd.read_csv(profiles, sep="\t", index_col=0)
    bulk = pd.read_csv(mixture, sep="\t", index_col=0)
    estimates = pd.read_csv(scores, sep="\t", index_col="Mixture")
    sv = numeric_matrix(sig, "Signature")
    pv = numeric_matrix(pro, "Inferred profiles")
    numeric_matrix(bulk, "Bulk mixture")
    expected = ["Myeloid DC", "B cells", "T cells", "Tregs", TARGET]
    require(list(sig.columns) == expected and list(cls.index) == expected, "Target/classes not retained")
    require(cls.isin([1, 2]).to_numpy().all() and cls.shape == (5, 25), "Invalid phenotype classes")
    require(((cls == 1).sum(axis=1) == 5).all() and ((cls == 1).sum(axis=0) == 1).all(),
            "Unexpected inferred-profile membership")
    require(pro.shape == (20674, 25) and len(sig) == 2076, "Unexpected selected job34 dimensions")
    require("CCL13" not in pro.index and "CCL13" not in sig.index, "Unexpected CCL13 feature")
    target_profiles = pv[:, (cls.loc[TARGET].to_numpy() == 1)]
    unique_profiles = np.unique(target_profiles.T, axis=0).shape[0]
    require(unique_profiles == 2, "Selected inferred targets should have only two distinct profiles")
    shared = sig.index.intersection(bulk.index)
    require(len(shared) == 1975, "Signature/bulk coverage differs from selected job35 inputs")
    require(bulk.shape[1] == 24 and len(estimates) == 24 and estimates.index.is_unique,
            "Expected 24 unique bulk samples")
    require(set(estimates.index) == set(bulk.columns), "Score/bulk sample join mismatch")
    values = numeric_matrix(estimates[expected], "Absolute scores")
    total = values.sum(axis=1)
    require((total > 0).all(), "Cannot normalize a zero component total")
    shares = values / total[:, None]
    require(np.allclose(shares.sum(axis=1), 1, rtol=0, atol=1e-12), "Normalized shares do not sum to 1")
    require(np.isfinite(estimates[["P-value", "Correlation", "RMSE", "Absolute score (sig.score)"]]).all().all(),
            "Invalid fit diagnostics")
    delta = float(np.max(np.abs(total - estimates["Absolute score (sig.score)"].to_numpy())))
    require(delta < 1e-6, "Component sum does not match absolute-score total")
    target_values = shares[:, -1]
    report = dict(analysis_label=LABEL, passed=True, target_retained=True,
                  signature_features=len(sig), feature_overlap_with_bulk=len(shared),
                  matrix_rank=int(np.linalg.matrix_rank(sv)), condition_number=float(np.linalg.cond(sv)),
                  overlap_rank=int(np.linalg.matrix_rank(sig.loc[shared].to_numpy())),
                  overlap_condition_number=float(np.linalg.cond(sig.loc[shared].to_numpy())),
                  independent_target_cells=2, target_columns_uploaded=20,
                  inferred_target_profiles=5, unique_inferred_target_profiles=int(unique_profiles),
                  n_samples=24, nonzero_target=int(np.count_nonzero(values[:, -1])),
                  all_zero_target=bool(not np.any(values[:, -1])),
                  target_share_range=[float(target_values.min()), float(target_values.max())],
                  maximum_absolute_sum_difference=delta,
                  scores_sha256=sha256(scores),
                  interpretation="Exploratory fitted-component shares; no independent biological validation.")
    require(report["matrix_rank"] == report["overlap_rank"] == 5, "Signature is rank deficient")
    output.mkdir(parents=True)
    pd.DataFrame(shares, index=estimates.index, columns=expected).to_csv(output / "normalized_shares.csv")
    (output / "service_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def check_logs(signature_log, deconvolution_log, signature_errors, deconvolution_errors):
    expected34 = ["replicates: 5", "sampling: 0.5", "fraction: 0.75", "k.max: 999",
                  "q.value: 0.01", "G.min: 300", "G.max: 500", "filter: FALSE", "QN: FALSE"]
    expected35 = ["perm: 100", "absolute: TRUE", "rmbatchSmode: TRUE", "QN: FALSE"]
    for log, expected in ((signature_log, expected34), (deconvolution_log, expected35)):
        text = log.read_text()
        require(all("[Options] " + token in text for token in expected), "Recorded service settings differ")
        require("All done." in text, "Service completion not recorded")
    require(signature_errors.stat().st_size == deconvolution_errors.stat().st_size == 0,
            "Service error log is not empty; review it before using outputs")
    return dict(analysis_label=LABEL, passed=True, signature_job=34, deconvolution_job=35,
                settings_match=True, error_logs_empty=True, execution="CIBERSORTx web GUI",
                internal_service_seed="not exposed")


def main():
    parser = argparse.ArgumentParser(description=LABEL)
    sub = parser.add_subparsers(dest="command", required=True)
    reference = sub.add_parser("reference")
    for name in ("original", "derived", "selected-cells", "column-mapping", "output"):
        reference.add_argument("--" + name, type=Path, required=True)
    remote = sub.add_parser("service")
    for name in ("signature", "classes", "profiles", "mixture", "scores", "output"):
        remote.add_argument("--" + name, type=Path, required=True)
    logs = sub.add_parser("logs")
    for name in ("signature-log", "deconvolution-log", "signature-errors", "deconvolution-errors", "output"):
        logs.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    require_external_output(args.output)
    if args.command == "service":
        result = service(args.signature, args.classes, args.profiles, args.mixture, args.scores, args.output)
    else:
        require(not args.output.exists(), "Output must be a new file")
        if args.command == "reference":
            result = readback(args.original, args.derived, args.selected_cells, args.column_mapping)
        else:
            result = check_logs(args.signature_log, args.deconvolution_log,
                                args.signature_errors, args.deconvolution_errors)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
