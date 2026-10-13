# Cologne deconvolution: selected manuscript workflow

The selected Cologne component comes from **CIBERSORTx web GUI signature job 34
and deconvolution job 35**. No API or container was used for these two jobs. The
code here recreates the upload reference and imports the recorded service
downloads. It does not submit jobs or require service credentials.

All Cologne reference outputs should be described as **“oversampled two-cell AND
reference—exploratory.”** The signature was derived from only two distinct
CCL8-positive **and** CCL13-positive cells, with ten exact copies of each. This
adds no independent cells, donors, or biological replication. Both cells share
one source-sample label. Clinical model confidence intervals do not incorporate
this reference uncertainty.

## External inputs and exact reference preparation

Obtain the annotated EcoTyper Seurat object and fixed selected-cell map from the
study's reproducibility archive. Their hashes and source references belong to
the input manifest; neither the object nor cell/sample mappings are distributed
in this public repository. Keep generated mappings and patient-level output
outside the checkout. The manifest's `reference` group contains these optional
reference-rebuild and service-audit inputs; its `analysis` group contains the
inputs needed for manuscript numerical reproduction. Reference commands must
write to an external directory: their mapping and normalized-share outputs are
local analysis material and must not be committed to the public repository.

The object already contains `predicted_celltype_MI` and `orig.ident` annotations.
The reference merges NK cells into T cells and Monocytes/cDC and pDC into Myeloid
DC. The target gate uses raw RNA counts `CCL8 > 0 AND CCL13 > 0`, restricted to
Myeloid DC. This gating rule is separate from the expression matrix supplied to
CIBERSORTx: the upload uses Kelvin's **SCT data** layer. **CCL13 is absent from
that SCT feature set**, although it is present in the raw-RNA gate.

From the repository root, with external absolute paths supplied by the user:

```sh
Rscript src/validation/manuscript_2026/reference/prepare_reference.R \
  --rds /external/reference/annotated_ecotyper.rds \
  --selected-cells /external/reference/selected_cells_LOCAL_ONLY.csv \
  --output /external/rebuild/pre_reference

python src/validation/manuscript_2026/reference/oversample.py \
  --reference /external/rebuild/pre_reference/Kelvin_DP_SCT_reference.txt \
  --selected-cells /external/reference/selected_cells_LOCAL_ONLY.csv \
  --output /external/rebuild/oversampled_reference

python src/validation/manuscript_2026/reference/validate.py reference \
  --original /external/rebuild/pre_reference/Kelvin_DP_SCT_reference.txt \
  --derived /external/rebuild/oversampled_reference/AND_2cell_oversampled10x_exploratory.txt \
  --selected-cells /external/reference/selected_cells_LOCAL_ONLY.csv \
  --column-mapping /external/rebuild/oversampled_reference/column_mapping_LOCAL_ONLY.csv \
  --output /external/rebuild/reference_readback.json
```

Each output directory/file must be new. The map, rather than a fresh random
sample, fixes the original selected cells and their order. The September 30
reference selection used seed 20260930; Kelvin's earlier sampling seed is not
known. The selected object has 20,674 SCT features and 3,759 selected cells.
Replacing each of the two target columns in place with ten copies gives 3,777
columns: 3,757 unchanged background cells and 20 target columns. Features,
background order and every expression string are checked independently on
readback. There is no expression jitter or additional normalization. The
historical preparation included a negative-value clamp, but the selected SCT
matrix has no negatives; the current code rejects a changed input instead of
silently altering it.

The selected pre-reference SHA-256 is
`2cad4990b056751cfa56857164e47b4886ed00fd28b06651d75d41a1f925c06d`.
The actual oversampled upload SHA-256 is
`db37330fe9082caac34e4a4143d3d1139fa637ed0cb5d5a64c232bf1da7b88f3`.

## Recorded GUI execution

| Stage | Selected settings |
|---|---|
| Signature job 34 | Five replicates; sampling 0.5; minimum expression 0.75; kappa 999; q-value 0.01; 300–500 genes per class; quantile normalization off; non-hematopoietic filtering off. |
| Deconvolution job 35 | Unchanged 24-sample Cologne/Stahl CPM matrix; job 34 signature; absolute scoring; S-mode with the same oversampled reference; quantile normalization off; 100 permutations. |

The job 34 runtime log records `[Options] fraction: 0.75`; an expanded website
summary rounded this field to 1. The recorded submission and runtime settings
are authoritative. The target was required to appear in both the downloaded
signature and inferred phenotype classes before job 35 was submitted. No
automatic retries, altered duplication counts or setting search are performed
by this module. CIBERSORTx does not expose its internal random seed, so a future
GUI rerun may differ. Exact reproduction of the manuscript's downstream result
uses the archived downloads, rather than assuming deterministic service output.

Required downloaded outputs are the job 34 signature (`*.bm.K999.txt`), phenotype
classes (`*inferred_phenoclasses.txt`), inferred profiles
(`*inferred_refsample.txt`), and job 35 `CIBERSORTx_Job35_Adjusted.txt`. Preserve
both jobs' runtime and empty error logs externally as the execution record.

```sh
python src/validation/manuscript_2026/reference/validate.py service \
  --signature /external/cologne/job34_signature.tsv \
  --classes /external/cologne/job34_classes.tsv \
  --profiles /external/cologne/job34_profiles.tsv \
  --mixture /external/cologne/Stahl_24_CPM.tsv \
  --scores /external/cologne/job35_scores.tsv \
  --output /external/rebuild/service_validation

python src/validation/manuscript_2026/reference/validate.py logs \
  --signature-log /external/cologne/job34_runtime_log.txt \
  --deconvolution-log /external/cologne/job35_runtime_log.txt \
  --signature-errors /external/cologne/job34_error_log.txt \
  --deconvolution-errors /external/cologne/job35_error_log.txt \
  --output /external/rebuild/service_log_validation.json
```

## Import, checks and interpretation

The retained components are Myeloid DC, B cells, T cells, Tregs, and CCL8+CCL13+
Myeloid DC. The imported target share is its absolute score divided by the sum
of **all five retained component scores**. It is an exploratory share of fitted
components, not an independently validated whole-tissue cell fraction. This
normalization must not be restricted to the myeloid components.

Validation checks 24 unique sample joins, finite nonnegative scores, positive
component totals, normalized sums, and agreement with the service's total
absolute score. It flags an all-zero target as uninformative. The selected
signature has 2,076 genes, 1,975 overlapping the bulk matrix, rank 5 and condition
number approximately 5.1733 (5.1982 on shared genes). The five inferred target
profiles contain only **two distinct profiles**. Signature q-values, fit
permutation P values, profile count and matrix conditioning are technical
diagnostics; they do not establish biological replication or validity.

The clinical pipeline uses 22 response-evaluable patients after excluding the
two radiotherapy-associated DR* cases from the 24 available bulk samples. It
performs that clinical join separately from reference construction and verifies
the exclusions. This module itself does not use clinical outcomes to define or
tune the reference and does not fit survival models.

Historical high-resolution expression inference is a separate Kelvin workflow.
Do not substitute these absolute component scores for those inferred expression
values, or use an older archived component whose reference identity is not
established. Superseded references and exploratory Cologne KM analyses are
excluded from the default manuscript reproduction.

## Dependencies and local verification

Reference preparation requires R, Seurat, SeuratObject and jsonlite. Python
reference validation requires Python 3.10+, NumPy and pandas. The recorded
service-output validation environment used R 4.6.0, Python 3.13.7, NumPy 2.2.6
and pandas 2.3.3, as recorded in the October 7 archive's selected-job
`audit/dependencies.json`. The repository dependency record includes versions
used for the current reproduction. Reference preparation writes an R session record. Synthetic guard
checks can be run without external data:

```sh
python -m unittest discover \
  -s src/validation/manuscript_2026/reference -p 'test_reference.py'
```
