# 2026 manuscript validation

This module reproduces the selected manuscript validation analyses from the
recorded October 7, 2026 reproducibility inputs. It adds LymphoMAP CAR-T and
Cologne while preserving the established clinical definitions and original
chemotherapy OS method. It does not rerun discovery analyses or alter manuscript
text.

## Run

Install the [dependencies](DEPENDENCIES.md), obtain the external input bundle as
described in [PROVENANCE.md](PROVENANCE.md), and run from the repository root:

```sh
python3 src/validation/manuscript_2026/stage_inputs.py \
  --bundle /path/to/extracted/Reproducibility_20261007 \
  --output-dir /path/to/external_inputs
python3 src/validation/manuscript_2026/run.py \
  --input-dir /path/to/external_inputs \
  --output-dir /path/to/new_results
```

The first command stages the required inputs from the extracted source bundle.
The second prepares the analysis populations, runs R models, independently checks
the results and renders standalone plots from the newly generated numerical
tables. Use a new output directory. Keep external inputs and generated
patient-level outputs outside the public checkout.

The output directory contains `prepared/` (analysis inputs and exclusions),
`numerical/` (fitted models, pools and KM tables), `verification/` (checks),
`plots/` (standalone exports) and `run_manifest.json` (run provenance).

The default workflow uses recorded CIBERSORTx outputs. It requires no browser,
service credentials, network access or Illustrator. It does not start or change
CIBERSORTx jobs.

To also rebuild the annotated EcoTyper reference, verify exact two-cell
oversampling and check the recorded signature, scores and service logs, include
the external reference files during staging and add `--rebuild-reference`:

```sh
python3 src/validation/manuscript_2026/stage_inputs.py \
  --bundle /path/to/extracted/Reproducibility_20261007 \
  --output-dir /path/to/external_inputs_with_reference \
  --include-reference
python3 src/validation/manuscript_2026/run.py \
  --input-dir /path/to/external_inputs_with_reference \
  --output-dir /path/to/new_results_with_reference \
  --rebuild-reference
```

This runs the same clinical analysis plus reference validation under
`reference/` in the new output directory. It still uses recorded GUI outputs
and submits no service jobs. It requires the optional Seurat dependencies and
additional memory for the RDS/reference matrix. See [CIBERSORTx.md](CIBERSORTx.md)
for each reference step and service settings. Use `--rscript /path/to/Rscript`
if the recorded R environment is not the default on your PATH. Every staged
input and generated output directory must be new and outside this checkout.

## Analyses and figure mapping

Figure numbers below identify the approved figure release; the manuscript may
subsequently be renumbered.

| Output | Analysis | Approved figure use |
|---|---|---|
| LymphoMAP CCL8+ KM | CAR-T EFS, marker-positive myeloid cells / all measured cells | Figure 3d |
| LymphoMAP CCL13+ KM | Same eligibility, CCL13 gate | Figure 4c |
| LymphoMAP CCL8+CCL13+ KM | Both markers positive in the same myeloid cell | Figure 5d |
| LymphoMAP CCL18+ KM | Same eligibility, CCL18 gate | Figure 6c |
| Chemotherapy OS forest | CHOP I, CHOP II and R-CHOP I; continuous per-IQR Cox effects | Figure 5e |
| CAR-T durability forest | Axi-cel/Moffitt, ZUMA-7 CAR-T, Cologne and LymphoMAP; continuous per-IQR logistic effects | Figure 5f |

Standalone plots are editable SVG/PDF and 600-dpi PNG. Native figure assembly,
panel rearrangements and unchanged molecular/clinical artwork are separate
Illustrator assets; this module does not reconstruct complete manuscript figures.

## Definitions that must remain fixed

**LymphoMAP:** one eligible biopsy immediately before CAR-T per patient, giving
22 patients. Published Mac, Mo, DC and pDC annotations identify myeloid cells.
Raw RNA count > 0 defines marker positivity; the AND group requires both markers
in the same cell. Positive myeloid-cell counts are divided by all measured cells,
not only myeloid cells. Single-marker groups can overlap. Valid samples with no
positive cells receive zero fractions. EFS years are converted to months.

KM splits use High > median and Low ≤ median among endpoint-evaluable patients.
The High/Low counts are 11/11 (CCL8), 9/13 (CCL13), 11/11 (CCL18) and 7/15 (AND).
Zero medians and ties are preserved. The analysis uses log-rank tests, high-versus-low
Cox models with Efron ties, and 95% log–log KM confidence limits. Curves, confidence
limits, censor marks, risk counts and model statistics are checked against the
approved tables. These fractions differ from Kelvin's single-marker bulk KMs,
which use CIBERSORTx-inferred myeloid/DC gene expression.

**Cologne:** all 24 bulk samples are deconvolved; two DR* cases that received
radiotherapy within six months are excluded from response models, leaving
22 patients (5 durable, 17 non-durable). The predictor is the AND absolute score
divided by the sum of all five retained component scores. This is a normalized
fitted-component share, not a within-myeloid fraction or a directly measured cell
count. The selected reference contains only two unique AND cells, each copied
ten times. It is exploratory; duplication adds no biological replication, and
clinical model confidence intervals do not include reference uncertainty.

**Forests:** predictor = (fraction − cohort median) / cohort IQR, using R's type-7
quartiles. Chemotherapy scaling follows Kelvin's calculation on the joined
predictor/clinical records before Cox complete-case omission. CAR-T scaling uses
the endpoint-evaluable population. A zero IQR or non-estimable model is flagged,
not replaced with a different scale or a null effect. Log-HRs and log-odds ratios
are pooled separately by REML with normal-theory Wald inference. Forests display
natural logs: zero is the null, negative is better outcome, positive is worse.

The CAR-T logistic outcome is non-durable = 1. Axi-cel/Moffitt uses study-defined
DR; Cologne uses six-month complete remission; ZUMA-7 and LymphoMAP use
EFS-derived six-month durability. Events at or before six months are non-durable;
follow-up beyond six months or event-free censoring exactly at six months is
durable. Earlier event-free censoring and invalid endpoints are excluded. These
related endpoints are not a uniform measured six-month response visit.

## Expected selected results and boundaries

| Forest | Patients | Cohorts | P | I² |
|---|---:|---:|---:|---:|
| Chemotherapy OS | 543 | 3 | 8.44127802277423 × 10⁻⁵ | 8.47614587470678% (displayed 8%) |
| CAR-T durability | 213 | 4 | 0.0443313761026975 | 0% |

The default output excludes LymphoMAP frontline analyses, exploratory Cologne
KMs, superseded forest versions and ZUMA-1 best response. ZUMA-1 contributes
separate manuscript CR/PD boxplots, not the durability pool. Its 57-patient
best-response population is cohort accounting, not an additional current forest.

The repository's existing scripts remain historical source material. In
particular, the September total-mixture forest using an archived, unverified
Stahl component and the pooled response analysis that included both ZUMA-7 arms
must not be substituted for this selected analysis. See
[PROVENANCE.md](PROVENANCE.md) for the source-to-workflow mapping.

See the [validation record](VALIDATION.md) for independent checks, synthetic-test
commands, exact reproduction results and export-review limitations.
