# Provenance and external inputs

## Source of the selected workflow

The 2026 validation module is adapted from the study's October 7, 2026
reproducibility addendum, `Reproducibility_20261007`. That addendum preserves
the approved analysis scripts, clinical/deconvolution inputs, LymphoMAP cell
annotations and gates, reference object, CIBERSORTx uploads/downloads and
validation records. It is an internal author archive, not a public data
deposition.

The source archive SHA-256 is
`c41cc6cabc7dd388f66a3624a2b32258c8e3eb2741542505c9152bf1fd3b01dc`.
The [input manifest](../../data/manuscript_2026/input_manifest.json) identifies
the required relative paths, roles and checksums. Checksums verify the identity
of an obtained input; they do not provide access to it.

| Manuscript result or preparation | October 7 source within the extracted addendum | Interpretation |
|---|---|---|
| LymphoMAP cell counting and clinical selection | `analysis/base_total_mixture/pipeline/prepare.py` and `source_data/` | Direct marker-positive myeloid-cell fractions among all measured cells |
| Four LymphoMAP CAR-T EFS KMs | `analysis/base_total_mixture/pipeline/analyze.R` | CAR-T branch only; earlier frontline calculations are outside this module |
| Preserved chemotherapy OS forest | `analysis/Figure5e/replay_original_OS.R` | IQR calculated before Cox complete-case omission; do not replace with later base-stage OS calculations |
| Selected CAR-T durability forest | `analysis/selected_Figure5f/analyze.R` | `oversampled_two_cell_AND_exploratory` result; historical comparison rows are not selected |
| EcoTyper preparation | `reference/annotations/` and `reference/frozen_AND_reference/` | Unchanged annotated source object; raw RNA marker gates and SCT reference expression |
| Selected Cologne deconvolution | `deconvolution/selected_AND_jobs34_35/` | Web GUI signature job 34 and fraction job 35, recorded downloads |
| Standalone KM styling | `rendering/compact_LymphoMAP/pipeline/render_KMs_only.R` | Numerical curves and statistics remain fixed; final native layouts are separate |

The October 7 `selected_Figure5f` directory was originally created as a
sensitivity analysis and then selected for the manuscript. Only its selected
four-cohort result is current. Earlier results using an archived Stahl column
are not interchangeable: Kelvin's exact historical signature, original job
settings and referenced original v1 file were unavailable. No claim is made
that this archived column was a verified AND component.

## Obtaining inputs

Public expression/data sources and study-specific clinical supplement references
are listed in [COHORTS.md](COHORTS.md). To reproduce the exact selected results,
request the annotated reference and processed October 7 input bundle from the
study authors through the contact listed in the repository README. Access is
subject to the source studies' terms and the authors' release arrangements.
The public accession for the EcoTyper source does not guarantee byte identity
with the study's annotated Seurat object.

The public repository includes code, documentation, aggregate expected results
and synthetic tests. Patient-level clinical records, cell annotations/gates,
the annotated RDS and large service matrices/downloads remain external. The
staging command copies the specified inputs into an external directory and
validates their identities. It does not recover missing clinical data from
figures or invent annotations.

## Reproduction boundaries

- Local reproduction starts from recorded deconvolution estimates and the
  supplied LymphoMAP annotations/gates. Rebuilding primary sequencing processing,
  discovery annotation or Kelvin's historical high-resolution CIBERSORTx
  expression workflow is outside this addition.
- Reference reconstruction reproduces the exact two-cell AND upload from the
  annotated object. Both cells share one source-sample label; CCL13 is absent
  from the SCT expression features. Twenty target columns represent only two
  unique cells. See [CIBERSORTx.md](CIBERSORTx.md).
- Recorded jobs used the CIBERSORTx **web GUI**, not an API or command-line
  container. The runner imports their outputs; it does not submit new jobs.
  A new web run can differ because the service did not expose its random seed.
- Model confidence intervals condition on the estimated predictor and do not
  propagate reference or deconvolution uncertainty. Numerical reproduction is
  not independent biological validation of this exploratory component.
- Native Illustrator assembly is separate from numerical analysis and
  standalone plot rendering. SVG identifiers, PDF metadata and platform font
  metrics can vary; validation compares scientific values and plotted geometry
  rather than requiring byte-identical newly rendered graphics.

No credentials, signed download URLs, personal service identifiers or internal
cloud links are needed by the public workflow.
