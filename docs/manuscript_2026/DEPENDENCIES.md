# Dependencies

Use a separate Python environment and an R installation available as `Rscript`.
The workflow is local: no browser, service credentials or
Illustrator is required for analysis and standalone rendering.

## R versions recorded in the approved replay

The October 7 reproducibility archive records the following versions. These are
tested versions, not claims that other versions produce different answers.

| Component | Version | Role |
|---|---|---|
| R | 4.6.0 | Statistical models and rendering |
| survival | 3.8-6 | Kaplan–Meier, log-rank and Cox models |
| metafor | 5.0-1 | REML pooling |
| jsonlite | 2.0.0 | Validation records |
| svglite | 2.2.2 | Editable SVG export |
| ragg | 1.5.2 | 600-dpi PNG export |
| grid / grDevices / stats | R 4.6.0 base libraries | Layout, standard PDF device and models |
| Seurat | 5.5.1 | Optional annotated-reference rebuild |
| SeuratObject | 5.4.0 | Optional annotated-reference rebuild |

The source plotting session also records systemfonts 1.3.2 and textshaping 1.0.5.
Use Helvetica or a metrically compatible sans-serif font and inspect exports
after font substitution. PDF uses the standard R PDF device; Cairo is not
required. Install only the reference packages if reconstructing the annotated
EcoTyper reference rather than replaying the recorded deconvolution output.

The version evidence is in the source bundle's
`validation/replay/{FINAL_Figure5e,FINAL_Figure5f}/R_sessionInfo.txt`,
`validation/replay/LymphoMAP_KMs/audit/R_sessionInfo.txt` and
`reference/annotations/R_sessionInfo.txt`. The workflow records the environment
used for each new run.

## Python

Input preparation and independent verification require Python with NumPy,
pandas, SciPy, statsmodels and openpyxl. This public module was tested with
Python 3.13.7, NumPy 2.2.6, pandas 2.3.3, SciPy 1.16.2, statsmodels 0.14.5 and
openpyxl 3.1.5. The October 7 source bundle also records the same Python,
NumPy, pandas and SciPy versions in
`deconvolution/selected_AND_jobs34_35/audit/dependencies.json`.
The statsmodels and openpyxl pins are verified from the October 13 environment;
the source analysis dependency list names them without version pins.

```sh
python3.13 -m venv /path/to/validation_environment
. /path/to/validation_environment/bin/activate
python3 -m pip install -r docs/manuscript_2026/requirements.txt
```

Reference construction must accommodate the 15,935-cell Seurat object and the
20,674-feature selected matrix. Recorded numerical replay is substantially
smaller and does not need to load the RDS. Numerical and synthetic checks should
be run after changing an environment; agreement with the expected coefficients,
confidence intervals, pooled results and KM tables is the acceptance criterion.
