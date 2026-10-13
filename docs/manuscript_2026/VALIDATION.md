# Validation record — 2026-10-13

The complete entry command, including `--rebuild-reference`, passed from an
isolated copy of this module with only the 31 declared external inputs. No
repository data, authoring tools, network requests or service submissions were
needed. Patient-level outputs and detailed logs remain in the local author
archive; this record contains aggregate results only.

| Check | Result |
|---|---|
| LymphoMAP eligibility | 22 unique immediately-pre-CAR-T patients |
| High/Low splits | CCL8 11/11; CCL13 9/13; CCL18 11/11; AND 7/15 |
| KMs | All four statistics, 85 curve rows and 32 risk-table rows match approved tables |
| Independent survival verification | Direct product-limit/Greenwood calculations, event/censor counts, log-rank and Python Cox checks pass |
| Cologne | 24 unique samples; two DR* exclusions; 22 response-evaluable |
| Component normalization | All five components; independent R/Python agreement within 1 × 10⁻¹² |
| Four-cohort durability pool | n = 213; P = 0.0443313761026975; I² = 0% |
| Chemotherapy OS pool | n = 543; P = 8.44127802277423 × 10⁻⁵; I² = 8.47614587470678% |
| Coefficients, SEs and intervals | Seven cohort fits match approved values and independent Python models |
| Predictor scaling and direction | Unscaled coefficients × IQR agree; outcome/predictor reversal confirms sign |
| Reference reconstruction | R pre-reference and Python oversampled upload are byte-identical to archived files |
| Reference support | Two unique target cells, ten copies each, one source-sample label; 3,757 unchanged background cells |
| Service record | GUI jobs 34–35, retained target, 24 sample joins and recorded settings verified; error logs empty |
| Synthetic guards | Eleven clinical/analysis tests and five reference tests pass |
| Public file scan | No patient-level files, specimen identifiers, credentials, internal cloud links or machine-specific paths in proposed additions |

The independent REML optimizer confirms the CAR-T boundary estimate τ² = 0.
For chemotherapy it reaches a slightly different stopping point from the
preserved R defaults: independent P = 8.442251209245725 × 10⁻⁵ and
I² = 8.4775464%. Field-specific tolerances account for this optimization
precision; the selected R output is checked directly against the historical
result with tighter tolerances. The historical procedure is not changed to
match a different optimizer.

All four 600-dpi KM PNGs and their four-panel comparison sheet are byte-identical
to approved compact exports in the tested environment. Both forests were
inspected against their regenerated tables, including weighted squares,
confidence limits, pooled diamonds, zero null and better-left direction.
Seven PDFs contain editable text of at least six points within page bounds;
SVG text and 600-dpi PNG dimensions were checked. Poppler bounds checks and
visual review were used; an automatic text/stroke collision audit was not
available. Platform fonts and vector metadata can differ on another computer.

These checks establish numerical and rendering reproduction, not independent
biological validation of the exploratory Cologne reference. Duplication does
not add biological replication, and the clinical confidence intervals do not
propagate reference uncertainty.

## Synthetic checks without study inputs

After installing the documented dependencies, run from the repository root:

```sh
python3 -m unittest discover -s tests/manuscript_2026 -v
python3 -m unittest discover \
  -s src/validation/manuscript_2026/reference -p 'test_reference.py' -v
```

Fixtures are synthetic and are created in temporary external directories. The
R integration guard requires R and the analysis dependencies. It confirms that
zero-IQR cohorts are reported as non-estimable and no partial pool is emitted.
