# Validation cohort roster

The manuscript uses nine treatment-defined validation cohorts. CHOP II and
R-CHOP I share GSE10846; the ZUMA-7 CAR-T and standard-of-care arms share
GSE248835. These are not nine independent accession series. Available profiles,
baseline specimens and endpoint-evaluable patients are different denominators;
no aggregate baseline total is asserted here.

| Cohort | Source / accession | Platform | Treatment | Available and analysis populations | Endpoint and current use |
|---|---|---|---|---|---|
| CHOP I | Monti et al., Cancer Cell 2012; [GSE34171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34171) | Microarray | CHOP | 176 expression profiles; 129 OS-evaluable | OS; chemotherapy forest and historical KMs |
| CHOP II | Lenz et al., NEJM 2008; [GSE10846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846) | Microarray | CHOP | 181 with expression and OS | OS; separate forest cohort; combined with R-CHOP I in specified KMs |
| R-CHOP I | Lenz et al., NEJM 2008; [GSE10846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846) | Microarray | R-CHOP | 233 with expression and OS, disjoint from CHOP II | OS; separate forest cohort; 414 combined CHOP II/R-CHOP I patients in specified KMs |
| Axi-cel / Moffitt | Jain et al., Blood 2021; [GSE153439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153439) | NanoString | Axicabtagene ciloleucel | 36 expression profiles; 35 response-evaluable (17 DR, 18 NDR); one insufficient-follow-up exclusion | Study-defined DR/NDR; CAR-T durability forest and boxplots |
| ZUMA-1 | Scholler et al., Nature Medicine 2022; [GSE197977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197977) | NanoString | Axicabtagene ciloleucel | Source table: 74 expression profiles, 73 clinically annotated; baseline subset 58; known best response 57; current CR/PD boxplots 40 (34 CR, 6 PD) | Best response; excluded from durability meta-analysis |
| ZUMA-7 CAR-T | Locke et al., Nature Medicine 2024; [GSE248835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248835) | NanoString | Axicabtagene ciloleucel | 134 screening specimens with EFS; 134 evaluable for six-month EFS-derived durability | EFS KMs and CAR-T durability forest |
| LymphoMAP CAR-T | Li et al., Cancer Cell 2025; [CELLxGENE collection](https://cellxgene.cziscience.com/collections/4429f25e-5c91-4981-9cf2-6a1044c36732), clinical supplementary Table S1B | Single-nucleus RNA-seq | 17 axi-cel; 5 lisocabtagene maraleucel | 22 patients, one eligible biopsy immediately before CAR-T per patient; all EFS/durability-evaluable | Four EFS KMs; six-month EFS-derived durability forest |
| Cologne | Stahl et al., Cancer Cell 2025; [Zenodo 10.5281/zenodo.15280550](https://doi.org/10.5281/zenodo.15280550) | Bulk RNA-seq | Available: 16 tisa-cel, 5 zamto-cel, 3 MB-CART20.1 | 24 bulk samples; 22 response-evaluable after excluding two DR* radiotherapy-associated cases; 5 DR, 17 NDR | Six-month complete remission; CAR-T durability forest using selected exploratory AND component |
| HDT-ASCT | ZUMA-7 standard of care; Locke et al.; [GSE248835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248835) | NanoString | Salvage chemotherapy followed, in responders, by high-dose therapy and autologous transplantation | 122 screening specimens with EFS | Specified EFS KMs; excluded from CAR-T durability forest |

## Endpoint and treatment notes

- Axi-cel, liso-cel and tisa-cel target CD19; zamto-cel targets CD19/CD20;
  MB-CART20.1 targets CD20. The HDT-ASCT row represents the randomized
  standard-of-care arm, not a claim that every patient underwent transplantation.
- Moffitt required at least six months of follow-up. Its published study DR/NDR
  labels are used; this criterion is not an individual EFS landmark.
- Cologne DR* denotes durable responders who received radiotherapy within six
  months after CAR-T. Those two cases are excluded; evaluable products are
  14 tisa-cel, 5 zamto-cel and 3 MB-CART20.1. The deconvolution uses bulk data,
  not Stahl's separate single-cell dataset.
- LymphoMAP includes one axi-cel plus radiotherapy and one liso-cel plus
  ibrutinib case. They meet its biopsy-to-next-treatment eligibility rule.
  Cologne's DR* exclusion is not imposed on a different cohort.
- ZUMA-1 overall response (OR) is CR/PR; non-response (NR) is SD/PD. One unknown
  response is excluded from the 58 baseline samples, yielding 57. Manuscript
  boxplots further restrict to CR versus PD; PR and SD are not in those plots.
- The supplementary CCL8 ZUMA-7 KM includes both treatment arms (134 + 122).
  It is not the CAR-T-only cohort used in the durability forest.
- CHOP denotes cyclophosphamide, doxorubicin, vincristine and prednisone;
  R-CHOP adds rituximab. OS, overall survival; EFS, event-free survival;
  DR/NDR, durable/non-durable response; CR/PR/SD/PD, complete response,
  partial response, stable disease and progressive disease.

LymphoMAP patient-atlas access is through the CELLxGENE collection above. GEO
accessions associated with Li et al. functional experiments do not identify the
patient atlas. The EcoTyper reference source is separate from these cohorts:
[GSE182436](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182436) and
[Lymphoma EcoTyper](https://ecotyper.stanford.edu/lymphoma).
