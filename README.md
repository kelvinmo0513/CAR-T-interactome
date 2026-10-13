# CAR-T-interactome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-pending-blue.svg)](https://doi.org/pending)

Single-cell interactome analysis reveals CCL8⁺CCL13⁺ tumor-associated macrophages as drivers of CAR T cell resistance in large B cell lymphoma. This repository contains computational pipelines for identifying resistance mechanisms through on-treatment single-cell profiling and cross-cohort validation.

This repository contains computational pipelines and analysis code for identifying resistance mechanisms through on-treatment single-cell profiling and cross-cohort validation in CAR T cell therapy for large B cell lymphoma (LBCL).
<br>
<br>
<img width="5244" height="5243" alt="Cancer Cell Graphical Abstract" src="https://github.com/user-attachments/assets/b1397334-068d-4ac9-9178-0c24b7fc129c" />
<br>
<br>
We performed the first on-treatment single-cell characterization of the tumor microenvironment following CAR T cell therapy in large B cell lymphoma (LBCL). Through comprehensive single-cell multi-omic profiling and interactome analysis, we identified CCL8⁺CCL13⁺ tumor-associated macrophages (TAMs) as key mediators of CAR T resistance.
<br>
<br>
<img width="3106" height="1125" alt="image" src="https://github.com/user-attachments/assets/8e1e1470-75ba-43e1-ba87-81d72581fb24" />
<br>
<br>

## Key Findings

- **CCL8⁺CCL13⁺ TAMs** are highly proliferative and immunosuppressive, consistently predicting poor outcomes across CAR T therapy and chemoimmunotherapy
- **IRF4/IRF8-driven dysfunction** in CAR T cells operates independently of canonical TOX-mediated exhaustion pathways
- **CCL18⁺ TAMs** promote immunosuppression through interactions with CCR8⁺ regulatory T cells
- Resistance mechanisms are conserved between CAR T cell therapy and conventional chemoimmunotherapy
- Enhanced prognostic accuracy compared to conventional macrophage markers (CD68⁺ and CD163⁺ TAMs)

## Repository Structure
```
CAR-T-interactome/
├── data/                    # Data processing scripts and references
├── src/                     # Core analysis pipelines
│   ├── preprocessing/       # Single-cell data QC and normalization
│   ├── interactome/         # REMI, CellPhoneDB, CellChat analyses
│   ├── deconvolution/       # CIBERSORTx bulk RNA-seq analysis
│   └── validation/          # Cross-cohort validation scripts
├── figures/                 # Figure generation scripts
└── README.md
```

## Current manuscript validation workflow

The [2026 manuscript validation guide](docs/manuscript_2026/README.md) documents the
LymphoMAP CAR-T analyses and the Cologne bulk deconvolution addition. Its portable
workflow reproduces four LymphoMAP EFS Kaplan–Meier plots, the selected four-cohort
durability meta-analysis and the preserved three-cohort chemotherapy OS forest.
It uses the approved total-mixture/cell-fraction definitions and recorded CIBERSORTx
outputs; it does not submit service jobs.

```sh
python3 src/validation/manuscript_2026/stage_inputs.py \
  --bundle /path/to/extracted/Reproducibility_20261007 \
  --output-dir /path/to/external_inputs
python3 src/validation/manuscript_2026/run.py \
  --input-dir /path/to/external_inputs \
  --output-dir /path/to/new_results
```

Install the [documented dependencies](docs/manuscript_2026/DEPENDENCIES.md) first.
The extracted input bundle is available from the study authors; it is not included
in this public repository. See [provenance and input access](docs/manuscript_2026/PROVENANCE.md).
The existing R Markdown scripts and figure examples remain available as historical
workflows; they do not, by themselves, reproduce the selected 2026 additions.
  
## Dataset

### Discovery Cohort
- **Patients**: 18 LBCL + 1 B-ALL receiving CAR T therapy
- **Timing**: Biopsies collected day 2 (axi-cel) or day 7 (CD19/CD22-CAR) post-infusion
- **Cell count**: 31,808 high-quality cells passing QC metrics
- **Technologies**: 
  - scRNA-seq (whole transcriptome)
  - scTCR-seq (TCR clonotyping)
  - CITE-seq (surface epitope profiling)
- **Platform**: 10X Genomics 5' Immune Profiling with Feature Barcoding

### Validation Cohorts
- **Nine treatment-defined cohorts**: CHOP I, CHOP II, R-CHOP I, Axi-cel/Moffitt,
  ZUMA-1, ZUMA-7 CAR-T, LymphoMAP CAR-T, Cologne and ZUMA-7 standard of care
  (labeled HDT-ASCT in the figures).
- **Data types**: Bulk RNA-seq, microarray, NanoString and single-nucleus RNA-seq.
- **Counts and endpoints**: Available profiles, baseline specimens and
  endpoint-evaluable patients are different populations. The
  [cohort table](docs/manuscript_2026/COHORTS.md) gives accession-specific counts,
  exclusions, treatments and the analyses to which each cohort contributes.
  The nine treatment-defined rows are not nine independent accession series.

## Installation
```bash
# Clone the repository
git clone https://github.com/kelvinmo0513/CAR-T-interactome.git
cd CAR-T-interactome
```

## Key Analysis Components

### Interactome Inference Methods

1. **REMI** (Regularized Microenvironment Interactome)
   - Graph-based conditional probability modeling
   - Community detection in ligand-receptor networks
   - Optimal for small sample sizes with high-dimensional data

2. **CellPhoneDB**
   - Curated ligand-receptor interaction database
   - Permutation-based statistical testing
   - Ranked interaction pairs across cell populations

3. **CellChat**
   - Mass action modeling of communication probability
   - Accounts for complex composition of LR pairs
   - Network visualization and differential analysis
  
### Cell Type Annotations

Cell populations were annotated using:
- **Azimuth**: Reference-based projection (Bone Marrow reference)
- **SingleR**: Automated annotation (Monaco Immune reference)
- Manual curation based on canonical marker expression

Major cell types analyzed:
- Conventional T cells (Tconv)
- Regulatory T cells (Treg)
- B cells
- Myeloid/Dendritic cells

## Data Availability

### Raw Single-Cell Data
- **Repository**: Gene Expression Omnibus (GEO) and dbGaP
- **Accession**: Will be available upon manuscript acceptance
- **Contents**: scRNA-seq, scTCR-seq, and CITE-seq data from 19 patients

### Reference Single-Cell Dataset
- **Source**: Lymphoma EcoTyper study
- **URL**: https://ecotyper.stanford.edu/lymphoma
- **GEO Accession**: GSE182436
- **Used for**: CIBERSORTx deconvolution of baseline samples

The public EcoTyper source is distinct from the study's annotated Seurat object
and exact CIBERSORTx upload files. Those inputs remain external to this code
repository. [Reference provenance and access](docs/manuscript_2026/PROVENANCE.md)
and the [CIBERSORTx execution record](docs/manuscript_2026/CIBERSORTx.md) describe
the selected Cologne workflow, including its exploratory two-cell reference.

### Validation Cohorts (Publicly Available)
- **GSE34171**: CHOP I; 176 expression profiles, 129 OS-evaluable patients.
- **GSE10846**: CHOP II (181) and R-CHOP I (233), analyzed separately in the OS forest.
- **GSE153439**: Axi-cel/Moffitt; 36 expression profiles, 35 response-evaluable patients.
- **GSE197977**: ZUMA-1; 74 expression profiles, 58 baseline samples and 57 with known
  best response; the manuscript CR/PD boxplots use 40 patients.
- **GSE248835**: ZUMA-7; 134 CAR-T and 122 standard-of-care screening specimens with EFS.
- **LymphoMAP**: [CELLxGENE collection 4429f25e-5c91-4981-9cf2-6a1044c36732](https://cellxgene.cziscience.com/collections/4429f25e-5c91-4981-9cf2-6a1044c36732);
  22 eligible CAR-T patients; clinical annotations in Li et al. supplementary Table S1B.
- **Cologne**: [Zenodo 10.5281/zenodo.15280550](https://doi.org/10.5281/zenodo.15280550);
  24 bulk samples, 22 response-evaluable patients; clinical annotations in Stahl et al.

See the [full cohort table](docs/manuscript_2026/COHORTS.md) for source studies,
analysis-specific counts and endpoint definitions. Public source access does not
imply that all processed study inputs are deposited in this repository.

### Clinical Trial Data
- **ZUMA-1**: Provided by Kite Pharma (NCT02348216)
- **CD19/CD22-CAR**: Stanford clinical trial (NCT03233854)

## Clinical Significance

Our findings establish CCL8⁺CCL13⁺ TAMs as:
- **Biomarkers** for treatment stratification
- **Therapeutic targets** for combination strategies
- **Mechanistic links** between myeloid immunosuppression and CAR T resistance

Potential therapeutic implications:
- Selective TAM depletion strategies
- Chemokine axis blockade (CCL8/CCL13)
- CAR T engineering to resist TAM-mediated suppression
- Modulation of IRF4/IRF8 transcriptional programs
  
## Citation

If you use this code or data in your research, please cite:
```bibtex
@article{mo2025ccl8ccl13,
  title={CCL8⁺CCL13⁺ tumor-associated macrophages mark early resistance to CAR T cell therapy in large B cell lymphoma},
  author={Mo, Kelvin C and Kramer, Anne Marijn and Yeh, Christine Yiwen and Hamilton, Mark P and Spiegel, Jay Y and Desai, Moksha H and Ehlinger, Zachary J and Yang, Eric and Ozawa, Michael G and Chen, Yiyun and Rietberg, Skyler and Tsui, Kristin CY and Prabhu, Snehit and Lee, Caroline and Frank, Matthew J and Muffly, Lori and Claire, Gursharan K and Bharadwaj, Sushma and Dahiya, Saurabh and Kong, Katherine A and Sotillo, Elena and Sahaf, Bita and Plevritis, Sylvia K and Miklos, David B and Mackall, Crystal L and Good, Zinaida},
  journal={Pending},
  year={2025},
  doi={Pending}
}
```

## Support

For questions or issues:
- **Technical issues**: Open an issue on GitHub
- **Scientific questions**: Contact Dr. Zinaida Good (zinaida@stanford.edu)
- **Collaboration inquiries**: Contact study leadership
  
## Authors

### Co-First Authors
- **Kelvin C. Mo** - Computational analysis and pipeline development
- **Anne Marijn Kramer** - Data interpretation and manuscript preparation

### Senior Authors
- **Zinaida Good, PhD** - Lead contact, study conception and supervision
- **Crystal L. Mackall, MD** - Study oversight and clinical expertise
- **David B. Miklos, MD, PhD** - Clinical trial leadership

### Institutional Affiliations
- Stanford University School of Medicine
- Center for Cancer Cell Therapy, Stanford Cancer Institute
- Parker Institute for Cancer Immunotherapy
- Weill Cancer Hub West (UCSF-Stanford)

See manuscript for complete author list and contributions.
  
## Acknowledgments

This work was supported by:
- Stanford Center for Cancer Systems Biology (NIH/NCI U54-CA209971)
- American Cancer Society and Stanford Cancer Institute (IRG-23-1074369-01-IRG)
- Sponsored research agreement with Kite Pharma (Gilead Sciences)
- Parker Institute for Cancer Immunotherapy
- NIH/NCI Pathway to Independence Award (K99CA293149, R00CA293149)

We thank all patients and their families for their participation, and the clinical teams at Stanford Hospital for their dedication.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
