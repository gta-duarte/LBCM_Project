# LBCM_Project
Copy-number variations (CNVs) are major drivers of genomic instability in cancer. CNVs lead to dosage imbalances that often disrupt gene expression in deleterious ways to cellular balance. This creates context-specific dependencies that might serve as therapeutic targets for cancer treatment and biomarkers in precision oncology. A core challenge of mapping these dependencies is distinguishing driver genes within broad CNV events from co-altered passengers. 

Here, we analyse recent pre-processed whole-exome sequencing CNV data and high-throughput CRISPR-Cas9 knockout datasets, released through the CellModelPassports and DepMap portals. Data from around 570 cancer cell models were used to systematically map CNV dependencies and prioritise candidate drivers within broad co-altered CNV regions. To achieve this, we clustered CNV profiles by identity and carried out an analysis of covariance (ANCOVA) to test CNV cluster -- essentiality associations. We recovered established signals in the literature such as ERBB2, CCND1 and MDM2 dependencies, suggesting events of oncogene addiction, with high statistical support. 

Our driver prioritisation analysis integrated transcriptomics to refine clusters from top ranked associations into candidate driver -- target gene pairs. We prioritised multiple putative co-drivers of described cancer gene dependencies like STARD3, MDM4 and CDK4, validating our methodological approach. This project highlights the importance of multi-omics integration for comprehensive cancer dependency mapping and biomarker discovery.

## Repository Structure
```text
├── R/           # Analysis scripts, functions, and pipelines
├── data/        # Raw and processed datasets
├── output/      # Figures, tables, and computational results
├── LICENSE      # GPLv3 license
└── README.md    # This file
```

## Requirements
- R ≥ 4.5.0
- Package dependencies will be provided via `renv.lock` or `DESCRIPTION` upon code release.

## License
This project is licensed under the GNU General Public License v3.0. See `LICENSE` for details.
