# Mosquito_olfaction_RNAseq_analysis

This repository contains the core R scripts used for the advanced statistical and bioinformatics analyses presented in the manuscript “Human scent primes the mosquito midgut for the metabolic demands of blood feeding”.

Included Scripts:
a_GO and KEGG for DEG scripts.R: Scripts for functional enrichment analyses.
b_GSEA script.R: Scripts for Gene Set Enrichment Analysis (GSEA)
c_WGCNA script.R: Scripts for co-expression network construction, module identification, and module-trait correlation.
d_Physiological assays scripts.R: Scripts for statistical testing (GLMM/LMM) of physiological and behavioral data.

Note on Variable Names:
Please note that variable/group names in the original analysis scripts may differ slightly from the finalized terminology in the manuscript (e.g., Human odor = Odor stimulation, Normal air = Passive exposure).

Data Availability in Code:
Please note that the raw data tables (e.g., RNA-seq count matrices) and intermediate files are not included in this repository due to file size limitations. The scripts provided here are designed to be run by replacing the input file paths with the corresponding datasets, which are available in the NCBI as stated in the manuscript.
