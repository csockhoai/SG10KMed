# SG10K_Med
The code here is used for analysis of the dataset and generation of all figures reported in the manuscript "Analysis of human disease variants from ancestrally diverse Asian genomes" by the SG10K_Health Med working group. The data files required for running the codes involve individual-level data, which will be made available for the purposes of replicating the results to researchers registered through the SG10K_Health Data Access Portal (https://www.npm.sg/collaborate/partners/sg10k/) and the National Precision Medicine (NPM) Project Data Access Committee (contact_npco@gis.a-star.edu.sg).

## General requirements

R 4.1.0 or above is required to run the code. Successful installation of the following R packages is required for execution of the code: tidyverse (1.3.1), RVAideMemoire (v 0.9-80), gplots (3.1.1), patchwork (1.1.1), janitor (2.1.0), data.table (1.14.0). Some older versions might also work.

Each script is standalone and can be performed independently of other scripts in this repository. The data files required to run the code in each script is specified within the script.

## Contents

### List of scripts

**demographics.R**
To tabulate basic demographics and overall pathogenic (P/LP) variant carriers, as displayed in Supplementary Table 2

**vartype.R**
To tabulate the list and number of pathogenic (P/LP) variants by consequence type, as described in Results and Supplementary Data 2

**acmgSF.R**
To analyse carriers of ACMG secondary findings (SF) genes in SG10K_Health, as displayed in Supplementary Table 3, Figure 1A and Supplementary Figure 2

**carrierdist_bygene.R**
To analyse carrier distribution of pathogenic/likely pathogenic (P/LP) variant carriers in dominant and recessive conditions genes in SG10K_Health, as displayed in Figure 1B-C, Supplementary Figure 1

**cnvdel.R**
To generate summary figure of CNV deletions identified by size of deletion and carrier frequency

**admixture.R**
To perform admixture analysis, as displayed in Figures 2A-B, Supplementary Figure 4, Supplementary Tables 4-5

**VUS-FP.R**
To perform analysis for identifying and validation of VUS-favour pathogenic variants (VUS-FP), as displayed in Figures 2C-D

**PGXanalysis.R**
To analyze the fraction of individuals with a risk allele in one of the 23 pharmacogenes and the fractions of individuals with concurrent P/LP variants in a CDC Tier1 genetic condition and a PGX risk allele relevant to their condition, as displayed in Figure 3

**PGXfreq.R**
This codes calculates the carrier frequency of pharmacophenotypes for 23 pharmacogenes with high confidence gene-drug interactions, allele frequencies and diplotype frequencies for genes with star allele nomenclature, as displayed in Table 2, Supplementary Table 6 and Supplementary Data 7,10

**PGXQC.R**
This code is to compare the correlation of 15X vs 30X WGS samples in PGX genes for potential batch effects as displayed in Supplementary Figure 5

### Troubleshooting guide

The code here have been tested in RStudio (1.4.1717). Below are potential troubleshooting guides.

**Error type 1: invalid graphics state**
When encountering 'invalid graphics state', please clear your session plot history and re-run the code.

**Error type 2: figure margins too large**
This error can be addressed by enlarging the plot window.
