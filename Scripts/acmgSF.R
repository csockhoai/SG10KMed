# This code is to analyse carriers of ACMG secondary findings (SF) genes in SG10K_Health
# Generated outcomes are displayed in Supplementary Table 3, Figure 1A and Supplementary Figure 2
# This script has the following sections: 
# 1) analysis of ACMG SF version 3 (73 genes) 
# 2) analysis of ACMG SF version 2 (59 genes)
# 3) generation of Figure1A, Supplementary Figure 2
# 4) summary of carriers in newly added genes in ACMG73 (SF v3.0) compared to ACMG59 (SF v2.0)
# 5) quantification of difference in MYBPC3 gene-level carrier burden between Indian and Chinese
# 6) re-calculation of carrier frequency for ACMG AD genes for SG10K_Med cohort age 50 and younger (Supplementary Data 3)
# 7) testing for statistical significance in proportions of carrier frequency across ancestry for highlighted variants

# Data files required:
# Full_Sample_Info_9051_edit_v2.txt
# All_Tier1_Pathogenic_Variants_r5.3_20211117.txt
# All_Tier1_Pathogenic_Variants_r5.3_20211117_varonly.txt
# Z0_ACMG73list.csv
# ExclusionList_01_APOB_LOFv.txt
# Variant_Level_r5.3_20211117_HBOC.txt

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)
library(RVAideMemoire)

# ================================
# Pre-02 : Load data files
# ================================
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
# individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")
# ACMG SF version 3 (v.3) genes list
acmg73 <- read.csv("Z0_ACMG73list.csv", header = TRUE)
# exclusion list of APOB loss-of-function variants: APOB LOF variants are not associated with autosomal dominant familial hypercholesterolemia
apoblof <- read.table("ExclusionList_01_APOB_LOFv.txt", header = TRUE, sep = "\t")
# list of PLP variants
PLP_var <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117_varonly.txt", header = TRUE, sep = "\t")
# subset of variant-level aggregate data for HBOC genes (BRCA1, BRCA2, PALB2)
PLP_var_HBOC <- read.table("Variant_Level_r5.3_20211117_HBOC.txt", header = TRUE, sep = "\t")

# =================================
# Pre-03 : Key variables
# =================================
# consolidated individual demographic and variant data
indv <- left_join(dem, PLP_indv)

# consolidated number of individuals by genetic ancestry
popn_ancestry <- dem %>% group_by(genetic_ethnicity) %>% count()
popn_ancestry <- as.data.frame(popn_ancestry)

# ========================================================
# Section 1.0 : Analysis of ACMG SF version 3 (73 genes)
# ========================================================
# note1: genotype code - 1=heterozygous, 2=homozygous, 3=hemizygous
# note2: inheritance abbreviation - AD=autosomal dominant, AR=autosomal recessive, XL=X-linked

# =====================================================================
# Section 1.1 : Analysis of ACMG SF version 3 (73 genes) - All genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene <- indv %>% filter(gene_symbol %in% acmg73$Gene) %>% group_by(npm_research_id,
                                                                                genetic_sex,
                                                                                genetic_ethnicity,
                                                                                gene_symbol,
                                                                                genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv <- indv_summary.gene %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# =====================================================================
# Section 1.2 : Analysis of ACMG SF version 3 (73 genes) - AD genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene.AD <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD"] &
                                          !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id,
                                                                                      genetic_sex,
                                                                                      genetic_ethnicity,
                                                                                      gene_symbol,
                                                                                      genotype_code) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AD <- indv_summary.gene.AD %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv.AD %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AD %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AD %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and number of variants per individual
indv_summary_count.per.indv.AD %>% group_by(genetic_ethnicity, Count) %>% summarise(Count = n())

# summary of PLP carriers by disease domains (three main domains: Cancer, Cardiovascular (CVD), Lipid disorder)
# disease domain: Cancer
acmg73_AD.cancer <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$PhenotypeCat3 == "Cancer" & acmg73$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity, gene_symbol, genotype_code) %>% summarise(Count = n())
acmg73_AD.cancer_by.indv <- acmg73_AD.cancer %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())
acmg73_AD.cancer_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())

# disease domain: CVD
acmg73_AD.CVD <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$PhenotypeCat3 == "Cardiovascular" & acmg73$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity, gene_symbol, genotype_code) %>% summarise(Count = n())
acmg73_AD.CVD_by.indv <- acmg73_AD.CVD %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())
acmg73_AD.CVD_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())

# disease domain: Lipid disorder
acmg73_AD.lipid <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$PhenotypeCat3 == "CardioMetabolic" & acmg73$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity, gene_symbol, genotype_code) %>% summarise(Count = n())
acmg73_AD.lipid_by.indv <- acmg73_AD.lipid %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())
acmg73_AD.lipid_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())

# =====================================================================
# Section 1.3 : Analysis of ACMG SF version 3 (73 genes) - AR genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene.AR <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AR"] |
                                          hgvs_c %in% apoblof$hgvs_c) %>% group_by(npm_research_id,
                                                                                   genetic_sex,
                                                                                   genetic_ethnicity,
                                                                                   gene_symbol,
                                                                                   genotype_code) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AR <- indv_summary.gene.AR %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv.AR %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AR %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AR %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# =====================================================================
# Section 1.4 : Analysis of ACMG SF version 3 (73 genes) - XL genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene.XL <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "XL"]) %>% group_by(npm_research_id,
                                                                                                               genetic_sex,
                                                                                                               genetic_ethnicity,
                                                                                                               gene_symbol,
                                                                                                               genotype_code) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.XL <- indv_summary.gene.XL %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv.XL %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.XL %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.XL %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# ==========================================================================================================
# Section 1.5 : Analysis of ACMG SF version 3 (73 genes) - Statistical testing of PLP carrier proportions
# ==========================================================================================================
# To test for statistical significance (Fisher's exact test) of proportions of PLP carriers across ancestry (Supplementary Table 3)

#(a) ACMG SF v3.0 : All genes
acmg73_overall <- indv_summary_count.per.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_overall <- left_join(popn_ancestry,acmg73_overall) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_overall <- tabl_overall %>% select(Count, noncarrier)
ctabl_overall.t <- t(ctabl_overall)
fisher.multcomp(ctabl_overall.t, p.method = "BH")

#(b) ACMG SF v3.0 : AD genes
acmg73_AD <- indv_summary_count.per.indv.AD %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD <- left_join(popn_ancestry,acmg73_AD) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD <- tabl_AD %>% select(Count, noncarrier)
ctabl_AD.t <- t(ctabl_AD)
fisher.multcomp(ctabl_AD.t, p.method = "BH")

#(c) ACMG SF v3.0 : AD genes by disease domain
acmg73_AD_cancer <- acmg73_AD.cancer_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD_cancer <- left_join(popn_ancestry,acmg73_AD_cancer) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD_cancer <- tabl_AD_cancer %>% select(Count, noncarrier)
ctabl_AD_cancer.t <- t(ctabl_AD_cancer)
fisher.multcomp(ctabl_AD_cancer.t, p.method = "BH")

acmg73_AD_CVD <- acmg73_AD.CVD_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD_CVD <- left_join(popn_ancestry,acmg73_AD_CVD) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD_CVD <- tabl_AD_CVD %>% select(Count, noncarrier)
ctabl_AD_CVD.t <- t(ctabl_AD_CVD)
fisher.multcomp(ctabl_AD_CVD.t, p.method = "BH")

acmg73_AD_lipid <- acmg73_AD.lipid_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD_lipid <- left_join(popn_ancestry,acmg73_AD_lipid) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD_lipid <- tabl_AD_lipid %>% select(Count, noncarrier)
ctabl_AD_lipid.t <- t(ctabl_AD_lipid)
fisher.multcomp(ctabl_AD_lipid.t, p.method = "BH")

#(d) ACMG SF v3.0 : AR genes
acmg73_AR <- indv_summary_count.per.indv.AR %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AR <- left_join(popn_ancestry,acmg73_AR) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AR <- tabl_AR %>% select(Count, noncarrier)
ctabl_AR.t <- t(ctabl_AR)
fisher.multcomp(ctabl_AR.t, p.method = "BH")

#(e) ACMG SF v3.0 : XL genes
acmg73_XL <- indv_summary_count.per.indv.XL %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_XL <- left_join(popn_ancestry,acmg73_XL) %>% replace(is.na(.),0)  %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_XL <- tabl_XL %>% select(Count, noncarrier)
ctabl_XL.t <- t(ctabl_XL)
fisher.multcomp(ctabl_XL.t, p.method = "BH")

# ========================================================
# Section 2.0 : Analysis of ACMG SF version 2 (59 genes)
# ========================================================
# note1: genotype code - 1=heterozygous, 2=homozygous, 3=hemizygous
# note2: inheritance abbreviation - AD=autosomal dominant, AR=autosomal recessive, XL=X-linked

# =====================================================================
# Section 2.1 : Analysis of ACMG SF version 2 (59 genes) - All genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3"]) %>% group_by(npm_research_id,
                                                                                                          genetic_sex,
                                                                                                          genetic_ethnicity,
                                                                                                          gene_symbol,
                                                                                                          genotype_code) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.v2 <- indv_summary.gene.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# =====================================================================
# Section 2.2 : Analysis of ACMG SF version 2 (59 genes) - AD genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene.AD.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" &
                                                                          acmg73$Inheritance == "AD"] &
                                             !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id,
                                                                                         genetic_sex,
                                                                                         genetic_ethnicity,
                                                                                         gene_symbol,
                                                                                         genotype_code) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AD.v2 <- indv_summary.gene.AD.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and number of variants per individual
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_ethnicity, Count) %>% summarise(Count = n())

# summary of PLP carriers by disease domains (three main domains: Cancer, Cardiovascular (CVD), Lipid disorder)
# disease domain: Cancer
acmg59_AD.cancer <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" & acmg73$PhenotypeCat3 == "Cancer" & acmg73$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity, gene_symbol, genotype_code) %>% summarise(Count = n())
acmg59_AD.cancer_by.indv <- acmg59_AD.cancer %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())
acmg59_AD.cancer_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())

# disease domain: CVD
acmg59_AD.CVD <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" & acmg73$PhenotypeCat3 == "Cardiovascular" & acmg73$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity, gene_symbol, genotype_code) %>% summarise(Count = n())
acmg59_AD.CVD_by.indv <- acmg59_AD.CVD %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())
acmg59_AD.CVD_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())

# disease domain: Lipid disorder
acmg59_AD.lipid <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" & acmg73$PhenotypeCat3 == "CardioMetabolic" & acmg73$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity, gene_symbol, genotype_code) %>% summarise(Count = n())
acmg59_AD.lipid_by.indv <- acmg59_AD.lipid %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())
acmg59_AD.lipid_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())

# =====================================================================
# Section 2.3 : Analysis of ACMG SF version 2 (59 genes) - AR genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene.AR.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" &
                                                                          acmg73$Inheritance == "AR"] |
                                             hgvs_c %in% apoblof$hgvs_c) %>% group_by(npm_research_id,
                                                                                      genetic_sex,
                                                                                      genetic_ethnicity,
                                                                                      gene_symbol,
                                                                                      genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AR.v2 <- indv_summary.gene.AR.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv.AR.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AR.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AR.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# =====================================================================
# Section 2.4 : Analysis of ACMG SF version 2 (59 genes) - XL genes
# =====================================================================
# individual-level summary of PLP variant carriers by gene and number of variants
indv_summary.gene.XL.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" &
                                                                          acmg73$Inheritance == "XL"]) %>% group_by(npm_research_id,
                                                                                                                    genetic_sex,
                                                                                                                    genetic_ethnicity,
                                                                                                                    gene_symbol,
                                                                                                                    genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.XL.v2 <- indv_summary.gene.XL.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 3
# by genetic ancestry
indv_summary_count.per.indv.XL.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.XL.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.XL.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# ==========================================================================================================
# Section 2.5 : Analysis of ACMG SF version 2 (59 genes) - Statistical testing of PLP carrier proportions
# ==========================================================================================================
# To test for statistical significance (Fisher's exact test) of proportions of PLP carriers across ancestry (Supplementary Table 3)

#(a) ACMG SF v2.0 : All genes
acmg59_overall <- indv_summary_count.per.indv.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_overall_v2 <- left_join(popn_ancestry,acmg59_overall) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_overall_v2 <- tabl_overall_v2 %>% select(Count, noncarrier)
ctabl_overall_v2.t <- t(ctabl_overall_v2)
fisher.multcomp(ctabl_overall_v2.t, p.method = "BH")

#(b) ACMG SF v2.0 : AD genes
acmg59_AD <- indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD_v2 <- left_join(popn_ancestry,acmg59_AD) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD_v2 <- tabl_AD_v2 %>% select(Count, noncarrier)
ctabl_AD_v2.t <- t(ctabl_AD_v2)
fisher.multcomp(ctabl_AD_v2.t, p.method = "BH")

#(c) ACMG SF v2.0 : AD genes by disease domain
acmg59_AD_cancer <- acmg59_AD.cancer_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD_cancer_v2 <- left_join(popn_ancestry,acmg59_AD_cancer) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD_cancer_v2 <- tabl_AD_cancer_v2 %>% select(Count, noncarrier)
ctabl_AD_cancer_v2.t <- t(ctabl_AD_cancer_v2)
fisher.multcomp(ctabl_AD_cancer_v2.t, p.method = "BH")

acmg59_AD_CVD <- acmg59_AD.CVD_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD_CVD_v2 <- left_join(popn_ancestry,acmg59_AD_CVD) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD_CVD_v2 <- tabl_AD_CVD_v2 %>% select(Count, noncarrier)
ctabl_AD_CVD_v2.t <- t(ctabl_AD_CVD_v2)
fisher.multcomp(ctabl_AD_CVD_v2.t, p.method = "BH")

acmg59_AD_lipid <- acmg59_AD.lipid_by.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AD_lipid_v2 <- left_join(popn_ancestry,acmg59_AD_lipid) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AD_lipid_v2 <- tabl_AD_lipid_v2 %>% select(Count, noncarrier)
ctabl_AD_lipid_v2.t <- t(ctabl_AD_lipid_v2)
fisher.multcomp(ctabl_AD_lipid_v2.t, p.method = "BH")

#(d) ACMG SF v2.0 : AR genes
acmg59_AR <- indv_summary_count.per.indv.AR.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_AR_v2 <- left_join(popn_ancestry,acmg59_AR) %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_AR_v2 <- tabl_AR_v2 %>% select(Count, noncarrier)
ctabl_AR_v2.t <- t(ctabl_AR_v2)
fisher.multcomp(ctabl_AR_v2.t, p.method = "BH")

#(e) ACMG SF v2.0 : XL genes
acmg59_XL <- indv_summary_count.per.indv.XL.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
tabl_XL_v2 <- left_join(popn_ancestry,acmg59_XL) %>% replace(is.na(.),0)  %>% rowwise() %>% mutate(noncarrier = (n - Count))
ctabl_XL_v2 <- tabl_XL_v2 %>% select(Count, noncarrier)
ctabl_XL_v2.t <- t(ctabl_XL_v2)
fisher.multcomp(ctabl_XL_v2.t, p.method = "BH")

# ====================================================================
# Section 3.0 : Generation of Figure1A and Supplementary Figure 2
# ====================================================================

# =======================================
# Section 3.1 : Generation of Figure1A
# =======================================
# consolidate the carrier counts of PLP variant carriers of AD genes by ancestry group into a dataframe
acmg73_AD_prevalence <- list(popn_ancestry,acmg73_AD, acmg73_AD_cancer, acmg73_AD_CVD, acmg73_AD_lipid) %>% reduce(left_join, by = "genetic_ethnicity")
colnames(acmg73_AD_prevalence) <- c("Ancestry", "Total", "Overall", "Cancer", "CVD", "Lipid")
acmg73_AD_prevalence

# calculate fraction of carriers by ancestry group
acmg73_AD_prevalence.frac <- acmg73_AD_prevalence %>% mutate(Overall_frac = Overall/Total,
                                                             Cancer_frac = Cancer/Total,
                                                             CVD_frac = CVD/Total,
                                                             Lipid_frac = Lipid/Total) %>% select(Ancestry, Overall_frac, Cancer_frac, CVD_frac, Lipid_frac)
acmg73_AD_prevalence.frac

# plot the adjusted carrier frequency of ACMG SF v3 AD genes by disease domain (Figure 1A)
acmg73_AD_prevalence.frac.t <- data.frame(t(acmg73_AD_prevalence.frac[-1]))
Disease <- matrix(c("Overall", "Cancer", "CVD", "Lipid"), nrow = 4, ncol = 1)
prev_tabl <- cbind.data.frame(Disease, acmg73_AD_prevalence.frac.t)
colnames(prev_tabl) <- c("Disease", "CH", "IND", "MY")
prev_tabl.long <- pivot_longer(prev_tabl, cols = 2:4, names_to = "Ancestry", values_to = "AdjCF")
prev_tabl.long$Ancestry <- as.factor(prev_tabl.long$Ancestry)
prev_tabl.long$Disease <- factor(prev_tabl.long$Disease, levels = unique(prev_tabl.long$Disease))
ancestrycol <- c('#d73027', '#4575b4', '#fdae61')

F1A <- ggplot(prev_tabl.long, aes(x = Disease, y = AdjCF, fill = Ancestry)) + geom_bar(stat = "identity", position = position_dodge()) + scale_fill_manual(values = ancestrycol) + theme(panel.background = element_rect(
  fill = "white",
  colour = "white",
  size = 0.5
),
panel.grid.major = element_line(size = 0.5, linetype = "solid"),
panel.border = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_line(
  size = 0.5,
  linetype = "solid",
  colour = "black"
),
axis.line.x = element_blank(),
axis.text.y = element_text(vjust = 0.2, hjust = 0.95, size = rel(0.85))
) + scale_y_continuous(labels = scales::percent) + ylab("Adjusted carrier frequency (%)")
F1A

# =====================================================
# Section 3.2 : Generation of Supplementary Figure 2
# =====================================================
# List of all PLP variants in HBOC genes (BRCA1, BRCA2, PALB2)
PLP_var_HBOC
# replace NA with 0 for carrier frequency dataframe for PLP variants in HBOC genes
PLP_var_HBOC[PLP_var_HBOC == 0] <- NA

# plot adjusted carrier frequency of HBOC genes (BRCA1, BRCA2, PALB2) by variant and ancestry group (Supplementary Figure 2)
PLP_var_HBOC.a <- pivot_longer(PLP_var_HBOC, cols = 10:12, names_to = "Ancestry", values_to = "Carrier_Freq")
PLP_var_HBOC.a$Ancestry <- as.factor(PLP_var_HBOC.a$Ancestry)
PLP_var_HBOC.a$Variant <- factor(PLP_var_HBOC.a$Variant, levels = rev(unique(PLP_var_HBOC.a$Variant)))
ancestrycol <- c('#d73027', '#4575b4', '#fdae61')

F15 <-
  ggplot(PLP_var_HBOC.a, aes(x = Variant, y = Carrier_Freq, na.rm = TRUE)) + geom_point(aes(shape = Ancestry, color = Ancestry),
                                                                                        shape = 16, size = 3,
                                                                                        alpha = 0.65) + scale_color_manual(values = ancestrycol)

F15 + theme(
  panel.background = element_rect(
    fill = "white",
    colour = "white",
    size = 0.5
  ),
  panel.grid.major = element_line(size = 0.5, linetype = "solid"),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_line(
    size = 0.1,
    linetype = 'dotted',
    colour = "#b3cbcb"
  ),
  panel.grid.major.y = element_line(
    size = 0.1,
    linetype = 'dotted',
    colour = "#b3cbcb"
  ),
  axis.line = element_line(
    size = 0.5,
    linetype = "solid",
    colour = "black"
  ),
  axis.text.x = element_text(
    vjust = 1,
    size = rel(0.85),
    angle = 90
  ),
  legend.background = element_blank(),
  legend.key = element_blank(),
  axis.title.y = element_text(vjust = 3)
) + ylab("Normalized carrier frequency") + coord_flip()

# ======================================================================================================================
# Section 4.0 : Carriers in newly added genes in ACMG73 (SF v3.0) compared to ACMG59 (SF v2.0) - mentioned in Results
# ======================================================================================================================
# List of all PLP variants
PLP_var
# PLP variants in ACMG SF v2 list of AD genes
PLP_AD_varcount.acmg59 <- PLP_var %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD" & acmg73$Version != 3] & !(hgvs_c %in% apoblof$hgvs_c))
# number of unique genes with PLP variants
length(unique(PLP_AD_varcount.acmg59$gene_symbol))
#write.table(PLP_AD_varcount.acmg59, "List_PLP_AD_ACMG59_variants.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# PLP variants in AD genes that are newly added to ACMG SF v3
PLP_AD_varcount.acmg73_only <- PLP_var %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD" & acmg73$Version == 3] & !(hgvs_c %in% apoblof$hgvs_c))
# number of unique genes with PLP variants
length(unique(PLP_AD_varcount.acmg73_only$gene_symbol))
#write.table(PLP_AD_varcount.acmg73_only, "List_PLP_AD_ACMG73only_variants.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# number of TTN P/LP carriers
TTN_carriers <- PLP_indv %>% filter(hgvs_c %in% PLP_AD_varcount.acmg73_only$hgvs_c[PLP_AD_varcount.acmg73_only$gene_symbol == "TTN"])
length(TTN_carriers$npm_research_id)

# number of PALB2 P/LP carriers
PALB2_carriers <- PLP_indv %>% filter(hgvs_c %in% PLP_AD_varcount.acmg73_only$hgvs_c[PLP_AD_varcount.acmg73_only$gene_symbol == "PALB2"])
length(PALB2_carriers$npm_research_id)

# =============================================================================================================================
# Section 5.0 : To quantify difference in MYBPC3 gene-level carrier burden between Indian and Chinese - mentioned in Results
# =============================================================================================================================
# List of all MYBPC3 carriers
MYBPC3_var.carriers <- PLP_indv %>% filter(gene_symbol == "MYBPC3") %>% select(hgvs_c, npm_research_id)
MYBPC3_carriers <- left_join(MYBPC3_var.carriers, dem, by = "npm_research_id") %>% group_by(genetic_ethnicity) %>% count() %>% rename(mybpc3_carrier = n) %>% left_join(., popn_ancestry) %>% rowwise() %>% mutate(frac_carrier = mybpc3_carrier/n)

# MYBPC3 carriers by ancestry group
MYBPC3_c <- MYBPC3_carriers$frac_carrier[MYBPC3_carriers$genetic_ethnicity == "C"]
MYBPC3_i <- MYBPC3_carriers$frac_carrier[MYBPC3_carriers$genetic_ethnicity == "I"]

# Quantify fold difference between Indian and Chinese groups
MYBPC3_IoverC <- MYBPC3_i/MYBPC3_c
MYBPC3_IoverC

# ==========================================================================================================
# Section 6.0 : To evaluate for potential survivorship bias influencing carrier frequency in ACMG AD genes 
# ==========================================================================================================
# To re-calculate carrier frequency for ACMG AD genes for SG10K_Med cohort age 50 and younger (Supplementary Data 3)

# ===========================================================================
# Section 6.1 : Consolidate and summarize SG10K_Health cohort age <= 50y 
# ===========================================================================
# subset demographic data for individuals age 50 and younger
dem_less50 <- dem %>% filter(!(is.na(age)) & age <= 50)

# summarize age distribution by ancestry
dem_less50_ancestry_sum <- dem_less50 %>% group_by(genetic_ethnicity) %>% summarise(Age_median = median(age), Age_min = min(age), Age_max = max(age))
dem_less50_ancestry <- dem_less50 %>% select(genetic_ethnicity, age) %>% group_by(genetic_ethnicity)
dem_less50_ancestry_total <- dem_less50_ancestry %>% group_by(genetic_ethnicity) %>% count(name = "Count")

# displaying distribution of cohort age by genetic ancestry and testing for difference in age distribution
# statistical test using ANOVA followed by Tukey test
table(dem_less50_ancestry$genetic_ethnicity)
myaov <- aov(dem_less50_ancestry$age ~ dem_less50_ancestry$genetic_ethnicity)
summary(myaov)
TukeyHSD(myaov)
plot(density(dem_less50_ancestry$age[which(dem_less50_ancestry$genetic_ethnicity == "C")]), col = "red")
lines(density(dem_less50_ancestry$age[which(dem_less50_ancestry$genetic_ethnicity == "I")]), col = "green", add = T)
lines(density(dem_less50_ancestry$age[which(dem_less50_ancestry$genetic_ethnicity == "M")]), col = "blue", add = T)

# subset PLP individual-level dataset for individuals age 50 and younger
PLP_indv_less50 <- PLP_indv %>% filter(npm_research_id %in% dem_less50$npm_research_id)
PLP_indv_less50_1a <- PLP_indv_less50 %>% select(gene_symbol, hgvs_g, hgvs_c, hgvs_p, npm_research_id, genotype_code)
dem_less50_sub <- dem_less50 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, age, study, target_depth)
PLP_indv_less50_1b <- left_join(PLP_indv_less50_1a, dem_less50_sub, by = "npm_research_id")

# ==========================================================
# Section 6.2 : Summary of gene-level carrier frequency  
# ==========================================================
# calculate carrier frequency by gene by ancestry
PLP_indv_less50_gene <-
  PLP_indv_less50_1b %>% select(gene_symbol, genetic_ethnicity) %>% group_by(gene_symbol, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = genetic_ethnicity, values_from = n) %>% replace(is.na(.), 0) %>% mutate(
    Total = (C + I + M),
    C_Frac = (C/dem_less50_ancestry_total$Count[dem_less50_ancestry_total$genetic_ethnicity == "C"]),
    I_Frac = (I/dem_less50_ancestry_total$Count[dem_less50_ancestry_total$genetic_ethnicity == "I"]),
    M_Frac = (M/dem_less50_ancestry_total$Count[dem_less50_ancestry_total$genetic_ethnicity == "M"])
  ) %>% rowwise() %>% mutate(max_Frac = pmax(C_Frac, I_Frac, M_Frac))

# filter to include only genes that have carrier frequency >= 0.1% in at least 1 ancestry group
PLP_indv_less50_gene_0.1percentandmore <- PLP_indv_less50_gene %>% filter(max_Frac >= 0.001)

# multiple pairwise comparisons for gene-level carrier frequencies
pval_list_less50 <- c()
fisher_multcomp_results_less50 <- c()

dem_less50_pop_sum <- dem_less50_ancestry_total

for(i in 1:nrow(PLP_indv_less50_gene_0.1percentandmore)){
  myrow <- PLP_indv_less50_gene_0.1percentandmore[i,]
  c_carrier <- as.integer(myrow$C)
  c_noncarrier <- as.integer(dem_less50_pop_sum$Count[which(dem_less50_pop_sum$genetic_ethnicity == "C")]) - as.integer(c_carrier)
  i_carrier <- as.integer(myrow$I)
  i_noncarrier <- as.integer(dem_less50_pop_sum$Count[which(dem_less50_pop_sum$genetic_ethnicity == "I")]) - as.integer(i_carrier)
  m_carrier <- as.integer(myrow$M)
  m_noncarrier <- as.integer(dem_less50_pop_sum$Count[which(dem_less50_pop_sum$genetic_ethnicity == "M")]) - as.integer(m_carrier)
  
  #Create matrix
  carrier_matrix <- matrix(c(c_carrier, c_noncarrier, i_carrier, i_noncarrier, m_carrier, m_noncarrier), 2,3, dimnames = list(carrier = c("Carrier", "Noncarrier"), ethnicity = c("Chinese", "Indian", "Malay")))
  pval <- fisher.test(carrier_matrix, workspace = 2e8)$p.value
  pval_list_less50 <- c(pval_list_less50, pval)
  
  fisherout <- fisher.multcomp(carrier_matrix, p.method = "BH")
  fisher_multcomp_results_less50 <- rbind(fisher_multcomp_results_less50, fisherout$p.value)
}
pval_list_adjusted_less50 <- p.adjust(pval_list_less50, method = "BH")
pval_table_less50 <- cbind.data.frame(PLP_indv_less50_gene_0.1percentandmore, pval_list_less50, pval_list_adjusted_less50, fisher_multcomp_results_less50)

# subset for ACMG SF v3.0 AD genes
pval_table_less50.acmgAD <- pval_table_less50 %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD"])

# ==========================================================
# Section 6.3 : Summary of variant-level carrier frequency  
# ==========================================================
#calculate carrier frequency by variant by ancestry
PLP_indv_less50_var <-
  PLP_indv_less50_1b %>% select(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity) %>% group_by(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = genetic_ethnicity, values_from = n) %>% replace(is.na(.), 0) %>% mutate(
    Total = (C + I + M),
    C_Frac = (C/dem_less50_ancestry_total$Count[dem_less50_ancestry_total$genetic_ethnicity == "C"]),
    I_Frac = (I/dem_less50_ancestry_total$Count[dem_less50_ancestry_total$genetic_ethnicity == "I"]),
    M_Frac = (M/dem_less50_ancestry_total$Count[dem_less50_ancestry_total$genetic_ethnicity == "M"])
  ) %>% rowwise() %>% mutate(max_Frac = pmax(C_Frac, I_Frac, M_Frac))

# filter to include only variants that have carrier frequency >= 0.1% in at least 1 ancestry group
PLP_indv_less50_var_0.1percentandmore <- PLP_indv_less50_var %>% filter(max_Frac >= 0.001)

# multiple pairwise comparisons for variant-level carrier frequencies
pval_list_less50_var <- c()
fisher_multcomp_results_less50_var <- c()

dem_less50_pop_sum <- dem_less50_ancestry_total

for(i in 1:nrow(PLP_indv_less50_var_0.1percentandmore)){
  myrow <- PLP_indv_less50_var_0.1percentandmore[i,]
  c_carrier <- as.integer(myrow$C)
  c_noncarrier <- as.integer(dem_less50_pop_sum$Count[which(dem_less50_pop_sum$genetic_ethnicity == "C")]) - as.integer(c_carrier)
  i_carrier <- as.integer(myrow$I)
  i_noncarrier <- as.integer(dem_less50_pop_sum$Count[which(dem_less50_pop_sum$genetic_ethnicity == "I")]) - as.integer(i_carrier)
  m_carrier <- as.integer(myrow$M)
  m_noncarrier <- as.integer(dem_less50_pop_sum$Count[which(dem_less50_pop_sum$genetic_ethnicity == "M")]) - as.integer(m_carrier)
  
  #Create matrix
  carrier_matrix <- matrix(c(c_carrier, c_noncarrier, i_carrier, i_noncarrier, m_carrier, m_noncarrier), 2,3, dimnames = list(carrier = c("Carrier", "Noncarrier"), ethnicity = c("Chinese", "Indian", "Malay")))
  pval <- fisher.test(carrier_matrix, workspace = 2e8)$p.value
  pval_list_less50_var <- c(pval_list_less50_var, pval)
  
  fisherout <- fisher.multcomp(carrier_matrix, p.method = "BH")
  fisher_multcomp_results_less50_var <- rbind(fisher_multcomp_results_less50_var, fisherout$p.value)
}
pval_list_adjusted_less50_var <- p.adjust(pval_list_less50_var, method = "BH")
pval_table_less50_var <- cbind.data.frame(PLP_indv_less50_var_0.1percentandmore, pval_list_less50_var, pval_list_adjusted_less50_var, fisher_multcomp_results_less50_var)

# subset for ACMG SF v3.0 AD genes
pval_table_less50.acmgAD_var <- pval_table_less50_var %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD"])
write.table(pval_table_less50.acmgAD_var, "Variant_Level_r5.3_under50_withPvalue.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# =============================================================================================================================
# Section 6.4 : [Extra] Summary of PLP variant-level carrier frequency by ancestry for the entire cohort (9051 individuals)
# =============================================================================================================================
# subset demographics data for entire cohort
dem_subset <- dem %>% select(npm_research_id, genetic_ethnicity, age, study, target_depth)
dem_subset_sum <- dem_subset %>% group_by(genetic_ethnicity) %>% count(name = "Count")

# consolidate individual-level PLP data with demographics
PLP_indv_all_var <- PLP_indv %>% select(gene_symbol, hgvs_g, hgvs_c, hgvs_p, npm_research_id, genotype_code)
PLP_indv_all_var.a <- left_join(PLP_indv_all_var, dem_subset, by = "npm_research_id")

# calculate variant-level carrier frequency by ancestry
PLP_indv_all_var.b <-
  PLP_indv_all_var.a %>% select(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity) %>% group_by(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = genetic_ethnicity, values_from = n) %>% replace(is.na(.), 0) %>% mutate(
    Total = (C + I + M),
    C_Frac = (C/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "C"]),
    I_Frac = (I/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "I"]),
    M_Frac = (M/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "M"])
  ) %>% rowwise() %>% mutate(max_Frac = pmax(C_Frac, I_Frac, M_Frac))

# include only variants that have carrier frequency >= 0.1% in at least 1 ancestry group
PLP_indv_all_var_0.1percentandmore <- PLP_indv_all_var.b %>% filter(max_Frac >= 0.001)

# multiple pairwise comparisons for VARIANT-level carrier frequencies
pval_list_all_var <- c()
fisher_multcomp_results_all_var <- c()

dem_all_pop_sum <- dem_subset_sum

for(i in 1:nrow(PLP_indv_all_var_0.1percentandmore)){
  myrow <- PLP_indv_all_var_0.1percentandmore[i,]
  c_carrier <- as.integer(myrow$C)
  c_noncarrier <- as.integer(dem_all_pop_sum$Count[which(dem_all_pop_sum$genetic_ethnicity == "C")]) - as.integer(c_carrier)
  i_carrier <- as.integer(myrow$I)
  i_noncarrier <- as.integer(dem_all_pop_sum$Count[which(dem_all_pop_sum$genetic_ethnicity == "I")]) - as.integer(i_carrier)
  m_carrier <- as.integer(myrow$M)
  m_noncarrier <- as.integer(dem_all_pop_sum$Count[which(dem_all_pop_sum$genetic_ethnicity == "M")]) - as.integer(m_carrier)
  
  #Create matrix
  carrier_matrix <- matrix(c(c_carrier, c_noncarrier, i_carrier, i_noncarrier, m_carrier, m_noncarrier), 2,3, dimnames = list(carrier = c("Carrier", "Noncarrier"), ethnicity = c("Chinese", "Indian", "Malay")))
  pval <- fisher.test(carrier_matrix, workspace = 2e8)$p.value
  pval_list_all_var <- c(pval_list_all_var, pval)
  
  fisherout <- fisher.multcomp(carrier_matrix, p.method = "BH")
  fisher_multcomp_results_all_var <- rbind(fisher_multcomp_results_all_var, fisherout$p.value)
}
pval_list_adjusted_all_var <- p.adjust(pval_list_all_var, method = "BH")
pval_table_all_var <- cbind.data.frame(PLP_indv_all_var_0.1percentandmore, pval_list_all_var, pval_list_adjusted_all_var, fisher_multcomp_results_all_var)

# subset for ACMG SF v3.0 AD genes
pval_table_all.acmgAD_var <- pval_table_all_var %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD"])

# ========================================================================================================================================================
# Section 7.0 : Testing for statistical significance in proportions of carrier frequency across ancestry for highlighted variants - mentioned in Results
# ========================================================================================================================================================
# statistical significance of carriers for MYBPC3 Arg597Gln
MYBPC3_R597Q_carriers <- pval_table_all.acmgAD_var %>% filter(hgvs_c == "NM_000256.3:c.1790G>A")
MYBPC3_R597Q_carriers
MYBPC3_R597Q_carriers$`Chinese:Indian`

# statistical significance of carriers for BRCA1 Asn909Lysfs*6
BRCA1_N909fs_carriers <- pval_table_all.acmgAD_var %>% filter(hgvs_c == "NM_007294.4:c.2726dup")
BRCA1_N909fs_carriers
BRCA1_N909fs_carriers$`Chinese:Malay`
