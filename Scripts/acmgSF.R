# This code is to analyse carriers of ACMG secondary findings (SF) genes in SG10K_Health
# Generated outcomes are displayed in Supplementary Table 5, Figure 1A and Extended Data Figure 2
# This script has three sections: 
# 1) analysis of ACMG SF version 3 (73 genes), 
# 2) analysis of ACMG SF version 2 (59 genes),
# 3) generation of Figure1A, Extended Data Figure 2

# Data files required:
# Full_Sample_Info_9051_edit_v2.txt
# All_Tier1_Pathogenic_Variants_r5.3_20211117.txt
# Z0_ACMG73list.csv
# ExclusionList_01_APOB_LOFv.txt
# Variant_Level_r5.3_20211117_HBOC.txt

library(tidyverse)
library(RVAideMemoire)

# load data files:
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
# individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")
# ACMG SF version 3 (v.3) genes list
acmg73 <- read.csv("Z0_ACMG73list.csv", header = TRUE)
# exclusion list of APOB loss-of-function variants: APOB LOF variants are not associated with autosomal dominant familial hypercholesterolemia
apoblof <- read.table("ExclusionList_01_APOB_LOFv.txt", header = TRUE, sep = "\t")

# consolidate individual demographic and variant data
indv <- left_join(dem, PLP_indv)
# consolidate number of individuals by genetic ancestry
popn_ancestry <- dem %>% group_by(genetic_ethnicity) %>% count()
popn_ancestry <- as.data.frame(popn_ancestry)

# Section 1: ACMG SF version 3 (73 genes)
# note1: genotype code - 1=heterozygous, 2=homozygous, 3=hemizygous
# note2: inheritance abbreviation - AD=autosomal dominant, AR=autosomal recessive, XL=X-linked

# individual-level summary of PLP variant carriers by gene and number of variants: All genes in ACMG SF v3
indv_summary.gene <- indv %>% filter(gene_symbol %in% acmg73$Gene) %>% group_by(npm_research_id,
                                                                                genetic_sex,
                                                                                genetic_ethnicity,
                                                                                gene_symbol,
                                                                                genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv <- indv_summary.gene %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 5: Feature = All genes
# by genetic ancestry
indv_summary_count.per.indv %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())


# individual-level summary of PLP variant carriers by gene and number of variants: AD genes in ACMG SF v3
indv_summary.gene.AD <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD"] &
                                          !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id,
                                                                                      genetic_sex,
                                                                                      genetic_ethnicity,
                                                                                      gene_symbol,
                                                                                      genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AD <- indv_summary.gene.AD %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 5: Feature = AD genes
# by genetic ancestry
indv_summary_count.per.indv.AD %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AD %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AD %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and number of variants per individual
indv_summary_count.per.indv.AD %>% group_by(genetic_ethnicity, Count) %>% summarise(Count = n())
# by disease domains (three main domains: Cancer, Cardiovascular (CVD), Lipid disorder)
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


# individual-level summary of PLP variant carriers by gene and number of variants: AR genes in ACMG SF v3
indv_summary.gene.AR <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AR"] |
                                          hgvs_c %in% apoblof$hgvs_c) %>% group_by(npm_research_id,
                                                                                   genetic_sex,
                                                                                   genetic_ethnicity,
                                                                                   gene_symbol,
                                                                                   genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AR <- indv_summary.gene.AR %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())


# number of carriers as tabulated in Supplementary Table 5: Feature = AR genes
# by genetic ancestry
indv_summary_count.per.indv.AR %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AR %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AR %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())


# individual-level summary of PLP variant carriers by gene and number of variants: XL genes in ACMG SF v3
indv_summary.gene.XL <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "XL"]) %>% group_by(npm_research_id,
                                                                                                               genetic_sex,
                                                                                                               genetic_ethnicity,
                                                                                                               gene_symbol,
                                                                                                               genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.XL <- indv_summary.gene.XL %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 5: Feature = XL genes
# by genetic ancestry
indv_summary_count.per.indv.XL %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.XL %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.XL %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())


# To test for statistical significance (Fisher's exact test) of proportions of PLP carriers across ancestry (Supplementary Table 5)
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


# Section 2: ACMG SF version 2 (59 genes)

# individual-level summary of PLP variant carriers by gene and number of variants: All genes in ACMG SF v2
indv_summary.gene.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3"]) %>% group_by(npm_research_id,
                                                                                                          genetic_sex,
                                                                                                          genetic_ethnicity,
                                                                                                          gene_symbol,
                                                                                                          genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.v2 <- indv_summary.gene.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 5: Feature = All genes
# by genetic ancestry
indv_summary_count.per.indv.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by gene and number of variants: AD genes in ACMG SF v2
indv_summary.gene.AD.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" &
                                                                          acmg73$Inheritance == "AD"] &
                                             !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id,
                                                                                         genetic_sex,
                                                                                         genetic_ethnicity,
                                                                                         gene_symbol,
                                                                                         genotype_code) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AD.v2 <- indv_summary.gene.AD.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 5: Feature = AD genes
# by genetic ancestry
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and number of variants per individual
indv_summary_count.per.indv.AD.v2 %>% group_by(genetic_ethnicity, Count) %>% summarise(Count = n())
# by disease domains (three main domains: Cancer, Cardiovascular (CVD), Lipid disorder)
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

# individual-level summary of PLP variant carriers by gene and number of variants: AR genes in ACMG SF v2
indv_summary.gene.AR.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" &
                                                                          acmg73$Inheritance == "AR"] |
                                             hgvs_c %in% apoblof$hgvs_c) %>% group_by(npm_research_id,
                                                                                      genetic_sex,
                                                                                      genetic_ethnicity,
                                                                                      gene_symbol,
                                                                                      genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.AR.v2 <- indv_summary.gene.AR.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 5: Feature = AR genes
# by genetic ancestry
indv_summary_count.per.indv.AR.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.AR.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.AR.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())

# individual-level summary of PLP variant carriers by gene and number of variants: XL genes in ACMG SF v2
indv_summary.gene.XL.v2 <- indv %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Version != "3" &
                                                                          acmg73$Inheritance == "XL"]) %>% group_by(npm_research_id,
                                                                                                                    genetic_sex,
                                                                                                                    genetic_ethnicity,
                                                                                                                    gene_symbol,
                                                                                                                    genotype_code) %>% summarise(Count = n())
# individual-level summary of PLP variant carriers by number of PLP variants per individual
indv_summary_count.per.indv.XL.v2 <- indv_summary.gene.XL.v2 %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity) %>% summarise(Count = n())

# number of carriers as tabulated in Supplementary Table 5: Feature = XL genes
# by genetic ancestry
indv_summary_count.per.indv.XL.v2 %>% group_by(genetic_ethnicity) %>% summarise(Count = n())
# by genetic sex
indv_summary_count.per.indv.XL.v2 %>% group_by(genetic_sex) %>% summarise(Count = n())
# by genetic ancestry and genetic sex
indv_summary_count.per.indv.XL.v2 %>% group_by(genetic_ethnicity, genetic_sex) %>% summarise(Count = n())


# To test for statistical significance (Fisher's exact test) of proportions of PLP carriers across ancestry (Supplementary Table 5)
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


# Section 3: Code for Figure1A and Extended Data Figure 2

# consolidate into a table the carrier counts of PLP variant carriers of AD genes by ancestry group
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


# plot adjusted carrier frequency of HBOC genes (BRCA1, BRCA2, PALB2) by variant and ancestry group (Extended Data Figure 2)
PLP_var_HBOC <- read.table("Variant_Level_r5.3_20211117_HBOC.txt", header = TRUE, sep = "\t")
PLP_var_HBOC[PLP_var_HBOC == 0] <- NA
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


# To consolidate list of PLP AD condition variants:
PLP_AD_varcount <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117_varonly.txt", header = TRUE, sep = "\t")
# PLP variants in ACMG SF v2 list of AD genes
PLP_AD_varcount.acmg59 <- PLP_AD_varcount %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD" & acmg73$Version != 3] & !(hgvs_c %in% apoblof$hgvs_c))
# PLP variants in AD genes that are newly added to ACMG SF v3
PLP_AD_varcount.acmg73_only <- PLP_AD_varcount %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD" & acmg73$Version == 3] & !(hgvs_c %in% apoblof$hgvs_c))


