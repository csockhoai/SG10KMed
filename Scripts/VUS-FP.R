# This code is to perform analysis for identifying and validation of VUS-favour pathogenic (VUS-FP) variants
# Generated outcomes are displayed in Figures 2C-D
# This script has two sections:
# 1) To identify and consolidate list of VUS-FP variants
# 2) To compare LDL cholesterol levels across LDLR variant carrier status

# Data files required:
# Tier2HighREVEL.csv
# Tier2SpliceAI_edit.csv
# All_Tier1_Pathogenic_Variants_r5.3_20211117_varonly.txt
# VUS-LP_exclusionlist.txt
# 20210711_LDLR-boxplot_v2.csv
# Full_Sample_Info_9051_edit_v2.txt
# Z0_SG10K_GeneList.csv
# Z0_MANElist.csv
# Z0_ACMG73list.csv
# Z0_ClinVar2plus_mutspec.csv
# ExclusionList_01_APOB_LOFv.txt

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)

# ================================
# Pre-02 : Load data files
# ================================
# VUS with high REVEL score (>0.7)
tier2ms <- read.csv("Tier2HighREVEL.csv", header = TRUE)
# VUS with high SpliceAI score (>0.8)
tier2spl <- read.csv("Tier2SpliceAI_edit.csv", header = TRUE)
# exclusion list of VUS that are (1) expert panel reviewed benign, OR (2) gnomAD homozygous carriers > 0
vuslp_excl <- read.table("VUS-LP_exclusionlist.txt", header = TRUE, sep = "\t")
# all identified pathogenic/likely pathogenic (PLP) variants in SG10K_Health (Tier 1 variants)
PLP_var <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117_varonly.txt", header = TRUE, sep = "\t")
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
# list of genes curated for this study
genes <- read.csv("Z0_SG10K_GeneList.csv", header = TRUE)
# list of MANE transcripts for genes
mane <- read.csv("Z0_MANElist.csv", header = TRUE)
# ACMG SF version 3 (v.3) genes list
acmg <- read.csv("Z0_ACMG73list.csv", header = TRUE)
# gene list with summary of variant consequence types for ClinVar two-star pathogenic/likely pathogenic variants
mutspec <- read.csv("Z0_ClinVar2plus_mutspec.csv", header = TRUE)
# exclusion list of APOB loss-of-function variants: APOB LOF variants are not associated with autosomal dominant familial hypercholesterolemia
apoblof <- read.table("ExclusionList_01_APOB_LOFv.txt", header = TRUE, sep = "\t")
# individual-level LDL cholesterol levels and LDLR variant carrier status
ldlr2 <- read.csv("20210711_LDLR-boxplot_v2.csv", header = TRUE)

# =======================================================================================
# Section 1.0 : Identify variants of uncertain significance-favour pathogenic (VUS-FP)
# =======================================================================================
# Identify list of VUS with prediction scores meeting deleterious threshold (i.e. variants with REVEL > 0.7 or Splice_AI > 0.8) in SG10K_Med gene list & affecting MANE transcript
VUSms_1 <- tier2ms %>% filter(SYMBOL %in% genes$Gene & Feature %in% mane$Transcript) %>% filter(!(HGVSc %in% PLP_var$hgvs_c))
VUSspl_1 <- tier2spl %>% filter(SYMBOL %in% genes$Gene & Feature %in% mane$Transcript) %>% filter(!(HGVSc %in% PLP_var$hgvs_c))
VUSms_1b <- VUSms_1 %>% filter(!(HGVSc %in% VUSspl_1$HGVSc))

# Count the total number of (missense and splice) VUS with prediction scores meeting deleterious threshold
total_VUSms <- length(VUSms_1b$HGVSc)
total_VUSspl <- length(VUSspl_1$HGVSc)
total_VUS <- sum(total_VUSms, total_VUSspl)
total_VUS

# To identify and consolidate VUS-FP variants
# Select missense variants with more than 2 pathogenic/likely pathogenic variant in the 25bp rolling window
tier2ms_vuslp.a <- VUSms_1b %>% filter(WINDOW_PLP_COUNT > 2)
# Select splice variants that occur in genes with more than 4 ClinVar two-star pathogenic/likely pathogenic nonsense/frameshift/splice variants
tier2spl_vuslp.a <- VUSspl_1 %>% filter(SYMBOL %in% mutspec$Gene[mutspec$PLP_NFS > 4])

# Combine missense and splice VUS-FP into a single list
cols.num <- c("REVEL_score", "SpliceAI_max", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
tier2ms_vuslp.a[cols.num] <- sapply(tier2ms_vuslp.a[cols.num], as.numeric)
tier2spl_vuslp.a[cols.num] <- sapply(tier2spl_vuslp.a[cols.num], as.numeric)
tier2ms_vuslp.b <- tier2ms_vuslp.a %>% mutate(MutType = case_when(REVEL_score > 0.7 ~ "MS"))
tier2spl_vuslp.b <- tier2spl_vuslp.a %>% mutate(MutType = case_when(SpliceAI_max >= 0.8 ~ "PTV"))
tier2_vuslp.all <- union(tier2ms_vuslp.b, tier2spl_vuslp.b) %>% mutate(VarType = case_when((MutType == "MS") ~ "Tier2_MS", (MutType == "PTV") ~ "Tier2_PTV"))
#tier2_vuslp.all <- union(tier2ms_vuslp.a, tier2spl_vuslp.a)

# Exclude variants that have (1) >0 number of homozygous carriers in SG10K_Med & gnomAD_exome & gnomAD_genome, or (2) expert-panel reviewed Benign status in ClinVar
tier2_vuslp <- tier2_vuslp.all %>% filter(!HGVSc %in% vuslp_excl$HGVSc) %>% filter(!(nHomAlt > 0) & !(gnomADe_controls_nhomalt > 0) & !(gnomADg_3_nhomalt > 0))

# Count number of VUS-FP
length(tier2_vuslp$HGVSc)
# Count number of VUS-FP not in ClinVar
length(tier2_vuslp$HGVSc[tier2_vuslp$ClinVar_CLNSIG == "."])

# ==============================================================================================
# Section 1.1 : Identify VUS-FP variants and carriers in ACMG SF v3.0 dominant condition genes
# ==============================================================================================
# Subset VUS-FP in ACMG SF v3.0 genes list
tier2_vuslp.acmg73 <- tier2_vuslp %>% filter(SYMBOL %in% acmg$Gene)
# Count number of VUS-FP in ACMG SF v3.0 genes list
length(tier2_vuslp.acmg73$HGVSc)

# Subset VUS-FP in ACMG SF v3.0 genes predisposing to autosomal dominant conditions
tier2_vuslp.acmg73_AD <- tier2_vuslp %>% filter(SYMBOL %in% acmg$Gene[acmg$Inheritance == "AD"])
# Identify carriers of VUS-FP in ACMG SF v3.0 AD genes
tier2_vuslp.acmg73_AD_indv <- tier2_vuslp.acmg73_AD %>% select(HetSamples)
tier2_vuslp.acmg73_AD_indv.2 <- tier2_vuslp.acmg73_AD_indv %>% separate_rows(HetSamples, sep = ",") %>% distinct(HetSamples)

# To compare with ACMG SF v3.0 AD gene P/LP vs VUS-FP carriers and identify number of additional at-risk individuals
# Identify carriers of PLP variants in ACMG SF v3.0 AD genes
PLP_acmg73_AD <- PLP_indv %>% filter(gene_symbol %in% acmg$Gene[acmg$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c))
PLP_acmg73_AD_indv <- PLP_acmg73_AD %>% select(npm_research_id) %>% distinct(npm_research_id)
# Identify non-overlapping individuals carrying PLP, VUS-FP variants
VUSFP_acmg73_AD_newindv <- tier2_vuslp.acmg73_AD_indv.2 %>% filter(!(HetSamples %in% PLP_acmg73_AD_indv$npm_research_id))
length(VUSFP_acmg73_AD_newindv$HetSamples)

# ====================================================================================================
# Section 1.2 : Compare distribution of PLP vs VUS-FP variants in ACMG SF v3.0 genes list (Figure 2C)
# ====================================================================================================
# note: in this section of code, tier1 = PLP variants, tier2 = VUS-FP
# Consolidate PLP variants in ACMG SF v3.0 by variant type (missense (MS), PTV)
tier1_ms_acmg <- PLP_var %>% filter(gene_symbol %in% acmg$Gene)
tier1_ms_acmg_2 <-
  tier1_ms_acmg %>% mutate(MutType = case_when((
    consequence == "missense_variant" |
      consequence == "missense_variant&splice_region_variant"
  ) ~ "MS",
  (
    consequence == "frameshift_variant" |
      consequence == "splice_acceptor_variant" |
      consequence == "splice_donor_variant" |
      consequence == "splice_donor_variant&intron_variant" |
      consequence == "stop_gained"
  ) ~ "PTV",
  TRUE ~ "Other"
  ))
tier1_var <- tier1_ms_acmg_2 %>% mutate(VarType = case_when((MutType == "MS") ~ "Tier1_MS",
                                                            (MutType == "PTV") ~ "Tier1_PTV",
                                                            (MutType == "Other") ~ "Tier1_Oth"))

# To summarize variant count of PLP, VUS-FP by variant type (MS, PTV)
tier1_var.a <- tier1_var %>% select(gene_symbol, VarType)
colnames(tier1_var.a) <- c("Gene", "VarType")
tier1_var.a$Gene <- factor(tier1_var.a$Gene, levels = unique(tier1_var.a$Gene))
tier1_var.a$VarType <- as.factor(tier1_var.a$VarType)
tier1_var.b <- tier1_var.a %>% group_by(Gene, VarType) %>% count(VarType) %>% filter(!VarType %in% c("Tier1_Oth")) %>% pivot_wider(names_from = "VarType", values_from = "n")

tier2_var.a <- tier2_vuslp %>% select(SYMBOL, VarType)
colnames(tier2_var.a) <- c("Gene", "VarType")
tier2_var.a$Gene <- factor(tier2_var.a$Gene, levels = unique(tier2_var.a$Gene))
tier2_var.a$VarType <- as.factor(tier2_var.a$VarType)
tier2_var.b <- tier2_var.a %>% group_by(Gene, VarType) %>% count(VarType) %>% pivot_wider(names_from = "VarType", values_from = "n")

# Plot distribution of PLP and VUS-LP by gene by variant type for genes in ACMG SF v3.0
acmg_domain <- acmg %>% select(Gene, PhenotypeCat3)
t1t2_var <- full_join(tier1_var.b, tier2_var.b, by = "Gene") 
t1t2_var.acmg <- left_join(acmg_domain, t1t2_var, by = "Gene")
t1t2_var.acmg[is.na(t1t2_var.acmg)] <- 0

t1t2_var.acmg_plot <- t1t2_var.acmg %>% pivot_longer(cols = 3:6, names_to = "Tier", values_to = "VarCount") %>% filter(VarCount > 0) %>% arrange(match(PhenotypeCat3, c("Cancer", "Cardiovascular", "CardioMetabolic", "Metabolic", "Other_ADgen", "Other_ARgen")), desc(VarCount))
t1t2_var.acmg_plot$Tier <- as.factor(t1t2_var.acmg_plot$Tier)
t1t2_var.acmg_plot$Gene <- factor(t1t2_var.acmg_plot$Gene, levels = unique(t1t2_var.acmg_plot$Gene))
tiercolor <- c('#ca0020', '#f4a582', '#0571b0', '#92c5de')

t1t2_var_acmg_summary <- t1t2_var.acmg_plot %>% group_by(Gene) %>% summarise(VarSum = sum(VarCount))
t1t2_var.acmg_plot_v2 <- data.frame()
for(i in 1:nrow(t1t2_var.acmg_plot)){
  myrow <- t1t2_var.acmg_plot[i,]
  tempvarsum <- t1t2_var_acmg_summary$VarSum[t1t2_var_acmg_summary$Gene == myrow$Gene]
  myrow <- cbind.data.frame(myrow, tempvarsum)
  t1t2_var.acmg_plot_v2 <- rbind.data.frame(t1t2_var.acmg_plot_v2, myrow)
}
colnames(t1t2_var.acmg_plot_v2)[ncol(t1t2_var.acmg_plot_v2)] <- "VarSum"

t1t2_var.acmg_plot_v3 <- t1t2_var.acmg_plot_v2 %>% arrange(match(PhenotypeCat3, c("Cancer", "Cardiovascular", "CardioMetabolic", "Metabolic", "Other_ADgen", "Other_ARgen")), desc(VarSum))
myorder <- match(t1t2_var.acmg_plot_v3$PhenotypeCat3, rev(c("Cancer", "Cardiovascular", "CardioMetabolic", "Metabolic", "Other_ADgen", "Other_ARgen")))
t1t2_var.acmg_plot_v3$Gene <- factor(t1t2_var.acmg_plot_v3$Gene, levels = unique(t1t2_var.acmg_plot_v3$Gene[order(myorder, t1t2_var.acmg_plot_v3$VarSum, decreasing =T)]))

F2C <-
  ggplot(t1t2_var.acmg_plot_v3, aes(x = Gene, y = VarCount, fill = Tier)) + geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + scale_fill_manual(values = tiercolor) + theme(
    panel.background = element_rect(
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
    axis.text.x = element_text(
      size = rel(0.85), vjust = 0.6, hjust = 1,
      angle = 45
    )
  ) + ylab("No. variants")
F2C

#Legend labels are modified for consistency of terminology using Adobe Illustrator (Tier1 = P/LP, Tier2 = VUS-FP)
#Dot matrix for disease domain are generated in Adobe Illustrator

# ==================================================================================
# Section 2.0 : Compare LDL cholesterol levels by LDLR variant carrier status 
# ==================================================================================
# individual-level LDL cholesterol levels and LDLR variant carrier status
ldlr2 <- read.csv("20210711_LDLR-boxplot_v2.csv", header = TRUE)

# To generate boxplot displaying distribution of LDL cholesterol levels by LDLR variant carrier status (Figure 2D)
ldlr2$LDLR_tier <- factor(ldlr2$LDLR_tier, levels = unique(ldlr2$LDLR_tier))
ggplot(ldlr2, aes(x = LDLR_tier, y = LDL)) + geom_boxplot() + theme(panel.background = element_blank(), panel.border = element_blank(), axis.line.x = element_line(linetype = "solid", size = 0.5), axis.line.y = element_line(linetype = "solid", size = 0.5)) + geom_hline(yintercept = 4.1, linetype = "dashed", size = 0.5) + ylim(0,8) + ylab("LDL cholesterol level, mmol/L") + xlab("LDLR variant status")

# To perform linear and logistic regression analysis of LDLR carriers to adjust for ancestry, sex, age (LDL_cat = captured lipid medication data)
# Consolidate demographic data to individual-level LDL/LDLR data
ldlr2_withdem <- ldlr2 %>% left_join(dem, by = c("Sample" = "npm_research_id"))

# Perform linear regression modeling to predict LDL cholesterol levels from LDLR variant status, correcting for covariates
baselinemodel <- glm(LDL~genetic_ethnicity+genetic_sex+age+lipid_med, data = ldlr2_withdem)
summary(glm(LDL~LDLR_tier+genetic_ethnicity+genetic_sex+age+lipid_med, data = ldlr2_withdem))
ldlr_linear <- glm(LDL~LDLR_tier+genetic_ethnicity+genetic_sex+age+lipid_med, data = ldlr2_withdem)

# Perform binomial logistic regression predicting LDL cholesterol levels from LDLR variant status, correcting for covariates
ldlr_logistic <- (glm(relevel(as.factor(LDL_cat), ref = "normal") ~LDLR_tier+genetic_ethnicity+genetic_sex+age, data = ldlr2_withdem, family = "binomial"))
summary(ldlr_logistic)
# obtain odds ratio:
exp(coef(ldlr_logistic))
# Obtain 95% confidence interval
exp(confint(ldlr_logistic))

