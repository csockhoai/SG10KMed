# This code is to perform analysis for identifying and validation of VUS-favour pathogenic variants
# Generated outcomes are displayed in Figures 2C-D
# This script has two sections:
# 1) To identify and consolidate VUS-FP variant vs PLP
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

library(tidyverse)

# load data files:
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


# Section 1: To identify and consolidate VUS-FP variant vs PLP

# To identify and consolidate VUS-FP variant list
# select missense variants with more than 2 pathogenic/likely pathogenic variant in the 25bp rolling window
tier2ms_vuslp.a <- tier2ms %>% filter(WINDOW_PLP_COUNT > 2)
# select splice variants that occur in genes with more than 4 ClinVar two-star pathogenic/likely pathogenic nonsense/frameshift/splice variants
tier2spl_vuslp.a <- tier2spl %>% filter(SYMBOL %in% mutspec$Gene[mutspec$PLP_NFS > 4])

# combine missense and splice VUS-LP into a single list
cols.num <- c("REVEL_score", "SpliceAI_max", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
tier2ms_vuslp.a[cols.num] <- sapply(tier2ms_vuslp.a[cols.num], as.numeric)
tier2spl_vuslp.a[cols.num] <- sapply(tier2spl_vuslp.a[cols.num], as.numeric)
tier2ms_vuslp.b <- tier2ms_vuslp.a %>% mutate(MutType = case_when(REVEL_score > 0.7 ~ "MS"))
tier2spl_vuslp.b <- tier2spl_vuslp.a %>% mutate(MutType = case_when(SpliceAI_max >= 0.8 ~ "PTV"))
tier2_vuslp.all <- union(tier2ms_vuslp.b, tier2spl_vuslp.b) %>% mutate(VarType = case_when((MutType == "MS") ~ "Tier2_MS", (MutType == "PTV") ~ "Tier2_PTV"))
#tier2_vuslp.all <- union(tier2ms_vuslp.a, tier2spl_vuslp.a)

# select for variants that overlap our genes list and are on MANE transcript
tier2_vuslp <- tier2_vuslp.all %>% filter(SYMBOL %in% genes$Gene) %>% filter(Feature %in% mane$Transcript)
# further select for variants in ACMG SF v3.0 genes list
tier2_vuslp.acmg73 <- tier2_vuslp %>% filter(SYMBOL %in% acmg$Gene)
# exclude variants that are (1) expert panel reviewed benign, OR (2) gnomAD homozygous carriers > 0
tier2_var <- tier2_vuslp.acmg73 %>% filter(!HGVSc %in% vuslp_excl$HGVSc)

# prepare curated pathogenic/likely pathogenic (PLP) variants from this study for comparison with VUS-FP
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

# To count number of variants for each category (PTV/MS) in PLP and VUS-LP lists
tier2_var.a <- tier2_var %>% select(SYMBOL, VarType)
colnames(tier2_var.a) <- c("Gene", "VarType")
tier2_var.a$Gene <- factor(tier2_var.a$Gene, levels = unique(tier2_var.a$Gene))
tier2_var.a$VarType <- as.factor(tier2_var.a$VarType)
tier2_var.b <- tier2_var.a %>% group_by(Gene, VarType) %>% count(VarType) %>% pivot_wider(names_from = "VarType", values_from = "n")

tier1_var.a <- tier1_var %>% select(gene_symbol, VarType)
colnames(tier1_var.a) <- c("Gene", "VarType")
tier1_var.a$Gene <- factor(tier1_var.a$Gene, levels = unique(tier1_var.a$Gene))
tier1_var.a$VarType <- as.factor(tier1_var.a$VarType)
tier1_var.b <- tier1_var.a %>% group_by(Gene, VarType) %>% count(VarType) %>% filter(!VarType %in% c("Tier1_Oth")) %>% pivot_wider(names_from = "VarType", values_from = "n")

# To plot distribution of PLP and VUS-LP by gene for ACMG SF v3 list
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


# Section 2: To compare LDL cholesterol levels across LDLR variant carrier status

# individual-level LDL cholesterol levels and LDLR variant carrier status
ldlr2 <- read.csv("20210711_LDLR-boxplot_v2.csv", header = TRUE)
# To generate boxplot comparing LDL cholesterol levels across LDLR variant carrier status (Figure 2D)
ldlr2$LDLR_tier <- factor(ldlr2$LDLR_tier, levels = unique(ldlr2$LDLR_tier))
ggplot(ldlr2, aes(x = LDLR_tier, y = LDL)) + geom_boxplot() + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, linetype = "solid", size = 0.5)) + geom_hline(yintercept = 4.1, linetype = "dashed", size = 0.5) + ylim(0,8) + ylab("LDL cholesterol level, mmol/L") + xlab("LDLR variant status")

# To perform logistic regression analysis of LDLR carriers to adjust for ancestry, sex, age (LDL_cat captured lipid med data)
# append demographic data to individual-level LDL/LDLR data
ldlr2_withdem <- ldlr2 %>% left_join(dem, by = c("Sample" = "npm_research_id"))
summary(glm(LDL~LDLR_tier+genetic_ethnicity+genetic_sex+age+lipid_med, data = ldlr2_withdem))
baselinemodel <- glm(LDL~genetic_ethnicity+genetic_sex+age+lipid_med, data = ldlr2_withdem)
ldlr_logistic <- (glm(relevel(as.factor(LDL_cat), ref = "normal") ~LDLR_tier+genetic_ethnicity+genetic_sex+age, data = ldlr2_withdem, family = "binomial"))
summary(ldlr_logistic)
# obtain odds ratio:
exp(coef(ldlr_logistic))
# obtain 95% confidence interval:
exp(confint(ldlr_logistic))
