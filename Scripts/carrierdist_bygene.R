# This code is to analyse carrier distribution of pathogenic/likely pathogenic (PLP) variant carriers in dominant and recessive conditions genes in SG10K_Health
# Generated outcomes are displayed in Figure 1B, Figure 1C, Supplementary Figure 3
# This script has three sections:
# 1) Distribution for dominant conditions
# 2) Distribution for recessive conditions
# 3) Batch effect analysis for concordance in carrier frequencies in AD and AR genes for a) 15X vs 30X samples, b) SG10K_Health vs gnomAD

# Data files required:
# Full_Sample_Info_9051_edit_v2.txt
# All_Tier1_Pathogenic_Variants_r5.3_20211117.txt
# Gene_Level_r5.3_20211117.txt
# Gene_Level_r5.3_20211117_WithPValues.txt
# Z0_ACMG73list.csv
# Z0_ACMGcarrierpanel_20211008.txt
# Z0_CommercialCarrierPanel_20211008.txt
# Z0_Gene_DominantSet_exclusionlist_20211110.txt
# Z0_GeneDiseaseAssoc_v2_20211110.txt
# Z0_GeneDiseaseDomain_list.txt
# Z0_Mackenzie_Genes.txt
# gnomAD_T1_AD-ARgenes_20211219_corr_v2.csv

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)
library(gplots)

# ================================
# Pre-02 : Load data files
# ================================
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
# individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")
# summary of gene-level carrier frequencies by ancestry group
PLP_gene <- read.table("Gene_Level_r5.3_20211117.txt", header = TRUE, sep = "\t")
# summary of gene-level carrier frequencies by ancestry group with adjusted p-values for pairwise comparison across ancestry groups
PLP_Rgene <- read.table("Gene_Level_r5.3_20211117_WithPValues.txt", header = TRUE, sep = "\t")
# summary of variant-level carrier frequencies by ancestry group
PLP_var <- read.table("Variant_Level_r5.3_20211117.txt", sep = "\t", header = TRUE)
# ACMG SF version 3 (v3.0) genes list
acmg73 <- read.csv("Z0_ACMG73list.csv", header = TRUE)
# gene list with mode of inheritance
genelist_moi <- read.table("Z0_GeneDiseaseAssoc_v2_20211110.txt", header = TRUE, sep = "\t")
# gene list with disease domain category for Figure 1B
domain_cat <- read.table("Z0_GeneDiseaseDomain_list.txt", header = TRUE, sep = "\t")
# exclusion list of AD genes: variants in these genes are not associated with AD conditions or phenotype
GeneAD_excl <- read.table("Z0_Gene_DominantSet_exclusionlist_20211110.txt", header = TRUE, sep = "\t")
# list of genes recommended by ACMG for carrier screening
ACMGcarrier <- read.table("Z0_ACMGcarrierpanel_20211008.txt", header = TRUE, sep = "\t")
# list of genes in current commercial carrier screening gene panels
commcarrier <- read.table("Z0_CommercialCarrierPanel_20211008.txt", header = TRUE, sep = "\t")
# list of curated genes associated with severe recessive disorders (reference PMID 32678339)
mackenzie <- read.table("Z0_Mackenzie_Genes.txt", header = TRUE, sep = "\t")
# gene-level carrier frequencies of top 10 dominant and recessive genes (Table 1) for SG10K_Health and gnomAD (EAS, SAS) populations
cor1 <- read.csv("gnomAD_T1_AD-ARgenes_20211219_corr_v2.csv", header = TRUE)

# =================================
# Pre-03 : Key variables
# =================================
# consolidate individual demographic and variant data
indv <- left_join(dem, PLP_indv)

# =========================================================================
# Section 1.0 : Analysis of carrier distribution for dominant conditions
# =========================================================================

# expand gene-level summary dataset to include a) maximum ancestry group carrier frequency, and b) ACMG SF v3 list membership
PLP_gene.a <- transform(PLP_gene, max_frac = pmax(C_Frac, I_Frac, M_Frac))
PLP_gene.b <- PLP_gene.a %>% mutate(ACMG73 = case_when((PLP_gene.a$gene_symbol %in% acmg73$Gene) ~ "Y", TRUE ~ "N"))

# subset genes that are associated with AD conditions AND 
# (1) in ACMG SF v3.0 list, or 
# (2) has a maximum ancestry group carrier frequency of 0.5% or greater
df_AD <-
  PLP_gene.b %>% filter((
    gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance == "AD" | genelist_moi$Inheritance == "AD, AR"] &
      !(gene_symbol %in% GeneAD_excl$gene_symbol) &
      (gene_symbol %in% acmg73$Gene)
  ) |
    (
      gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance == "AD" | genelist_moi$Inheritance == "AD, AR"] &
        !(gene_symbol %in% GeneAD_excl$gene_symbol) & max_frac >= 0.005
    )
  )
df_AD

# =========================================================================================================
# Section 1.1 : Consolidate top 10 AD genes with highest carrier frequency from each ancestry (Table 1)
# =========================================================================================================

df_AD_CH <- df_AD %>% slice_max(C_Frac, n = 10)
df_AD_IND <- df_AD %>% slice_max(I_Frac, n = 10)
df_AD_MY <- df_AD %>% slice_max(M_Frac, n = 10)
df_AD_top10 <- rbind.data.frame(df_AD_CH, df_AD_IND, df_AD_MY) %>% distinct() %>% arrange(desc(C_Frac)) %>% rowwise() %>% mutate(C_percent = (C_Frac*100),
                                                                                                                                 I_percent = (I_Frac*100),
                                                                                                                                 M_percent = (M_Frac*100))
df_AD_top10_table1 <- df_AD_top10 %>% select(gene_symbol, C, C_percent, I, I_percent, M, M_percent, ACMG73)  
#write.table(df_AD_top10_table1, "output_carfreq_Table1_top10_AD.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# =============================================================================================================
# Section 1.2 : Plot adjusted carrier frequency for all AD genes in ACMG SF v3.0 list by ancestry (Figure 1B)
# =============================================================================================================
# prepare data table to plot Figure 1B (adjusted carrier frequency of AD genes by ancestry group and disease domain)
# assign genes to domain category
df_AD.a <- left_join(df_AD, domain_cat) %>% arrange(Domain)

# calculate carrier frequency (per 100000) adjusted by each ancestry group size 
df_AD.b <- df_AD.a %>% mutate(CH = (C_Frac*100000)/3, IND = (I_Frac*100000)/3, MY = (M_Frac*100000)/3) %>% rowwise() %>% mutate(Total = sum(CH, IND, MY))

# sort by highest total adjusted carrier frequency
domains <- unique(df_AD.b$Domain)
domains_customsort <- c()
for(i in 1:length(domains)){
  domains_sub <- df_AD.b[which(df_AD.b$Domain == domains[i]),]
  domains_customsort <- rbind(domains_customsort, domains_sub[order(domains_sub$Total, decreasing = T),])
}
df_AD.c <- pivot_longer(domains_customsort, cols = 12:14, names_to = "Ancestry", values_to = "CarFreq100K")
df_AD.c$Ancestry <- as.factor(df_AD.c$Ancestry)
df_AD.c$gene_symbol <- factor(df_AD.c$gene_symbol, levels = rev(unique(df_AD.c$gene_symbol)))

# genes with total carrier frequency <1% are separated from genes with carrier frequency >1% into different figure sub-panel for better resolution
# plot AD genes with high carrier frequencies (greater than 1%) in Figure 1B top panel
df_AD.high <- filter(df_AD.c, Domain %in% c("1_CF_High"))
mycolor <- c('#d73027', '#4575b4', '#fdae61')
F1Btop <- ggplot(df_AD.high, aes(x = gene_symbol, y = CarFreq100K, fill = Ancestry)) + geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + scale_fill_manual(values = mycolor) + theme(panel.background = element_rect(
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
axis.text.y = element_text(vjust = 0.2, hjust = 0.95, size = rel(0.85), face = "italic")
) + ylab("Adjusted carrier frequency (per 100,000)") + coord_flip()
F1Btop

# plot remaining AD genes in Figure 1B lower panel
df_AD.all <- filter(df_AD.c,!Domain %in% c("1_CF_High"))
mycolor <- c('#d73027', '#4575b4', '#fdae61')
F1Blow <- ggplot(df_AD.all, aes(x = gene_symbol, y = CarFreq100K, fill = Ancestry)) + geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + scale_fill_manual(values = mycolor) + theme(panel.background = element_rect(
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
axis.text.y = element_text(vjust = 0.2, hjust = 0.95, size = rel(0.85), face = "italic")
) + ylab("Adjusted carrier frequency (per 100,000)") + coord_flip()
F1Blow

# top and bottom panels of Figure 1B were consolidated using Adobe Illustrator

# ===================================================================================================================
# Section 1.3 : Summary of PLP variant-level carrier frequency by ancestry for the entire cohort (9051 individuals)
# ===================================================================================================================
# subset demographics data for entire cohort
dem_subset <- dem %>% select(npm_research_id, genetic_ethnicity, age, study, target_depth)
dem_subset_sum <- dem_subset %>% group_by(genetic_ethnicity) %>% count(name = "Count")

# consolidate individual-level PLP data with demographics
PLP_indv_all_var <- PLP_indv %>% select(gene_symbol, hgvs_g, hgvs_c, hgvs_p, npm_research_id, genotype_code)
PLP_indv_all_var.a <- left_join(PLP_indv_all_var, dem_subset, by = "npm_research_id")

# generate allele frequency for each variant
PLP_indv_all_var.a.2 <- PLP_indv_all_var.a %>% mutate(allele_count = case_when(genotype_code == "1" | genotype_code == "3" ~ 1, TRUE ~ 2))
PLP_var_ac_sum <- PLP_indv_all_var.a.2 %>% select(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity, allele_count) %>% group_by(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity) %>% summarise(across(everything(), sum)) %>% pivot_wider(names_from = genetic_ethnicity, values_from = allele_count) %>% replace(is.na(.),0) %>% mutate(
  Total_AC = (C + I + M),
  AF_CH = (C/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "C"]),
  AF_IND = (I/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "I"]),
  AF_MY = (M/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "M"]),
  AF_Total = (Total_AC/9051)
) %>% rowwise() %>% mutate(max_AF = pmax(AF_CH, AF_IND, AF_MY))

# calculate variant-level carrier frequency by ancestry
PLP_indv_all_var.b <-
  PLP_indv_all_var.a %>% select(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity) %>% group_by(gene_symbol, hgvs_g, hgvs_c, hgvs_p, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = genetic_ethnicity, values_from = n) %>% replace(is.na(.), 0) %>% mutate(
    Total = (C + I + M),
    C_Frac = (C/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "C"]),
    I_Frac = (I/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "I"]),
    M_Frac = (M/dem_subset_sum$Count[dem_subset_sum$genetic_ethnicity == "M"])
  ) %>% rowwise() %>% mutate(max_Frac = pmax(C_Frac, I_Frac, M_Frac))

# filter for only variants that have carrier frequency >= 0.1% in at least 1 ancestry group
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

# obtain carrier frequency for highlighted variants mentioned in Results
# variant carrier frequency for NOTCH3 Arg544Cys
NOTCH3_R544C <- pval_table_all_var %>% filter(hgvs_c == "NM_000435.3:c.1630C>T")
NOTCH3_R544C

# variant carrier frequency for PRSS1 Gly208Ala
PRSS1_G208A <- pval_table_all_var %>% filter(hgvs_c == "NM_002769.5:c.623G>C")
PRSS1_G208A

# variant carrier frequency for CTRC Ala73Thr
CTRC_A73T <- pval_table_all_var %>% filter(hgvs_c == "NM_007272.3:c.217G>A")
CTRC_A73T

# =========================================================================
# Section 2.0 : Analysis of carrier distribution for recessive conditions
# =========================================================================

# expand gene-level summary dataset (with calculated p-values) to include a) maximum ancestry group carrier frequency, and b) gene lists membership
PLP_Rgene.a <- transform(PLP_Rgene, max_frac = pmax(C_Frac, I_Frac, M_Frac))
PLP_Rgene.b <- PLP_Rgene.a %>% mutate(ACMG73 = case_when(PLP_Rgene.a$gene_symbol %in% acmg73$Gene ~ "Y", TRUE ~ "N"),
                                      ACMGcarrierpan = case_when(PLP_Rgene.a$gene_symbol %in% ACMGcarrier$OMIM.gene.name ~ "Y", TRUE ~ "N"),
                                      CommCarrier = case_when(PLP_Rgene.a$gene_symbol %in% commcarrier$Gene ~ "Y", TRUE ~ "N"),
                                      Mackenzie = case_when(PLP_Rgene.a$gene_symbol %in% mackenzie$Gene ~ "Y", TRUE ~ "N"))

# subset genes that are associated with recessive conditions AND 
# (1) maximum ancestry group carrier frequency of 0.5% or greater, AND 
# (2) adjusted p-value less than 0.05
df_AR <-
  PLP_Rgene.b %>% filter(gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance != "AD"] &
                           max_frac >= 0.005 & pval_list_adjusted_r5p3 < 0.05)

# ============================================================================================================================
# Section 2.1 : Plot AR genes with significantly different carrier frequencies across ancestry groups in heatmap (Figure 1C)
# ============================================================================================================================

recessive.1 <- df_AR %>% select(gene_symbol, C_Frac, I_Frac, M_Frac)
rownames(recessive.1) <- recessive.1$gene_symbol
recessive.1 <- recessive.1[,-1]
recessive.1.scaled <- t(scale(t(recessive.1)))
recessive.1_scaled <- cbind.data.frame(df_AR$ACMGcarrierpan, df_AR$gene_symbol, recessive.1.scaled)
colnames(recessive.1_scaled)[1:5] <- c("ACMGcarrier", "Gene", "CH", "IND", "MY")
lab <- recessive.1_scaled$ACMGcarrier
lab[lab == "N"] <- "Black"
lab[lab == "Y"] <- "Red"
labcols <- lab
mycolor <- colorRampPalette(c("yellow", "black", "blue"))
recessive.1_scaled_matrix <- recessive.1_scaled[, 3:ncol(recessive.1_scaled)]
heatmap.2(as.matrix(recessive.1_scaled_matrix), scale = "none", trace = "none", density.info = "none", col = mycolor(100), colRow = labcols, key = TRUE, keysize = 1)

# additional aesthetic edits (e.g. figure re-sizing, typography edits were done in Adobe Illustrator)

# =========================================================================================================
# Section 2.2 : Consolidate top 10 AR genes with highest carrier frequency from each ancestry (Table 1)
# =========================================================================================================
# List of top 10 AR genes with highest carrier frequencies from each ancestry group
df_AR.2 <-
  PLP_Rgene.b %>% filter(gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance != "AD"] &
                           max_frac >= 0.005)

df_AR.2_CH <- df_AR.2 %>% slice_max(C_Frac, n = 10)
df_AR.2_IND <- df_AR.2 %>% slice_max(I_Frac, n = 10)
df_AR.2_MY <- df_AR.2 %>% slice_max(M_Frac, n = 10)
df_AR.2_top10 <-
  rbind.data.frame(df_AR.2_CH, df_AR.2_IND, df_AR.2_MY) %>% distinct() %>% select(gene_symbol, C, I, M, Total, C_Frac, I_Frac, M_Frac) %>% arrange(desc(C_Frac)) %>% rowwise() %>% mutate(
    C_percent = (C_Frac * 100),
    I_percent = (I_Frac *
                   100),
    M_percent = (M_Frac *
                   100)
  )
df_AR.2_top10_table1 <- df_AR.2_top10 %>% select(gene_symbol, C, C_percent, I, I_percent, M, M_percent)
write.table(df_AR.2_top10_table1, "output_carfreq_Table1_top10_AR.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# =================================================================================================================
# Section 2.3 : Identify high carrier frequency variants and genes of recessive conditions - mentioned in Results
# =================================================================================================================
# Identify genes with carrier frequencies exceeding 1%
df_AR.3 <- PLP_Rgene.b %>% filter(gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance != "AD"] & max_frac >= 0.01)

# Identify recurrent variants in genes with high carrier frequencies
# highest variant carrier frequencies for GJB2
GJB2_max <- pval_table_all_var %>% filter(gene_symbol == "GJB2" & max_Frac >= 0.01)
GJB2_max

# highest variant carrier frequencies for HFE
HFE_max <- pval_table_all_var %>% filter(gene_symbol == "HFE" & max_Frac >= 0.01)
HFE_max

# highest variant carrier frequencies for CFTR
CFTR_max <- pval_table_all_var %>% filter(gene_symbol == "CFTR" & max_Frac >= 0.01)
CFTR_max

# highest variant carrier frequencies for SLC25A13 in Chinese
SLC25A13_c_max <- pval_table_all_var %>% filter(gene_symbol == "SLC25A13") %>% slice_max(C_Frac)
SLC25A13_c_max

# highest variant carrier frequencies for GNE in Indian
GNE_i_max <- pval_table_all_var %>% filter(gene_symbol == "GNE") %>% slice_max(I_Frac)
GNE_i_max

# Variants of recessive genes with highest carrier frequency in Malay
var_AR_MY <- pval_table_all_var %>% filter(M_Frac >= 0.01 & C_Frac < 0.01 & I_Frac < 0.01) %>% arrange(desc(M_Frac))

# Variants of recessive genes with highest carrier frequency in Chinese
var_AR_CH <- pval_table_all_var %>% filter(C_Frac >= 0.01 & I_Frac < 0.01 & M_Frac < 0.01) %>% arrange(desc(C_Frac))

# Variants of recessive genes with highest carrier frequency in Indians
var_AR_IND <- pval_table_all_var %>% filter(I_Frac >= 0.01 & C_Frac < 0.01 & M_Frac < 0.01) %>% arrange(desc(I_Frac))

# =========================================================================================================
# Section 2.4 : List of recessive genes missed by carrier testing gene panels - mentioned in Results
# =========================================================================================================
# List of recessive genes with carrier frequency exceeding 0.5% and not in carrier testing panels
# number of recessive genes with carrier frequency > 0.5%
length(df_AR.2$gene_symbol)

# List of genes in ACMG recommendation for carrier screening
Rgenes_acmgcarrier <- df_AR.2 %>% filter(ACMGcarrierpan == "Y")
length(Rgenes_acmgcarrier$gene_symbol)

# List of genes not in ACMG recommendation for carrier screening but covered in commercial carrier testing panels
Rgenes_commcarrier <- df_AR.2 %>% filter(ACMGcarrierpan == "N" & CommCarrier == "Y")
length(Rgenes_commcarrier$gene_symbol)

# List of genes not in ACMG recommendation for carrier screening or commercial carrier testing panels
Rgenes_notcovered <- df_AR.2 %>% filter(ACMGcarrierpan == "N" & CommCarrier == "N")
length(Rgenes_notcovered$gene_symbol)

# List of severe recessive genes with carrier frequency exceeding 0.5% and not covered in carrier testing panels
sevRgenes_notcovered <- df_AR.2 %>% filter(ACMGcarrierpan == "N" & CommCarrier == "N" & Mackenzie == "Y")
sevRgenes_notcovered_count <- length(Rgenes_notcovered$gene_symbol)
sevRgenes_notcovered_count
#total number of severe recessive genes
sevRgenes_all <- df_AR.2 %>% filter(Mackenzie == "Y")
sevRgenes_all_count <- length(sevRgenes_all$gene_symbol)
sevRgenes_all_count
#fraction of missed Asian severe recessive genes
sevRgenes_notcovered_frac <- sevRgenes_notcovered_count/sevRgenes_all_count
sevRgenes_notcovered_frac

# =============================================================================================
# Section 3.0 : Batch effect analysis - concordance in carrier frequencies in AD and AR genes
# =============================================================================================
# This code block has two sub-sections: 
# A) to compare the concordance of carrier frequency (by ancestry group) between 15X and 30X sequenced samples (Supplementary Figure 5)
# B) to compare the concordance of carrier frequency of top 10 dominant and recessive condition genes (Table 1) between SG10K_Health and gnomAD, Chinese vs EAS and Indians vs SAS (Supplementary Figure 3)

# =============================================================================================================================================
# Section 3.1 : Batch effect analysis - concordance in carrier frequencies in AD and AR genes in 15X vs 30X samples (Supplementary Figure 5)
# =============================================================================================================================================
# subset individual demographics for 15X and 30X samples
dem_15x <- dem %>% filter(target_depth == "15x")
dem_30x <- dem %>% filter(target_depth == "30x")

# remove P/LP variants with maximal carrier frequency > 15%
PLP_var_b1 <- PLP_var %>% mutate(max_freq = pmax(C_Frac,I_Frac,M_Frac))
PLP_var_b1_more10percent <- PLP_var_b1 %>% filter(max_freq > 0.1)
PLP_indv_b2 <- PLP_indv %>% filter(!hgvs_c %in% PLP_var_b1_more10percent$hgvs_c)

# filter for variants in ACMG SF v3.0 only
PLP_indv_b2.acmg73 <- PLP_indv_b2 %>% filter(gene_symbol %in% acmg73$Gene)

# consolidate individual demographic and variant data : all PLP variants
indv_all <- left_join(dem, PLP_indv_b2)
indv_15x <- left_join(dem_15x, PLP_indv_b2)
indv_30x <- left_join(dem_30x, PLP_indv_b2)

# consolidate individual demographic and variant data : PLP variants in ACMG SF v3.0 genes only
indv_all.acmg <- left_join(dem, PLP_indv_b2.acmg73)
indv_15x.acmg <- left_join(dem_15x, PLP_indv_b2.acmg73)
indv_30x.acmg <- left_join(dem_30x, PLP_indv_b2.acmg73)

# summarize number of 15X vs 30X samples by ancestry group
dem_15x %>% group_by(genetic_ethnicity) %>% count()
dem_30x %>% group_by(genetic_ethnicity) %>% count()

# ===============================================================================================================================================
# Section 3.1.1 : Batch effect analysis - concordance in carrier frequencies in dominant genes for 15X vs 30X samples (Supplementary Figure 5A)
# ===============================================================================================================================================
# consolidate carrier frequencies by genetic ancestry for 15X samples
gene_15x <- indv_15x %>% group_by(gene_symbol, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = "genetic_ethnicity", values_from = "n")
col_order <- c("gene_symbol", "C", "I", "M")
gene_15x.a <- gene_15x[, col_order]
gene_15x.a <- as.data.frame(gene_15x.a)
gene_15x.a[is.na(gene_15x.a)] <- 0
gene_15x.a$C <- as.numeric(gene_15x.a$C)
gene_15x.a$I <- as.numeric(gene_15x.a$I)
gene_15x.a$M <- as.numeric(gene_15x.a$M)

gene_15x.b <- gene_15x.a %>% mutate(C_Frac = C/3927,
                                    I_Frac = I/1717,
                                    M_Frac = M/1257)

gene_15x.c <- transform(gene_15x.b, max_frac = pmax(C_Frac, I_Frac, M_Frac))
gene_15x.d <- gene_15x.c %>% mutate(ACMG73 = case_when((gene_15x.c$gene_symbol %in% acmg73$Gene) ~ "Y", TRUE ~ "N"))
df_AD_15x <- gene_15x.d %>% filter((gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance == "AD" | genelist_moi$Inheritance == "AD, AR"] & !(gene_symbol %in% GeneAD_excl$gene_symbol) & (gene_symbol %in% acmg73$Gene)) | (gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance == "AD" | genelist_moi$Inheritance == "AD, AR"] & !(gene_symbol %in% GeneAD_excl$gene_symbol) & max_frac >= 0.005))

# consolidate carrier frequencies by genetic ancestry for 30X samples
gene_30x <- indv_30x %>% group_by(gene_symbol, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = "genetic_ethnicity", values_from = "n")
col_order <- c("gene_symbol", "C", "I", "M")
gene_30x.a <- gene_30x[, col_order]
gene_30x.a <- as.data.frame(gene_30x.a)
gene_30x.a[is.na(gene_30x.a)] <- 0
gene_30x.a$C <- as.numeric(gene_30x.a$C)
gene_30x.a$I <- as.numeric(gene_30x.a$I)
gene_30x.a$M <- as.numeric(gene_30x.a$M)

gene_30x.b <- gene_30x.a %>% mutate(C_Frac = C/1575,
                                    I_Frac = I/224,
                                    M_Frac = M/351)

gene_30x.c <- transform(gene_30x.b, max_frac = pmax(C_Frac, I_Frac, M_Frac))
gene_30x.d <- gene_30x.c %>% mutate(ACMG73 = case_when((gene_30x.c$gene_symbol %in% acmg73$Gene) ~ "Y", TRUE ~ "N"))
df_AD_30x <- gene_30x.d %>% filter((gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance == "AD" | genelist_moi$Inheritance == "AD, AR"] & !(gene_symbol %in% GeneAD_excl$gene_symbol) & (gene_symbol %in% acmg73$Gene)) | (gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance == "AD" | genelist_moi$Inheritance == "AD, AR"] & !(gene_symbol %in% GeneAD_excl$gene_symbol) & max_frac >= 0.005))

# generate scatterplot and perform correlation test of carrier frequency for 15X vs 30X samples by ancestry group
# for Chinese:
gene_30x_subset.CH <- df_AD_30x %>% select(gene_symbol, C_Frac)
gene_15x_subset.CH <- df_AD_15x %>% select(gene_symbol, C_Frac) %>% rename(C_Frac_15x = C_Frac)
gene_AD_15x30x.CH <- left_join(gene_15x_subset.CH, gene_30x_subset.CH) %>% replace(is.na(.),0)

ggplot(gene_AD_15x30x.CH, aes(x = C_Frac_15x, y = C_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("CH carrier frequency (15X samples, AD genes)") + ylab("CH carrier frequency (30X samples)")
cor.test(gene_AD_15x30x.CH$C_Frac_15x, gene_AD_15x30x.CH$C_Frac)
cor.test(gene_AD_15x30x.CH$C_Frac_15x, gene_AD_15x30x.CH$C_Frac)$p.value

# for Indian:
gene_30x_subset.IND <- df_AD_30x %>% select(gene_symbol, I_Frac)
gene_15x_subset.IND <- df_AD_15x %>% select(gene_symbol, I_Frac) %>% rename(I_Frac_15x = I_Frac)
gene_AD_15x30x.IND <- left_join(gene_15x_subset.IND, gene_30x_subset.IND) %>% replace(is.na(.),0)

ggplot(gene_AD_15x30x.IND, aes(x = I_Frac_15x, y = I_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("IND carrier frequency (15X samples, AD genes)") + ylab("IND carrier frequency (30X samples)")
cor.test(gene_AD_15x30x.IND$I_Frac_15x, gene_AD_15x30x.IND$I_Frac)
cor.test(gene_AD_15x30x.IND$I_Frac_15x, gene_AD_15x30x.IND$I_Frac)$p.value

# for Malay:
gene_30x_subset.MY <- df_AD_30x %>% select(gene_symbol, M_Frac)
gene_15x_subset.MY <- df_AD_15x %>% select(gene_symbol, M_Frac) %>% rename(M_Frac_15x = M_Frac)
gene_AD_15x30x.MY <- left_join(gene_15x_subset.MY, gene_30x_subset.MY) %>% replace(is.na(.),0)

ggplot(gene_AD_15x30x.MY, aes(x = M_Frac_15x, y = M_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("MY carrier frequency (15X samples, AD genes)") + ylab("MY carrier frequency (30X samples)")
cor.test(gene_AD_15x30x.MY$M_Frac_15x, gene_AD_15x30x.MY$M_Frac)
cor.test(gene_AD_15x30x.MY$M_Frac_15x, gene_AD_15x30x.MY$M_Frac)$p.value

# ===============================================================================================================================================
# Section 3.1.2 : Batch effect analysis - concordance in carrier frequencies in recessive genes for 15X vs 30X samples (Supplementary Figure 5B)
# ===============================================================================================================================================
# consolidate carrier frequencies by genetic ancestry for 15X samples
gene_15x <- indv_15x %>% group_by(gene_symbol, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = "genetic_ethnicity", values_from = "n")
col_order <- c("gene_symbol", "C", "I", "M")
gene_15x.a <- gene_15x[, col_order]
gene_15x.a <- as.data.frame(gene_15x.a)
gene_15x.a[is.na(gene_15x.a)] <- 0
gene_15x.a$C <- as.numeric(gene_15x.a$C)
gene_15x.a$I <- as.numeric(gene_15x.a$I)
gene_15x.a$M <- as.numeric(gene_15x.a$M)

gene_15x.b <- gene_15x.a %>% mutate(C_Frac = C/3927,
                                    I_Frac = I/1717,
                                    M_Frac = M/1257)

gene_15x.c <- transform(gene_15x.b, max_frac = pmax(C_Frac, I_Frac, M_Frac))
gene_15x.R1 <- gene_15x.c %>% mutate(ACMG73 = case_when(gene_15x.c$gene_symbol %in% acmg73$Gene ~ "Y", TRUE ~ "N"),
                                     ACMGcarrierpan = case_when(gene_15x.c$gene_symbol %in% ACMGcarrier$OMIM.gene.name ~ "Y", TRUE ~ "N"),
                                     CommCarrier = case_when(gene_15x.c$gene_symbol %in% commcarrier$Gene ~ "Y", TRUE ~ "N"),
                                     Mackenzie = case_when(gene_15x.c$gene_symbol %in% mackenzie$Gene ~ "Y", TRUE ~ "N"))
df_AR_15x <- gene_15x.R1 %>% filter(gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance != "AD"] & max_frac >= 0.005)

# consolidate carrier frequencies by genetic ancestry for 30X samples
gene_30x <- indv_30x %>% group_by(gene_symbol, genetic_ethnicity) %>% count() %>% pivot_wider(names_from = "genetic_ethnicity", values_from = "n")
col_order <- c("gene_symbol", "C", "I", "M")
gene_30x.a <- gene_30x[, col_order]
gene_30x.a <- as.data.frame(gene_30x.a)
gene_30x.a[is.na(gene_30x.a)] <- 0
gene_30x.a$C <- as.numeric(gene_30x.a$C)
gene_30x.a$I <- as.numeric(gene_30x.a$I)
gene_30x.a$M <- as.numeric(gene_30x.a$M)

gene_30x.b <- gene_30x.a %>% mutate(C_Frac = C/1575,
                                    I_Frac = I/224,
                                    M_Frac = M/351)

gene_30x.c <- transform(gene_30x.b, max_frac = pmax(C_Frac, I_Frac, M_Frac))
gene_30x.R1 <- gene_30x.c %>% mutate(ACMG73 = case_when(gene_30x.c$gene_symbol %in% acmg73$Gene ~ "Y", TRUE ~ "N"),
                                     ACMGcarrierpan = case_when(gene_30x.c$gene_symbol %in% ACMGcarrier$OMIM.gene.name ~ "Y", TRUE ~ "N"),
                                     CommCarrier = case_when(gene_30x.c$gene_symbol %in% commcarrier$Gene ~ "Y", TRUE ~ "N"),
                                     Mackenzie = case_when(gene_30x.c$gene_symbol %in% mackenzie$Gene ~ "Y", TRUE ~ "N"))
df_AR_30x <- gene_30x.R1 %>% filter(gene_symbol %in% genelist_moi$Gene[genelist_moi$Inheritance != "AD"] & max_frac >= 0.005)

# generate scatterplot and perform correlation test of carrier frequency for 15X vs 30X samples by ancestry group
# for Chinese:
gene_30x_Rsubset.CH <- df_AR_30x %>% select(gene_symbol, C_Frac)
gene_15x_Rsubset.CH <- df_AR_15x %>% select(gene_symbol, C_Frac) %>% rename(C_Frac_15x = C_Frac)
gene_R_15x30x.CH <- left_join(gene_15x_Rsubset.CH, gene_30x_Rsubset.CH) %>% replace(is.na(.),0)

ggplot(gene_R_15x30x.CH, aes(x = C_Frac_15x, y = C_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("CH carrier frequency (15X samples, AR genes)") + ylab("CH carrier frequency (30X samples)")
cor.test(gene_R_15x30x.CH$C_Frac_15x, gene_R_15x30x.CH$C_Frac)
cor.test(gene_R_15x30x.CH$C_Frac_15x, gene_R_15x30x.CH$C_Frac)$p.value

# for Indian:
gene_30x_Rsubset.IND <- df_AR_30x %>% select(gene_symbol, I_Frac)
gene_15x_Rsubset.IND <- df_AR_15x %>% select(gene_symbol, I_Frac) %>% rename(I_Frac_15x = I_Frac)
gene_R_15x30x.IND <- left_join(gene_15x_Rsubset.IND, gene_30x_Rsubset.IND) %>% replace(is.na(.),0)

ggplot(gene_R_15x30x.IND, aes(x = I_Frac_15x, y = I_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("IND carrier frequency (15X samples, AR genes)") + ylab("IND carrier frequency (30X samples)")
cor.test(gene_R_15x30x.IND$I_Frac_15x, gene_R_15x30x.IND$I_Frac)
cor.test(gene_R_15x30x.IND$I_Frac_15x, gene_R_15x30x.IND$I_Frac)$p.value

# for Malay:
gene_30x_Rsubset.MY <- df_AR_30x %>% select(gene_symbol, M_Frac)
gene_15x_Rsubset.MY <- df_AR_15x %>% select(gene_symbol, M_Frac) %>% rename(M_Frac_15x = M_Frac)
gene_R_15x30x.MY <- left_join(gene_15x_Rsubset.MY, gene_30x_Rsubset.MY) %>% replace(is.na(.),0)

ggplot(gene_R_15x30x.MY, aes(x = M_Frac_15x, y = M_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("MY carrier frequency (15X samples, AR genes)") + ylab("MY carrier frequency (30X samples)")
cor.test(gene_R_15x30x.MY$M_Frac_15x, gene_R_15x30x.MY$M_Frac)
cor.test(gene_R_15x30x.MY$M_Frac_15x, gene_R_15x30x.MY$M_Frac)$p.value

# =========================================================================================================================================================
# Section 3.1.3 : [Extra] Batch effect analysis - concordance in carrier frequencies in AD and AR genes for 15X vs 30X samples (overall, not by ancestry)
# =========================================================================================================================================================
# analysis for ACMG SF v3.0 dominant condition gene variants only, comparing total carrier frequency for 15X vs 30X samples
df_AD_30x.b1 <- df_AD_30x %>% rowwise() %>% mutate(Total = (C+I+M), Total_frac = Total/2150)
df_AD_30x.b2 <- df_AD_30x.b1 %>% select(gene_symbol, Total_frac) %>% rename(Total_frac_30x = Total_frac)
df_AD_15x.b1 <- df_AD_15x %>% rowwise() %>% mutate(Total = (C+I+M), Total_frac = Total/6901)
df_AD_15x.b2 <- df_AD_15x.b1 %>% select(gene_symbol, Total_frac) %>% rename(Total_frac_15x = Total_frac)
df_AD_all <- left_join(df_AD_30x.b2, df_AD_15x.b2) %>% replace(is.na(.),0)

# generate scatterplot and perform correlation test
ggplot(df_AD_all, aes(x = Total_frac_15x, y = Total_frac_30x)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("Carrier frequency (15X samples, AD genes)") + ylab("Carrier frequency (30X samples)")
cor.test(df_AD_all$Total_frac_15x, df_AD_all$Total_frac_30x)

# analysis for ACMG SF v3.0 recessive condition gene variants only, comparing total carrier frequency for 15X vs 30X samples
df_AR_30x.b1 <- df_AR_30x %>% rowwise() %>% mutate(Total = (C+I+M), Total_frac = Total/2150)
df_AR_30x.b2 <- df_AR_30x.b1 %>% select(gene_symbol, Total_frac) %>% rename(Total_frac_30x = Total_frac)
df_AR_15x.b1 <- df_AR_15x %>% rowwise() %>% mutate(Total = (C+I+M), Total_frac = Total/6901)
df_AR_15x.b2 <- df_AR_15x.b1 %>% select(gene_symbol, Total_frac) %>% rename(Total_frac_15x = Total_frac)
df_AR_all <- left_join(df_AR_30x.b2, df_AR_15x.b2) %>% replace(is.na(.),0)

# generate scatterplot and perform correlation test
ggplot(df_AR_all, aes(x = Total_frac_15x, y = Total_frac_30x)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("Carrier frequency (15X samples, AR genes)") + ylab("Carrier frequency (30X samples)")
cor.test(df_AR_all$Total_frac_15x, df_AR_all$Total_frac_30x)

# analysis for ACMG SF v3.0 dominant condition gene variants with carrier frequency <5% only, comparing total carrier frequency for 15X vs 30X samples
df_AD_all_less5 <- df_AD_all %>% mutate(max_freq = pmax(Total_frac_30x, Total_frac_15x)) %>% filter(max_freq < 0.05)
ggplot(df_AD_all_less5, aes(x = Total_frac_15x, y = Total_frac_30x)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("Carrier frequency (15X samples, AD genes)") + ylab("Carrier frequency (30X samples)")
cor.test(df_AD_all_less5$Total_frac_15x, df_AD_all_less5$Total_frac_30x)

# analysis for ACMG SF v3.0 recessive condition gene variants with carrier frequency <5% only, comparing total carrier frequency for 15X vs 30X samples
df_AR_all_less5 <- df_AR_all %>% mutate(max_freq = pmax(Total_frac_30x, Total_frac_15x)) %>% filter(max_freq < 0.05)
ggplot(df_AR_all_less5, aes(x = Total_frac_15x, y = Total_frac_30x)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("Carrier frequency (15X samples, AR genes)") + ylab("Carrier frequency (30X samples)")
cor.test(df_AR_all_less5$Total_frac_15x, df_AR_all_less5$Total_frac_30x)

# ================================================================================================================================================
# Section 3.2 : Batch effect analysis - concordance in carrier frequencies in SG10K_Health (CH/IND) vs gnomAD (EAS/SAS) (Supplementary Figure 3)
# ================================================================================================================================================

# generate correlation plot and perform correlation test for AD and AR genes with CH/EAS/IND/SAS with max carrier frequency >= 0.5% (Supplementary Figure 3)
cor1_SG <- cor1 %>% select(MOI, Gene, CH, IND)
cor1_gn <- cor1 %>% select(MOI, Gene, EAS, SAS)
colnames(cor1_gn) <- c("MOI", "Gene", "CH", "IND")
cor2_SG <- cor1_SG %>% pivot_longer(cols = 3:4, names_to = "Ancestry", values_to = "CarrierFreq_SG")
cor2_gn <- cor1_gn %>% pivot_longer(cols = 3:4, names_to = "Ancestry", values_to = "CarrierFreq_gn")
cor3 <- full_join(cor2_SG, cor2_gn) %>% filter(!is.na(CarrierFreq_SG))

cor3$Ancestry <- as.factor(cor3$Ancestry)
ggplot(cor3, aes(x = CarrierFreq_SG, y = CarrierFreq_gn)) + geom_point(aes(color = Ancestry), size = 1.5) + scale_color_manual(values = c('#d73027', '#4575b4')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("SG10K carrier frequency") + ylab("gnomAD carrier frequency")
cor.test(cor3$CarrierFreq_SG, cor3$CarrierFreq_gn)

# generate correlation plot and perform correlation test for AD and AR genes with CH/EAS only
cor3_CH <- cor3 %>% filter(Ancestry == "CH")
cor3_CH$Ancestry <- as.factor(cor3_CH$Ancestry)
ggplot(cor3_CH, aes(x = CarrierFreq_SG, y = CarrierFreq_gn)) + geom_point(aes(color = Ancestry), size = 1.5) + scale_color_manual(values = c('#d73027', '#4575b4')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("SG10K carrier frequency") + ylab("gnomAD carrier frequency")
cor.test(cor3_CH$CarrierFreq_SG, cor3_CH$CarrierFreq_gn)

# generate correlation plot and perform correlation test for AD and AR genes with IND/SAS only
cor3_IND <- cor3 %>% filter(Ancestry == "IND")
cor3_IND$Ancestry <- as.factor(cor3_IND$Ancestry)
ggplot(cor3_IND, aes(x = CarrierFreq_SG, y = CarrierFreq_gn)) + geom_point(aes(color = Ancestry), size = 1.5) + scale_color_manual(values = c('#4575b4')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("SG10K carrier frequency") + ylab("gnomAD carrier frequency")
cor.test(cor3_IND$CarrierFreq_SG, cor3_IND$CarrierFreq_gn)


