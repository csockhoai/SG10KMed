# This code is to analyse carrier distribution of pathogenic/likely pathogenic (PLP) variant carriers in dominant and recessive conditions genes in SG10K_Health
# Generated outcomes are displayed in Figure 1B, Figure 1C, Supplementary Figure 1
# This script has three sections:
# 1) Distribution for dominant conditions
# 2) Distribution for recessive conditions
# 3) QC analysis for top 10 dominant/recessive condition genes

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

library(tidyverse)
library(gplots)

# load data files:
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
# individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")
# consolidate individual demographic and variant data
indv <- left_join(dem, PLP_indv)
# summary of gene-level carrier frequencies by ancestry group
PLP_gene <- read.table("Gene_Level_r5.3_20211117.txt", header = TRUE, sep = "\t")
# summary of gene-level carrier frequencies by ancestry group with adjusted p-values for pairwise comparison across ancestry groups
PLP_Rgene <- read.table("Gene_Level_r5.3_20211117_WithPValues.txt", header = TRUE, sep = "\t")

# ACMG SF version 3 (v.3) genes list
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


# Section 1: Carrier distribution for dominant conditions

# include two additional info (maximum ancestry group carrier frequency, ACMG SF v3 list membership) to gene-level summary dataset
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


# Section 2: Carrier distribution for recessive conditions

# include two additional info (maximum ancestry group carrier frequency, gene lists memberships) to gene-level summary dataset with p-values
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

# plot AR genes with significantly different carrier frequencies by ancestry group in heatmap (Figure 1C)
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


# Section 3: QC analysis for top 10 dominant/recessive condition genes
# This code block has two sub-sections: 
# A) to compare the concordance of carrier frequency (by ancestry group) between 15X and 30X sequenced samples,
# B) to compare the concordance of carrier frequency of top 10 dominant and recessive condition genes (Table 1) between SG10K_Health and gnomAD, Chinese vs EAS and Indians vs SAS

# Sub-section A: compare 15X vs 30X samples

# separate 15X and 30X samples
dem_15x <- dem %>% filter(target_depth == "15x")
dem_30x <- dem %>% filter(target_depth == "30x")
# consolidate individual demographic and variant data
indv_all <- left_join(dem, PLP_indv)
indv_15x <- left_join(dem_15x, PLP_indv)
indv_30x <- left_join(dem_30x, PLP_indv)

# summary number of 15X vs 30X samples by ancestry group
dem_15x %>% group_by(genetic_ethnicity) %>% count()
dem_30x %>% group_by(genetic_ethnicity) %>% count()

# generate dataframe for comparison for DOMINANT condition genes in 15X vs 30X samples
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

# generate scatterplot of carrier frequency for top 20 dominant condition genes in 15X vs 30X samples by ancestry group
# for Chinese:
gene_15x_top20.CH <- df_AD_15x %>% slice_max(C_Frac, n=20)
gene_30x_subset.CH <- df_AD_30x %>% select(gene_symbol, C_Frac)
gene_15x_top20.CH.sub <- gene_15x_top20.CH %>% select(gene_symbol, C_Frac) %>% rename(C_Frac_15x = C_Frac)
gene_top20.CH <- left_join(gene_15x_top20.CH.sub, gene_30x_subset.CH) %>% replace(is.na(.),0)

ggplot(gene_top20.CH, aes(x = C_Frac_15x, y = C_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("CH carrier frequency (15X samples, top20 genes)") + ylab("CH carrier frequency (30X samples)")
cor.test(gene_top20.CH$C_Frac_15x, gene_top20.CH$C_Frac)

# for Indian:
gene_15x_top20.IND <- df_AD_15x %>% slice_max(I_Frac, n=20)
gene_30x_subset.IND <- df_AD_30x %>% select(gene_symbol, I_Frac)
gene_15x_top20.IND.sub <- gene_15x_top20.IND %>% select(gene_symbol, I_Frac) %>% rename(I_Frac_15x = I_Frac)
gene_top20.IND <- left_join(gene_15x_top20.IND.sub, gene_30x_subset.IND) %>% replace(is.na(.),0)

ggplot(gene_top20.IND, aes(x = I_Frac_15x, y = I_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("IND carrier frequency (15X samples, top20 genes)") + ylab("IND carrier frequency (30X samples)")
cor.test(gene_top20.IND$I_Frac_15x, gene_top20.IND$I_Frac)

# for Malay:
gene_15x_top20.MY <- df_AD_15x %>% slice_max(M_Frac, n=20)
gene_30x_subset.MY <- df_AD_30x %>% select(gene_symbol, M_Frac)
gene_15x_top20.MY.sub <- gene_15x_top20.MY %>% select(gene_symbol, M_Frac) %>% rename(M_Frac_15x = M_Frac)
gene_top20.MY <- left_join(gene_15x_top20.MY.sub, gene_30x_subset.MY) %>% replace(is.na(.),0)

ggplot(gene_top20.MY, aes(x = M_Frac_15x, y = M_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("MY carrier frequency (15X samples, top20 genes)") + ylab("MY carrier frequency (30X samples)")
cor.test(gene_top20.MY$M_Frac_15x, gene_top20.MY$M_Frac)


# generate dataframe for comparison for RECESSIVE condition genes in 15X vs 30X samples
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

# generate scatterplot of carrier frequency for top 20 recessive condition genes in 15X vs 30X samples by ancestry group
# for Chinese:
gene_15x_Rtop20.CH <- df_AR_15x %>% slice_max(C_Frac, n=20)
gene_30x_Rsubset.CH <- df_AR_30x %>% select(gene_symbol, C_Frac)
gene_15x_Rtop20.CH.sub <- gene_15x_Rtop20.CH %>% select(gene_symbol, C_Frac) %>% rename(C_Frac_15x = C_Frac)
gene_Rtop20.CH <- left_join(gene_15x_Rtop20.CH.sub, gene_30x_Rsubset.CH) %>% replace(is.na(.),0)

ggplot(gene_Rtop20.CH, aes(x = C_Frac_15x, y = C_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("CH carrier frequency (15X samples, top20 genes)") + ylab("CH carrier frequency (30X samples)")
cor.test(gene_Rtop20.CH$C_Frac_15x, gene_Rtop20.CH$C_Frac)

# for Indian:
gene_15x_Rtop20.IND <- df_AR_15x %>% slice_max(I_Frac, n=20)
gene_30x_Rsubset.IND <- df_AR_30x %>% select(gene_symbol, I_Frac)
gene_15x_Rtop20.IND.sub <- gene_15x_Rtop20.IND %>% select(gene_symbol, I_Frac) %>% rename(I_Frac_15x = I_Frac)
gene_Rtop20.IND <- left_join(gene_15x_Rtop20.IND.sub, gene_30x_Rsubset.IND) %>% replace(is.na(.),0)

ggplot(gene_Rtop20.IND, aes(x = I_Frac_15x, y = I_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("IND carrier frequency (15X samples, top20 genes)") + ylab("IND carrier frequency (30X samples)")
cor.test(gene_Rtop20.IND$I_Frac_15x, gene_Rtop20.IND$I_Frac)

# for Malay:
gene_15x_Rtop20.MY <- df_AR_15x %>% slice_max(M_Frac, n=20)
gene_30x_Rsubset.MY <- df_AR_30x %>% select(gene_symbol, M_Frac)
gene_15x_Rtop20.MY.sub <- gene_15x_Rtop20.MY %>% select(gene_symbol, M_Frac) %>% rename(M_Frac_15x = M_Frac)
gene_Rtop20.MY <- left_join(gene_15x_Rtop20.MY.sub, gene_30x_Rsubset.MY) %>% replace(is.na(.),0)

ggplot(gene_Rtop20.MY, aes(x = M_Frac_15x, y = M_Frac)) + geom_point() + geom_smooth(method = lm) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ))  + xlab("MY carrier frequency (15X samples, top20 genes)") + ylab("MY carrier frequency (30X samples)")
cor.test(gene_Rtop20.MY$M_Frac_15x, gene_Rtop20.MY$M_Frac)


# Sub-section B: comparing SG10K_Health (CH/IND) vs gnomAD (EAS/SAS)

# generate correlation plot for AD and AR genes with CH/EAS/IND/SAS max frequency >= 0.5% (Supplementary Figure 1)
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


# generate separate plots for comparing CH/EAS and IND/SAS
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



