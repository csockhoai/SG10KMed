# This code is to perform admixture analysis
# Generated outcomes are displayed in Figures 2A-B, Supplementary Figure 4, Supplementary Tables 4-5
# This script has four sections:
# 1) ADMIXTURE analysis
# 2) Comparing genetic ancestry with self-reported race/ethnicity (R/E)
# 3) Summary of carriers of ancestry-specific P/LP variants in R/E matched and R/E mismatched groups
# 4) Calculations and statistical tests for comparing ancestral components of carriers/non-carriers of ancestry-specific P/LP variants

# Data files required:
# Full_Sample_Info_9051_edit_v2.txt
# All_Tier1_Pathogenic_Variants_r5.3_20211117.txt
# SG10KMed.r5.3.fam
# SG10KMed.r5.3.pruned.3.Q
# Variant_Level_r5.3_20211117.txt
# Combined_Spec_Var_Local_Ancestry_Based_r5p3_20220704.txt

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)
library(patchwork)

# ================================
# Pre-02 : Load data files
# ================================
admixture_fam_5.3 <- read.table("SG10KMed.r5.3.fam", stringsAsFactors = F, sep = " ")
full_sample_info_9051 <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
admixture_table_5.3 <- read.table("SG10KMed.r5.3.pruned.3.Q")
combined_spec_var_r5p3_main <- read.table("Combined_Spec_Var_Local_Ancestry_Based_r5p3_20220704.txt", sep = "\t", header = TRUE)

# =======================================
# Section 1: ADMIXTURE analysis
# =======================================

# To generate ADMIXTURE plot Figure 2A
# load admixture sample list (to obtain order of samples for the Q values table that will be loaded next)
admixture_fam_5.3 <- read.table("SG10KMed.r5.3.fam", stringsAsFactors = F, sep = " ")
admixture_samplelist_5.3 <- admixture_fam_5.3$V1

full_sample_info_9051 <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
admixture_unrel_indx_5.3 <- which(admixture_samplelist_5.3 %in% full_sample_info_9051$npm_research_id)

# Assign colors for plot
eth_sub_col_5.3 <- rep("gray",dim(full_sample_info_9051)[1]) # black
eth_sub_col_5.3[full_sample_info_9051$supplied_ethnicity =="C"] <- "#BB5566"
eth_sub_col_5.3[full_sample_info_9051$supplied_ethnicity =="I"] <- "#004488"
eth_sub_col_5.3[full_sample_info_9051$supplied_ethnicity =="M"] <- "#DDAA33"

admixture_table_5.3 <- read.table("SG10KMed.r5.3.pruned.3.Q") #Load Q-values
admixture_table_5.3 <- admixture_table_5.3[admixture_unrel_indx_5.3,] #Subset to retain only unrelated individuals

# Cluster and reorder samples
admixture.rows.cor_5.3 <- cor(t(admixture_table_5.3), use = "pairwise.complete.obs", method = "pearson") #Calculate all pairwise correlation
admixture.rows.hc_5.3 <- hclust(as.dist(1-admixture.rows.cor_5.3)) #Apply hierarchical clustering 
clust_admixture_5.3 <- admixture.rows.hc_5.3$order #Get ordering of samples by position in hierarchical clustering tree
admixture_table_5.3 <- admixture_table_5.3[clust_admixture_5.3,] #Reorder samples in Q-values table according to hierarchical clustering

# Plot admixture values and sample R/E labels
rownames(admixture_table_5.3) = NULL
par(mfrow=c(2,1))
barplot(t(as.matrix(admixture_table_5.3)), col=c("#DDAA33","#004488","#BB5566"), border=NA,space = c(0), ylab = "Q")
barplot(rep(1,length(eth_sub_col_5.3)), col=eth_sub_col_5.3[clust_admixture_5.3],space = c(0), border=NA, yaxt='n')
par(mfrow=c(1,1))

# Merge Q-value table with sample info table
joined_admix_5.3 <- full_sample_info_9051 %>% select(npm_research_id, genetic_ethnicity, supplied_ethnicity)
joined_admix_5.3 <- cbind.data.frame(joined_admix_5.3[clust_admixture_5.3,], admixture_table_5.3)
colnames(joined_admix_5.3) <- c("NPM_Research_ID", "Genetic_Ethnicity", "Supplied_Ethnicity", "Malay_Comp", "Indian_Comp", "Chinese_Comp")
joined_admix_5.3 <- joined_admix_5.3 %>% rowwise() %>% mutate(maxQ = max(Malay_Comp, Indian_Comp, Chinese_Comp))

# ==================================================================================
# Section 2: Comparing genetic ancestry with self-reported race/ethnicity (R/E) 
# ==================================================================================

# To generate barplot for comparing maxQ levels of R/E mismatched and R/E matched individuals (Supplementary Figure 4)
maxq <- joined_admix_5.3 %>% rowwise() %>% mutate(RE_Status = case_when(Genetic_Ethnicity == Supplied_Ethnicity ~ "Match",
                                                                        TRUE ~ "Mismatch"))
maxq %>% group_by(RE_Status) %>% count()
maxq$RE_Status <- factor(maxq$RE_Status, levels = rev(unique(maxq$RE_Status)))
F9 <- ggplot(maxq, aes(x = RE_Status, y = maxQ, fill = RE_Status)) + geom_violin(trim = FALSE) + geom_boxplot(width=0.05) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, linetype = "solid", size = 0.5)) + xlab("R/E vs Genetic Ancestry Concordance")
F9

# To calculate the median maxQ of R/E mismatch and R/E match group, and statistically test the difference
maxq %>% group_by(RE_Status) %>% summarise(Median = median(maxQ))
wilcox.test(maxq$maxQ~maxq$RE_Status)$p.value

# To count R/E mismatch & match individuals and summarise descriptive statistics (Supplementary Table 4)
maxq %>% group_by(RE_Status, Supplied_Ethnicity, Genetic_Ethnicity) %>% count()
maxq %>% group_by(RE_Status, Supplied_Ethnicity, Genetic_Ethnicity) %>% summarise(Median_MY = median(Malay_Comp), 
                                                                                  Median_IND = median(Indian_Comp), 
                                                                                  Median_CH = median(Chinese_Comp),
                                                                                  Q1_MY = quantile(Malay_Comp, 0.25),
                                                                                  Q3_MY = quantile(Malay_Comp, 0.75),
                                                                                  Q1_IND = quantile(Indian_Comp, 0.25),
                                                                                  Q3_IND = quantile(Indian_Comp, 0.75),
                                                                                  Q1_CH = quantile(Chinese_Comp, 0.25),
                                                                                  Q3_CH = quantile(Chinese_Comp, 0.75))

# To calculate fraction of CH & IND individuals with >20% MY ancestral component (cryptic admixture in R/E-matched)
cryptic_admix.my <- joined_admix_5.3 %>% rowwise() %>% mutate(RE_status = case_when(Genetic_Ethnicity == Supplied_Ethnicity ~ "matched", Genetic_Ethnicity != Supplied_Ethnicity ~ "mismatched")) %>% filter(RE_status == "matched" & Malay_Comp > 0.2 & Genetic_Ethnicity != "M")
cryptic_admix_2.my <- cryptic_admix.my %>% group_by(Genetic_Ethnicity) %>% count() %>% rename(cryp.admix_indv = n)
cryptic_admix_2.my$total_indv <- c(5502, 1941)
cryptic_admix_2.my$frac_cryp.admix <- cryptic_admix_2.my$cryp.admix_indv/cryptic_admix_2.my$total_indv
cryptic_admix_2.my

cryptic_admix.ind <- joined_admix_5.3 %>% rowwise() %>% mutate(RE_status = case_when(Genetic_Ethnicity == Supplied_Ethnicity ~ "matched", Genetic_Ethnicity != Supplied_Ethnicity ~ "mismatched")) %>% filter(RE_status == "matched" & Indian_Comp > 0.2 & Genetic_Ethnicity != "I")
cryptic_admix_2.ind <- cryptic_admix.ind %>% group_by(Genetic_Ethnicity) %>% count() %>% rename(cryp.admix_indv = n)
cryptic_admix_2.ind$total_indv <- c(5502, 1608)
cryptic_admix_2.ind$frac_cryp.admix <- cryptic_admix_2.ind$cryp.admix_indv/cryptic_admix_2.ind$total_indv
cryptic_admix_2.ind

cryptic_admix.ch <- joined_admix_5.3 %>% rowwise() %>% mutate(RE_status = case_when(Genetic_Ethnicity == Supplied_Ethnicity ~ "matched", Genetic_Ethnicity != Supplied_Ethnicity ~ "mismatched")) %>% filter(RE_status == "matched" & Chinese_Comp > 0.2 & Genetic_Ethnicity != "C")
cryptic_admix_2.ch <- cryptic_admix.ch %>% group_by(Genetic_Ethnicity) %>% count() %>% rename(cryp.admix_indv = n)
cryptic_admix_2.ch$total_indv <- c(1941, 1608)
cryptic_admix_2.ch$frac_cryp.admix <- cryptic_admix_2.ch$cryp.admix_indv/cryptic_admix_2.ch$total_indv
cryptic_admix_2.ch

# =============================================================================================================
# Section 3: Summary of carriers of ancestry-specific P/LP variants in R/E matched and R/E mismatched groups 
# =============================================================================================================

# read file of ancestry-specific variants defined by local ancestry
combined_spec_var_r5p3_main <- read.table("Combined_Spec_Var_Local_Ancestry_Based_r5p3_20220704.txt", sep = "\t", header = TRUE)
combined_spec_var_r5p3 <- combined_spec_var_r5p3_main %>% rename(gene_symbol = Symbol)

# subset ancestry-specific variant by ancestry
c_spec_var_r5p3 <- combined_spec_var_r5p3 %>% filter(LA_C > 0)
i_spec_var_r5p3 <- combined_spec_var_r5p3 %>% filter(LA_I > 0)
m_spec_var_r5p3 <- combined_spec_var_r5p3 %>% filter(LA_M > 0)

#read file of individual-level data for every P/LP variant
all_tier1_pathogenic_variants <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")
all_tier1_pathogenic_variants <- all_tier1_pathogenic_variants %>% left_join(full_sample_info_9051)

# To generate discordant carrier list per ancestry group
# First we create a function that will return a list of individuals with genetic ethnicity differing from that of the R/E specific variants
# For example, given a table with Chinese-specific variants, we find all individuals of non-Chinese genetic ancestry
get_discord_list <- function(spectable, ge){
  temp_discord_list <- c()
  for( i in 1:nrow(spectable))
  {
    #Get individuals 
    tempres <- all_tier1_pathogenic_variants %>% filter(hgvs_c == spectable$hgvs_c[i] & genetic_ethnicity != ge) %>% select(npm_research_id)
    temp_discord_list <- c(temp_discord_list, tempres$npm_research_id)
  }
  temp_discord_list <- unique(sort(temp_discord_list))
  temp_discord_list
}

# Call this function, once for each ethnicity
c_discord_list <- get_discord_list(c_spec_var_r5p3, "C")
m_discord_list <- get_discord_list(m_spec_var_r5p3, "M")
i_discord_list <- get_discord_list(i_spec_var_r5p3, "I")

# To compare maxQ in concordant vs discordant carriers (analysis of all ancestries combined) (Supplementary Table 5)
joined_admix_5.3_csub <- joined_admix_5.3 %>% filter(Genetic_Ethnicity == "C")
ingla_c <- joined_admix_5.3_csub$NPM_Research_ID %in% c(i_discord_list,m_discord_list)

joined_admix_5.3_isub <- joined_admix_5.3 %>% filter(Genetic_Ethnicity == "I")
ingla_i <- joined_admix_5.3_isub$NPM_Research_ID %in% c(m_discord_list, c_discord_list)

joined_admix_5.3_msub <- joined_admix_5.3 %>% filter(Genetic_Ethnicity == "M")
ingla_m <- joined_admix_5.3_msub$NPM_Research_ID %in% c(c_discord_list, i_discord_list)

joined_admix_5.3_cim <- joined_admix_5.3 %>% filter(Genetic_Ethnicity %in% c('C','I','M'))
ingla_total <- joined_admix_5.3_cim$NPM_Research_ID %in% unique(sort(c(c_discord_list, i_discord_list, m_discord_list)))
ingla_total_indv <- joined_admix_5.3_cim[ingla_total,]

# Get ancestry-specific max-cutoff
cim_lowquantiles <- joined_admix_5.3_cim %>% group_by(Genetic_Ethnicity) %>% summarise(lowquantile = quantile(maxQ, 0.25))
c_lowquantile <- joined_admix_5.3_cim %>% filter(Genetic_Ethnicity == "C" & maxQ < cim_lowquantiles$lowquantile[cim_lowquantiles$Genetic_Ethnicity == "C"])
i_lowquantile <- joined_admix_5.3_cim %>% filter(Genetic_Ethnicity == "I" & maxQ < cim_lowquantiles$lowquantile[cim_lowquantiles$Genetic_Ethnicity == "I"])
m_lowquantile <- joined_admix_5.3_cim %>% filter(Genetic_Ethnicity == "M" & maxQ < cim_lowquantiles$lowquantile[cim_lowquantiles$Genetic_Ethnicity == "M"])
lowquantile_ids <- unique(sort(c(c_lowquantile$NPM_Research_ID, i_lowquantile$NPM_Research_ID, m_lowquantile$NPM_Research_ID)))
lowquantile_status <- joined_admix_5.3_cim$NPM_Research_ID %in% lowquantile_ids
c_indx <- which(joined_admix_5.3_cim$Genetic_Ethnicity == "C")
i_indx <- which(joined_admix_5.3_cim$Genetic_Ethnicity == "I")
m_indx <- which(joined_admix_5.3_cim$Genetic_Ethnicity == "M")

# Calculate odds ratio of admixed vs non-admixed individuals carrying discordant P/LP variant and the adjusted (Benjamini-Hochberg) p-values
fisher.test(table(ingla_total[c_indx], lowquantile_status[c_indx]))
fisher.test(table(ingla_total[i_indx], lowquantile_status[i_indx]))
fisher.test(table(ingla_total[m_indx], lowquantile_status[m_indx]))
p.adjust(c(fisher.test(table(ingla_total[c_indx], lowquantile_status[c_indx]))$p.value, fisher.test(table(ingla_total[i_indx], lowquantile_status[i_indx]))$p.value, fisher.test(table(ingla_total[m_indx], lowquantile_status[m_indx]))$p.value), method = "BH")
#Pooled analysis
fisher.test(table(ingla_total, lowquantile_status))

# ==============================================================================================================================================
# Section 4: Calculations and statistical tests for comparing ancestral components of carriers/non-carriers of ancestry-specific P/LP variants
# ==============================================================================================================================================

# To calculate the fraction of ancestral component for discordant PLP carriers from each ancestry group, for all three ancestry-specific variant categories
# e.g. fraction of Chinese ancestral component in Malay and Indian individuals who are carrier or non-carrier for Chinese-specific PLP variants

# Malay carriers of c_spec_vars
c_spec_m_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% c_spec_var_r5p3$hgvs_c & genetic_ethnicity == "M")
c_spec_m_carriers_status <- joined_admix_5.3_msub$NPM_Research_ID %in% c_spec_m_carriers$npm_research_id
wilcox.test(joined_admix_5.3_msub$Chinese_Comp~c_spec_m_carriers_status)
boxplot(joined_admix_5.3_msub$Chinese_Comp~c_spec_m_carriers_status)

# Malay carriers of i_spec_vars
i_spec_m_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% i_spec_var_r5p3$hgvs_c & genetic_ethnicity == "M")
i_spec_m_carriers_status <- joined_admix_5.3_msub$NPM_Research_ID %in% i_spec_m_carriers$npm_research_id
wilcox.test(joined_admix_5.3_msub$Indian_Comp~i_spec_m_carriers_status)
boxplot(joined_admix_5.3_msub$Indian_Comp~i_spec_m_carriers_status)

# Indian carriers of c_spec_vars
c_spec_i_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% c_spec_var_r5p3$hgvs_c & genetic_ethnicity == "I")
c_spec_i_carriers_status <- joined_admix_5.3_isub$NPM_Research_ID %in% c_spec_i_carriers$npm_research_id
wilcox.test(joined_admix_5.3_isub$Chinese_Comp~c_spec_i_carriers_status)
boxplot(joined_admix_5.3_isub$Chinese_Comp~c_spec_i_carriers_status)

# Indian carriers of m_spec_vars
m_spec_i_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% m_spec_var_r5p3$hgvs_c & genetic_ethnicity == "I")
m_spec_i_carriers_status <- joined_admix_5.3_isub$NPM_Research_ID %in% m_spec_i_carriers$npm_research_id
wilcox.test(joined_admix_5.3_isub$Malay_Comp~m_spec_i_carriers_status)
boxplot(joined_admix_5.3_isub$Malay_Comp~m_spec_i_carriers_status)

# Chinese carriers of i_spec_vars
i_spec_c_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% i_spec_var_r5p3$hgvs_c & genetic_ethnicity == "C")
i_spec_c_carriers_status <- joined_admix_5.3_csub$NPM_Research_ID %in% i_spec_c_carriers$npm_research_id
wilcox.test(joined_admix_5.3_csub$Indian_Comp~i_spec_c_carriers_status)
boxplot(joined_admix_5.3_csub$Indian_Comp~i_spec_c_carriers_status)

# Chinese carriers of m_spec_vars
m_spec_c_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% m_spec_var_r5p3$hgvs_c & genetic_ethnicity == "C")
m_spec_c_carriers_status <- joined_admix_5.3_csub$NPM_Research_ID %in% m_spec_c_carriers$npm_research_id
wilcox.test(joined_admix_5.3_csub$Malay_Comp~m_spec_c_carriers_status)
boxplot(joined_admix_5.3_csub$Malay_Comp~m_spec_c_carriers_status)

# Others
c_spec_c_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% c_spec_var_r5p3$hgvs_c & genetic_ethnicity == "C")
c_spec_c_carriers_status <- joined_admix_5.3_csub$NPM_Research_ID %in% c_spec_c_carriers$npm_research_id
wilcox.test(joined_admix_5.3_csub$Chinese_Comp~c_spec_c_carriers_status)
m_spec_m_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% m_spec_var_r5p3$hgvs_c & genetic_ethnicity == "M")
m_spec_m_carriers_status <- joined_admix_5.3_msub$NPM_Research_ID %in% m_spec_m_carriers$npm_research_id
wilcox.test(joined_admix_5.3_msub$Malay_Comp~m_spec_m_carriers_status)
i_spec_i_carriers <- all_tier1_pathogenic_variants %>% filter(hgvs_c %in% i_spec_var_r5p3$hgvs_c & genetic_ethnicity == "I")
i_spec_i_carriers_status <- joined_admix_5.3_isub$NPM_Research_ID %in% i_spec_i_carriers$npm_research_id
wilcox.test(joined_admix_5.3_isub$Indian_Comp~i_spec_i_carriers_status)

# To tabulate median values of ancestral components calculated above (Supplementary Table 5)
spec_plot_df <- data.frame()
spec_plot_colnames <- c("NPM_Research_ID", "Ancestral_Fraction", "Carrier_Status", "Specific_Ancestry", "Genetic_Ancestry")

# C_Spec_C
c_c_temp <- joined_admix_5.3_csub %>% select(NPM_Research_ID, Chinese_Comp)
c_c_temp <- cbind.data.frame(c_c_temp, c_spec_c_carriers_status, rep("Chinese_Specific", nrow(c_c_temp)), rep("Chinese", nrow(c_c_temp)))
colnames(c_c_temp) <- spec_plot_colnames
c_c_carrier_median <- median(c_c_temp$Ancestral_Fraction[c_c_temp$Carrier_Status == TRUE])
c_c_noncarrier_median <- median(c_c_temp$Ancestral_Fraction[c_c_temp$Carrier_Status == FALSE])
c_c_carrier_mean <- mean(c_c_temp$Ancestral_Fraction[c_c_temp$Carrier_Status == TRUE])
c_c_noncarrier_mean <- mean(c_c_temp$Ancestral_Fraction[c_c_temp$Carrier_Status == FALSE])

# C_Spec_I
c_i_temp <- joined_admix_5.3_isub %>% select(NPM_Research_ID, Chinese_Comp)
c_i_temp <- cbind.data.frame(c_i_temp, c_spec_i_carriers_status, rep("Chinese_Specific", nrow(c_i_temp)), rep("Indian", nrow(c_i_temp)))
colnames(c_i_temp) <- spec_plot_colnames
c_i_carrier_median <- median(c_i_temp$Ancestral_Fraction[c_i_temp$Carrier_Status == TRUE])
c_i_noncarrier_median <- median(c_i_temp$Ancestral_Fraction[c_i_temp$Carrier_Status == FALSE])
c_i_carrier_mean <- mean(c_i_temp$Ancestral_Fraction[c_i_temp$Carrier_Status == TRUE])
c_i_noncarrier_mean <- mean(c_i_temp$Ancestral_Fraction[c_i_temp$Carrier_Status == FALSE])

# C_Spec_M
c_m_temp <- joined_admix_5.3_msub %>% select(NPM_Research_ID, Chinese_Comp)
c_m_temp <- cbind.data.frame(c_m_temp, c_spec_m_carriers_status, rep("Chinese_Specific", nrow(c_m_temp)), rep("Malay", nrow(c_m_temp)))
colnames(c_m_temp) <- spec_plot_colnames
c_m_carrier_median <- median(c_m_temp$Ancestral_Fraction[c_m_temp$Carrier_Status == TRUE])
c_m_noncarrier_median <- median(c_m_temp$Ancestral_Fraction[c_m_temp$Carrier_Status == FALSE])
c_m_carrier_mean <- mean(c_m_temp$Ancestral_Fraction[c_m_temp$Carrier_Status == TRUE])
c_m_noncarrier_mean <- mean(c_m_temp$Ancestral_Fraction[c_m_temp$Carrier_Status == FALSE])

# I_Spec_I
i_i_temp <- joined_admix_5.3_isub %>% select(NPM_Research_ID, Indian_Comp)
i_i_temp <- cbind.data.frame(i_i_temp, i_spec_i_carriers_status, rep("Indian_Specific", nrow(i_i_temp)), rep("Indian", nrow(i_i_temp)))
colnames(i_i_temp) <- spec_plot_colnames
i_i_carrier_median <- median(i_i_temp$Ancestral_Fraction[i_i_temp$Carrier_Status == TRUE])
i_i_noncarrier_median <- median(i_i_temp$Ancestral_Fraction[i_i_temp$Carrier_Status == FALSE])
i_i_carrier_mean <- mean(i_i_temp$Ancestral_Fraction[i_i_temp$Carrier_Status == TRUE])
i_i_noncarrier_mean <- mean(i_i_temp$Ancestral_Fraction[i_i_temp$Carrier_Status == FALSE])

# I_Spec_C
i_c_temp <- joined_admix_5.3_csub %>% select(NPM_Research_ID, Indian_Comp)
i_c_temp <- cbind.data.frame(i_c_temp, i_spec_c_carriers_status, rep("Indian_Specific", nrow(i_c_temp)), rep("Chinese", nrow(i_c_temp)))
colnames(i_c_temp) <- spec_plot_colnames
i_c_carrier_median <- median(i_c_temp$Ancestral_Fraction[i_c_temp$Carrier_Status == TRUE])
i_c_noncarrier_median <- median(i_c_temp$Ancestral_Fraction[i_c_temp$Carrier_Status == FALSE])
i_c_carrier_mean <- mean(i_c_temp$Ancestral_Fraction[i_c_temp$Carrier_Status == TRUE])
i_c_noncarrier_mean <- mean(i_c_temp$Ancestral_Fraction[i_c_temp$Carrier_Status == FALSE])

# I_Spec_M
i_m_temp <- joined_admix_5.3_msub %>% select(NPM_Research_ID, Indian_Comp)
i_m_temp <- cbind.data.frame(i_m_temp, i_spec_m_carriers_status, rep("Indian_Specific", nrow(i_m_temp)), rep("Malay", nrow(i_m_temp)))
colnames(i_m_temp) <- spec_plot_colnames
i_m_carrier_median <- median(i_m_temp$Ancestral_Fraction[i_m_temp$Carrier_Status == TRUE])
i_m_noncarrier_median <- median(i_m_temp$Ancestral_Fraction[i_m_temp$Carrier_Status == FALSE])
i_m_carrier_mean <- mean(i_m_temp$Ancestral_Fraction[i_m_temp$Carrier_Status == TRUE])
i_m_noncarrier_mean <- mean(i_m_temp$Ancestral_Fraction[i_m_temp$Carrier_Status == FALSE])

# M_Spec_M
m_m_temp <- joined_admix_5.3_msub %>% select(NPM_Research_ID, Malay_Comp)
m_m_temp <- cbind.data.frame(m_m_temp, m_spec_m_carriers_status, rep("Malay_Specific", nrow(m_m_temp)), rep("Malay", nrow(m_m_temp)))
colnames(m_m_temp) <- spec_plot_colnames
m_m_carrier_median <- median(m_m_temp$Ancestral_Fraction[m_m_temp$Carrier_Status == TRUE])
m_m_noncarrier_median <- median(m_m_temp$Ancestral_Fraction[m_m_temp$Carrier_Status == FALSE])
m_m_carrier_mean <- mean(m_m_temp$Ancestral_Fraction[m_m_temp$Carrier_Status == TRUE])
m_m_noncarrier_mean <- mean(m_m_temp$Ancestral_Fraction[m_m_temp$Carrier_Status == FALSE])

# M_Spec_C
m_c_temp <- joined_admix_5.3_csub %>% select(NPM_Research_ID, Malay_Comp)
m_c_temp <- cbind.data.frame(m_c_temp, m_spec_c_carriers_status, rep("Malay_Specific", nrow(m_c_temp)), rep("Chinese", nrow(m_c_temp)))
colnames(m_c_temp) <- spec_plot_colnames
m_c_carrier_median <- median(m_c_temp$Ancestral_Fraction[m_c_temp$Carrier_Status == TRUE])
m_c_noncarrier_median <- median(m_c_temp$Ancestral_Fraction[m_c_temp$Carrier_Status == FALSE])
m_c_carrier_mean <- mean(m_c_temp$Ancestral_Fraction[m_c_temp$Carrier_Status == TRUE])
m_c_noncarrier_mean <- mean(m_c_temp$Ancestral_Fraction[m_c_temp$Carrier_Status == FALSE])

# M_Spec_I
m_i_temp <- joined_admix_5.3_isub %>% select(NPM_Research_ID, Malay_Comp)
m_i_temp <- cbind.data.frame(m_i_temp, m_spec_i_carriers_status, rep("Malay_Specific", nrow(m_i_temp)), rep("Indian", nrow(m_i_temp)))
colnames(m_i_temp) <- spec_plot_colnames
m_i_carrier_median <- median(m_i_temp$Ancestral_Fraction[m_i_temp$Carrier_Status == TRUE])
m_i_noncarrier_median <- median(m_i_temp$Ancestral_Fraction[m_i_temp$Carrier_Status == FALSE])
m_i_carrier_mean <- mean(m_i_temp$Ancestral_Fraction[m_i_temp$Carrier_Status == TRUE])
m_i_noncarrier_mean <- mean(m_i_temp$Ancestral_Fraction[m_i_temp$Carrier_Status == FALSE])

spec_plot_df <- rbind.data.frame(c_c_temp, c_i_temp, c_m_temp, i_i_temp, i_c_temp, i_m_temp, m_m_temp, m_c_temp, m_i_temp)

# To consolidate number of carriers and range of ancestral fraction, by ancestry, of ancestry-specific variants for each ancestral component (Supplementary Table 5)
c_c_sum_carrier <- c_c_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Chinese") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
c_i_sum_carrier <- c_i_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Chinese") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
c_m_sum_carrier <- c_m_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Chinese") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
i_c_sum_carrier <- i_c_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Indian") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
i_i_sum_carrier <- i_i_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Indian") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
i_m_sum_carrier <- i_m_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Indian") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
m_c_sum_carrier <- m_c_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Malay") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
m_i_sum_carrier <- m_i_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Malay") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
m_m_sum_carrier <- m_m_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% count() %>% mutate(Var_ancestral_comp = "Malay") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
cim_sum_carrier <- rbind.data.frame(c_c_sum_carrier, c_i_sum_carrier, c_m_sum_carrier, i_c_sum_carrier, i_i_sum_carrier, i_m_sum_carrier, m_c_sum_carrier, m_i_sum_carrier, m_m_sum_carrier) %>% rename(Carrier_Count = n)

c_c_sum_ancestral.comp <- c_c_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Chinese") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
c_i_sum_ancestral.comp <- c_i_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Chinese") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
c_m_sum_ancestral.comp <- c_m_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Chinese") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
i_c_sum_ancestral.comp <- i_c_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Indian") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
i_i_sum_ancestral.comp <- i_i_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Indian") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
i_m_sum_ancestral.comp <- i_m_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Indian") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
m_c_sum_ancestral.comp <- m_c_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Malay") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
m_i_sum_ancestral.comp <- m_i_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Malay") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
m_m_sum_ancestral.comp <- m_m_temp %>% group_by(Carrier_Status, Genetic_Ancestry) %>% summarise(Min_ancestral_comp = min(Ancestral_Fraction), Max_ancestral_comp = max(Ancestral_Fraction), Median_ancestral_comp = median(Ancestral_Fraction), Mean_ancestral_comp = mean(Ancestral_Fraction)) %>% mutate(Var_ancestral_comp = "Malay") %>% relocate(Var_ancestral_comp, .before = Carrier_Status)
cim_sum_ancestral.comp <- rbind.data.frame(c_c_sum_ancestral.comp, c_i_sum_ancestral.comp, c_m_sum_ancestral.comp, i_c_sum_ancestral.comp, i_i_sum_ancestral.comp, i_m_sum_ancestral.comp, m_c_sum_ancestral.comp, m_i_sum_ancestral.comp, m_m_sum_ancestral.comp)
#write.table(cim_sum_ancestral.comp, "OUTPUT_admixture_suppt9_ancestralcomp.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# To compute adjusted p-values for spec_plot_df (Supplementary Table 5)
p_c_c <- wilcox.test(joined_admix_5.3_csub$Chinese_Comp~c_spec_c_carriers_status)$p.value
p_c_i <- wilcox.test(joined_admix_5.3_isub$Chinese_Comp~c_spec_i_carriers_status)$p.value
p_c_m <- wilcox.test(joined_admix_5.3_msub$Chinese_Comp~c_spec_m_carriers_status)$p.value
p_i_c <- wilcox.test(joined_admix_5.3_csub$Indian_Comp~i_spec_c_carriers_status)$p.value
p_i_i <- wilcox.test(joined_admix_5.3_isub$Indian_Comp~i_spec_i_carriers_status)$p.value
p_i_m <- wilcox.test(joined_admix_5.3_msub$Indian_Comp~i_spec_m_carriers_status)$p.value
p_m_c <- wilcox.test(joined_admix_5.3_csub$Malay_Comp~m_spec_c_carriers_status)$p.value
p_m_i <- wilcox.test(joined_admix_5.3_isub$Malay_Comp~m_spec_i_carriers_status)$p.value
p_m_m <- wilcox.test(joined_admix_5.3_msub$Malay_Comp~m_spec_m_carriers_status)$p.value
spec_var_wilcox_pvals <- c(p_c_c, p_c_i, p_c_m, p_i_c, p_i_i, p_i_m, p_m_c, p_m_i, p_m_m)
p.adjust(spec_var_wilcox_pvals, method = "BH")

# To generate barplot for ancestral component fractions of ancestry-specific variant carriers in each ancestry group (Figure 2B)
# Rearrange the order of columns and rename variant carrier status
admix <- spec_plot_df
col_order <- c("Specific_Ancestry", "NPM_Research_ID", "Genetic_Ancestry", "Carrier_Status", "Ancestral_Fraction")
admix_a <- admix[, col_order]
admix_b <- admix_a %>% mutate_at("Carrier_Status", str_replace, "FALSE", "Non-Carrier") %>% mutate_at("Carrier_Status", str_replace, "TRUE", "Carrier")

# Generate barplot of Chinese ancestral component fraction for carriers of Chinese-specific variants across ancestry groups: Figure 2B_left panel
c_spec <- admix_b %>% filter(Specific_Ancestry == "Chinese_Specific") %>% arrange(Genetic_Ancestry)
c_spec$Ancestral_Fraction <- as.numeric(c_spec$Ancestral_Fraction)
c_spec$Genetic_Ancestry <- factor(c_spec$Genetic_Ancestry, levels = unique(c_spec$Genetic_Ancestry))
c_spec$Carrier_Status <- factor(c_spec$Carrier_Status, levels = rev(unique(c_spec$Carrier_Status)))
c <- ggplot(c_spec, aes(x = Genetic_Ancestry, y = Ancestral_Fraction, fill = Carrier_Status)) + geom_boxplot() + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, linetype = "solid", size = 0.5), legend.position = "none")
c
# Generate barplot of Indian ancestral component fraction for carriers of Indian-specific variants across ancestry groups: Figure 2B_middle panel
i_spec <- admix_b %>% filter(Specific_Ancestry == "Indian_Specific") %>% arrange(Genetic_Ancestry)
i_spec$Ancestral_Fraction <- as.numeric(i_spec$Ancestral_Fraction)
i_spec$Genetic_Ancestry <- factor(i_spec$Genetic_Ancestry, levels = unique(i_spec$Genetic_Ancestry))
i_spec$Carrier_Status <- factor(i_spec$Carrier_Status, levels = rev(unique(i_spec$Carrier_Status)))
i <- ggplot(i_spec, aes(x = Genetic_Ancestry, y = Ancestral_Fraction, fill = Carrier_Status)) + geom_boxplot() + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, linetype = "solid", size = 0.5), legend.position = "none") + ylim(0,1.00)
i
# Generate barplot of Malay ancestral component fraction for carriers of Malay-specific variants across ancestry groups: Figure 2B_right panel
m_spec <- admix_b %>% filter(Specific_Ancestry == "Malay_Specific") %>% arrange(Genetic_Ancestry)
m_spec$Ancestral_Fraction <- as.numeric(m_spec$Ancestral_Fraction)
m_spec$Genetic_Ancestry <- factor(m_spec$Genetic_Ancestry, levels = unique(m_spec$Genetic_Ancestry))
m_spec$Carrier_Status <- factor(m_spec$Carrier_Status, levels = rev(unique(m_spec$Carrier_Status)))
m <- ggplot(m_spec, aes(x = Genetic_Ancestry, y = Ancestral_Fraction, fill = Carrier_Status)) + geom_boxplot() + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, linetype = "solid", size = 0.5), legend.background = element_blank(), legend.key = element_blank())
m

#To consolidate all three barplots into a single panel (Figure 2B)
c + i + m
#p-value brackets were added to the figure separately using Adobe Illustrator

