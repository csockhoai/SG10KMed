# This code is to compare the correlation of 15X vs 30X WGS samples in PGX genes for potential batch effects as displayed in Supplementary Figure 5
# Due to different variant calling tools used for PGX genes, comparison of 15X vs 30X was done by tool type. There are 4 independent comparisons:
# (1) for genes using CYRIUS caller (CYP2D6)
# (2) for genes using Aldy caller
# (3) for genes using force-call from VCF : except G6PD and HLA
# (4) for genes using force-call from VCF : G6PD and HLA-A, HLA-B
# Please note the code block for each gene is repeated, hence variables may be overwritten.

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)

# ============================================================
# Section 1.0 : For PGX genes using CYRIUS caller - CYP2D6
# ============================================================

# Create individual-level summary by read depth
indv <- read.table("r5.3_9105_UnrelatedSamples_Demographics_20211117.txt", header = TRUE, sep = "\t")
dem_15x <- indv %>% filter(target_depth == "15x")
dem_30x <- indv %>% filter(target_depth == "30x")

# Annotate by phenotype
df1 <- read.csv("R_CYP2D6_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_CYP2D6_SGallele_v3_final.csv", header = TRUE)
tabl_alle_b <- tabl_alle %>% select(Allele_SG, Activity_score)
df2 <- df1 %>% left_join(tabl_alle_b, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle_b, c("Allele_2" = "Allele_SG")) %>% rename(Activity_alle1 = Activity_score.x, Activity_alle2 = Activity_score.y) %>% mutate(phenotype_score = Activity_alle1 + Activity_alle2)
df2b <- df2 %>% mutate(phenotype = case_when(phenotype_score > 2 ~ "Ultrarapid_metabolizer",
                                             (phenotype_score > 1 & phenotype_score <= 2) ~ "Normal_metabolizer",
                                             (phenotype_score >= 0.5 & phenotype_score <= 1) ~ "Intermediate_metabolizer",
                                             (phenotype_score >= 0 & phenotype_score < 0.5) ~ "Poor_metabolizer",
                                             phenotype_score < 0 ~ "Uncallable",
                                             TRUE ~ "Indeterminate"))
df3 <- df2b %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Normal_metabolizer"] <- "CYP2D6_NM"
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "CYP2D6_IM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "CYP2D6_PM"
df3$phenotype[df3$phenotype == "Ultrarapid_metabolizer"] <- "CYP2D6_UM"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "CYP2D6_ID"
cyp2d6_15x <- left_join(dem_15x, df3)
cyp2d6_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cyp2d6_15x_1 <- cyp2d6_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Uncallable") %>% replace(is.na(.),0) %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x)
cyp2d6_15x_1 <- as.data.frame(cyp2d6_15x_1)
cyp2d6_15x_2 <- cyp2d6_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/2457, IND = I/1287, MY = M/885) %>% slice(1:(n() - 1))
cyp2d6_15x_3 <- cyp2d6_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cyp2d6_30x_1 <- cyp2d6_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Uncallable") %>% replace(is.na(.),0) %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x)
cyp2d6_30x_1 <- as.data.frame(cyp2d6_30x_1)
cyp2d6_30x_2 <- cyp2d6_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1413, IND = I/194, MY = M/239) %>% slice(1:(n() - 1))
cyp2d6_30x_3 <- cyp2d6_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cyp2d6_1 <- full_join(cyp2d6_15x_3, cyp2d6_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cyp2d6_1$ancestry <- as.factor(cyp2d6_1$ancestry)
ggplot(cyp2d6_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CYP2D6 CF (15X samples)") + ylab("CYP2D6 CF (30X samples)")

cor.test(cyp2d6_1$s_15x, cyp2d6_1$s_30x)
cor.test(cyp2d6_1$s_15x, cyp2d6_1$s_30x)$p.value

# ==================================================================================================
# Section 2.0 : PGX genes using Aldy Caller - CYP2B6, CYP2C9, CYP2C19, CYP3A4, CYP3A5, NAT2, DPYD*  
# ==================================================================================================

# Create individual-level summary by read depth
indv <- read.table("r5.3_9105_UnrelatedSamples_Demographics_20211117.txt", header = TRUE, sep = "\t")
dem_15x <- indv %>% filter(target_depth == "15x")
dem_30x <- indv %>% filter(target_depth == "30x")

# =======================================================
# Section 2.1 : PGX genes using Aldy Caller - CYP2B6
# =======================================================

# Annotate by phenotype
df1 <- read.csv("R_CYP2B6_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_CYP2B6_v1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_CYP2B6_v1_phenotypetable.csv", header = TRUE)
df2 <- df1 %>% left_join(tabl_alle, c("Allele_1" = "Allele")) %>% left_join(tabl_alle, c("Allele_2" = "Allele")) %>% rename(Func_alle1 = Allele_function.x, Func_alle2 = Allele_function.y)
df2b <- c()
for(i in 1:nrow(df2))
{
  myrow <- df2[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df2b <- rbind(df2b, unlist(c(myrow, phenotype)))
}
colnames(df2b) <- c('npm_research_id', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df2b <- data.frame(df2b)
df3 <- df2b %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "CYP2B6_IM"
df3$phenotype[df3$phenotype == "Normal_metabolizer"] <- "CYP2B6_NM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "CYP2B6_PM"
df3$phenotype[df3$phenotype == "Rapid_metabolizer"] <- "CYP2B6_RM"
df3$phenotype[df3$phenotype == "Ultrarapid_metabolizer"] <- "CYP2B6_UM"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "CYP2B6_ID"
cyp2b6_15x <- left_join(dem_15x, df3)
cyp2b6_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cyp2b6_15x_1 <- cyp2b6_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% replace(is.na(.),0) %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x)
cyp2b6_15x_1 <- as.data.frame(cyp2b6_15x_1)
cyp2b6_15x_2 <- cyp2b6_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3877, IND = I/1702, MY = M/1243) %>% slice(1:(n() - 1))
cyp2b6_15x_3 <- cyp2b6_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cyp2b6_30x_1 <- cyp2b6_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% replace(is.na(.),0) %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x)
cyp2b6_30x_1 <- as.data.frame(cyp2b6_30x_1)
cyp2b6_30x_2 <- cyp2b6_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
cyp2b6_30x_3 <- cyp2b6_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cyp2b6_1 <- full_join(cyp2b6_15x_3, cyp2b6_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cyp2b6_1$ancestry <- as.factor(cyp2b6_1$ancestry)
ggplot(cyp2b6_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CYP2B6 CF (15X samples)") + ylab("CYP2B6 CF (30X samples)")

cor.test(cyp2b6_1$s_15x, cyp2b6_1$s_30x)
cor.test(cyp2b6_1$s_15x, cyp2b6_1$s_30x)$p.value

# =======================================================
# Section 2.2 : PGX genes using Aldy Caller - CYP2C9
# =======================================================

# Annotate by phenotype
df1 <- read.csv("R_CYP2C9_v1.csv", header = TRUE)
df3 <- df1 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Intermediate_metabolizer "] <- "CYP2C9_IM"
df3$phenotype[df3$phenotype == "Normal_metabolizer "] <- "CYP2C9_NM"
df3$phenotype[df3$phenotype == "Poor_metabolizer "] <- "CYP2C9_PM"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "CYP2C9_ID"
cyp2c9_15x <- left_join(dem_15x, df3)
cyp2c9_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cyp2c9_15x_1 <- cyp2c9_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% replace(is.na(.),0) %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x)
cyp2c9_15x_1 <- as.data.frame(cyp2c9_15x_1)
cyp2c9_15x_2 <- cyp2c9_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3877, IND = I/1702, MY = M/1244) %>% slice(1:(n() - 1))
cyp2c9_15x_3 <- cyp2c9_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cyp2c9_30x_1 <- cyp2c9_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
cyp2c9_30x_1 <- as.data.frame(cyp2c9_30x_1)
cyp2c9_30x_2 <- cyp2c9_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
cyp2c9_30x_3 <- cyp2c9_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cyp2c9_1 <- full_join(cyp2c9_15x_3, cyp2c9_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cyp2c9_1$ancestry <- as.factor(cyp2c9_1$ancestry)
ggplot(cyp2c9_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CYP2C9 CF (15X samples)") + ylab("CYP2C9 CF (30X samples)")

cor.test(cyp2c9_1$s_15x, cyp2c9_1$s_30x)
cor.test(cyp2c9_1$s_15x, cyp2c9_1$s_30x)$p.value

# =======================================================
# Section 2.3 : PGX genes using Aldy Caller - CYP2C19
# =======================================================

# Annotate by phenotype
df1 <- read.csv("R_CYP2C19_v1.csv", header = TRUE)
tabl_phe <- read.csv("R_CYP2C19_v1_phenotypetable.csv", header = TRUE)
df2 <- df1 %>% left_join(tabl_phe, c("Diplotype_updated" = "Diplotype_CPIC"))
df3 <- df2 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "CYP2C19_IM"
df3$phenotype[df3$phenotype == "Normal_metabolizer "] <- "CYP2C19_NM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "CYP2C19_PM"
df3$phenotype[df3$phenotype == "Rapid_metabolizer "] <- "CYP2C19_RM"
df3$phenotype[df3$phenotype == "Ultrarapid_metabolizer "] <- "CYP2C19_UM"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "CYP2C19_ID"
cyp2c19_15x <- left_join(dem_15x, df3)
cyp2c19_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cyp2c19_15x_1 <- cyp2c19_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
cyp2c19_15x_1 <- as.data.frame(cyp2c19_15x_1)
cyp2c19_15x_2 <- cyp2c19_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3877, IND = I/1702, MY = M/1244) %>% slice(1:(n() - 1))
cyp2c19_15x_3 <- cyp2c19_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cyp2c19_30x_1 <- cyp2c19_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
cyp2c19_30x_1 <- as.data.frame(cyp2c19_30x_1)
cyp2c19_30x_2 <- cyp2c19_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
cyp2c19_30x_3 <- cyp2c19_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cyp2c19_1 <- full_join(cyp2c19_15x_3, cyp2c19_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cyp2c19_1$ancestry <- as.factor(cyp2c19_1$ancestry)
ggplot(cyp2c19_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CYP2C19 CF (15X samples)") + ylab("CYP2C19 CF (30X samples)")

cor.test(cyp2c19_1$s_15x, cyp2c19_1$s_30x)
cor.test(cyp2c19_1$s_15x, cyp2c19_1$s_30x)$p.value

# =======================================================
# Section 2.4 : PGX genes using Aldy Caller - CYP3A4
# =======================================================

# Annotate by phenotype
df1 <- read.csv("R_CYP3A4_v1.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(Diplotype == "*1/*1" ~ "Normal_metabolism",
                                            (Diplotype == "*1/*22" | Diplotype == "*22/*22") ~ "Decreased_metabolism",
                                            Diplotype == "Untyped" ~ "Untyped",
                                            TRUE ~ "Unknown"))
df3 <- df2 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Decreased_metabolism"] <- "CYP3A4_PM"
df3$phenotype[df3$phenotype == "Normal_metabolism"] <- "CYP3A4_NM"
df3$phenotype[df3$phenotype == "Unknown"] <- "CYP3A4_ID"
cyp3a4_15x <- left_join(dem_15x, df3)
cyp3a4_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cyp3a4_15x_1 <- cyp3a4_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
cyp3a4_15x_1 <- as.data.frame(cyp3a4_15x_1)
cyp3a4_15x_2 <- cyp3a4_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3877, IND = I/1702, MY = M/1244) %>% slice(1:(n() - 1))
cyp3a4_15x_3 <- cyp3a4_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cyp3a4_30x_1 <- cyp3a4_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
cyp3a4_30x_1 <- as.data.frame(cyp3a4_30x_1)
cyp3a4_30x_2 <- cyp3a4_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
cyp3a4_30x_3 <- cyp3a4_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cyp3a4_1 <- full_join(cyp3a4_15x_3, cyp3a4_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cyp3a4_1$ancestry <- as.factor(cyp3a4_1$ancestry)
ggplot(cyp3a4_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CYP3A4 CF (15X samples)") + ylab("CYP3A4 CF (30X samples)")

cor.test(cyp3a4_1$s_15x, cyp3a4_1$s_30x)
cor.test(cyp3a4_1$s_15x, cyp3a4_1$s_30x)$p.value

# =======================================================
# Section 2.5 : PGX genes using Aldy Caller - CYP3A5
# =======================================================

# Annotate by phenotype
df1 <- read.csv("R_CYP3A5_v1.csv", header = TRUE)
df1[c("Allele_1", "Allele_2")][is.na(df1[c("Allele_1", "Allele_2")])] <- 0
tabl_alle <- read.csv("R_CYP3A5_v1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_CYP3A5_v1_phenotypetable.csv", header = TRUE)
df2 <- df1 %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename("Func_alle1" = "Allele_func.x", "Func_alle2" = "Allele_func.y")
df2b <- c()
for(i in 1:nrow(df2))
{
  myrow <- df2[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df2b <- rbind(df2b, unlist(c(myrow, phenotype)))
}
colnames(df2b) <- c('npm_research_id', 'Diplotype', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df2b <- data.frame(df2b)
df3 <- df2b %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "CYP3A5_IM"
df3$phenotype[df3$phenotype == "Normal_metabolizer"] <- "CYP3A5_NM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "CYP3A5_PM"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "CYP3A5_ID"
cyp3a5_15x <- left_join(dem_15x, df3)
cyp3a5_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cyp3a5_15x_1 <- cyp3a5_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
cyp3a5_15x_1 <- as.data.frame(cyp3a5_15x_1)
cyp3a5_15x_2 <- cyp3a5_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3875, IND = I/1694, MY = M/1240) %>% slice(1:(n() - 1))
cyp3a5_15x_3 <- cyp3a5_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cyp3a5_30x_1 <- cyp3a5_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
cyp3a5_30x_1 <- as.data.frame(cyp3a5_30x_1)
cyp3a5_30x_2 <- cyp3a5_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1573, IND = I/223, MY = M/322) %>% slice(1:(n() - 1))
cyp3a5_30x_3 <- cyp3a5_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cyp3a5_1 <- full_join(cyp3a5_15x_3, cyp3a5_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cyp3a5_1$ancestry <- as.factor(cyp3a5_1$ancestry)
ggplot(cyp3a5_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CYP3A5 CF (15X samples)") + ylab("CYP3A5 CF (30X samples)")

cor.test(cyp3a5_1$s_15x, cyp3a5_1$s_30x)
cor.test(cyp3a5_1$s_15x, cyp3a5_1$s_30x)$p.value

# =======================================================
# Section 2.6 : PGX genes using Aldy Caller - NAT2
# =======================================================

# Annotate by phenotype
df1 <- read.csv("R_NAT2_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_NAT2_alleletable.csv", header = TRUE)
df2 <- df1 %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename("Func_alle1" = "Allele_func.x", "Func_alle2" = "Allele_func.y")
df2b <- df2 %>% mutate(phenotype = case_when((Func_alle1 == "rapid" & Func_alle2 == "rapid") ~ "Rapid_acetylator",
                                             (Func_alle1 == "slow" & Func_alle2 == "slow") ~ "Slow_acetylator",
                                             (Func_alle1 == "rapid" & Func_alle2 == "slow") ~ "Intermediate_acetylator",
                                             (Func_alle1 == "slow" & Func_alle2 == "rapid") ~ "Intermediate_acetylator",
                                             (Func_alle1 == "rapid" & Func_alle2 == "unknown") ~ "Indeterminate",
                                             (Func_alle1 == "unknown" & Func_alle2 == "rapid") ~ "Indeterminate",
                                             (Func_alle1 == "slow" & Func_alle2 == "unknown") ~ "Indeterminate",
                                             (Func_alle1 == "unknown" & Func_alle2 == "slow") ~ "Indeterminate",
                                             (Func_alle1 == "unknown" & Func_alle2 == "unknown") ~ "Indeterminate",
                                             TRUE ~ "Untyped"))
df3 <- df2b %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Intermediate_acetylator"] <- "NAT2_Int"
df3$phenotype[df3$phenotype == "Rapid_acetylator"] <- "NAT2_Rapid"
df3$phenotype[df3$phenotype == "Slow_acetylator"] <- "NAT2_Slow"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "NAT2_ID"
nat2_15x <- left_join(dem_15x, df3)
nat2_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
nat2_15x_1 <- nat2_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
nat2_15x_1 <- as.data.frame(nat2_15x_1)
nat2_15x_2 <- nat2_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3877, IND = I/1702, MY = M/1244) %>% slice(1:(n() - 1))
nat2_15x_3 <- nat2_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
nat2_30x_1 <- nat2_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
nat2_30x_1 <- as.data.frame(nat2_30x_1)
nat2_30x_2 <- nat2_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/223, MY = M/351) %>% slice(1:(n() - 1))
nat2_30x_3 <- nat2_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

nat2_1 <- full_join(nat2_15x_3, nat2_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
nat2_1$ancestry <- as.factor(nat2_1$ancestry)
ggplot(nat2_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("NAT2 CF (15X samples)") + ylab("NAT2 CF (30X samples)")

cor.test(nat2_1$s_15x, nat2_1$s_30x)
cor.test(nat2_1$s_15x, nat2_1$s_30x)$p.value

# =======================================================
# Section 2.7 : PGX genes using Aldy Caller - DPYD
# =======================================================

# Annotate by phenotype
df1 <- read.csv("R_DPYD_v1.csv", header = TRUE)
df3 <- df1 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "DPYD_IM"
df3$phenotype[df3$phenotype == "Normal_metabolizer"] <- "DPYD_NM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "DPYD_PM"
dpyd_15x <- left_join(dem_15x, df3)
dpyd_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
dpyd_15x_1 <- dpyd_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
dpyd_15x_1 <- as.data.frame(dpyd_15x_1)
dpyd_15x_2 <- dpyd_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3872, IND = I/1692, MY = M/1240) %>% slice(1:(n() - 1))
dpyd_15x_3 <- dpyd_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
dpyd_30x_1 <- dpyd_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
dpyd_30x_1 <- as.data.frame(dpyd_30x_1)
dpyd_30x_2 <- dpyd_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1573, IND = I/223, MY = M/322) %>% slice(1:(n() - 1))
dpyd_30x_3 <- dpyd_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

dpyd_1 <- full_join(dpyd_15x_3, dpyd_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
dpyd_1$ancestry <- as.factor(dpyd_1$ancestry)
ggplot(dpyd_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("DPYD CF (15X samples)") + ylab("DPYD CF (30X samples)")

cor.test(dpyd_1$s_15x, dpyd_1$s_30x)
cor.test(dpyd_1$s_15x, dpyd_1$s_30x)$p.value

# ======================================================================
# Section 2.8 : Combine all PGX genes using Aldy Caller (15X vs 30X)
# ======================================================================

# Combine all Aldy-called genes carrier frequency into one dataframe (15x vs 30x)
aldy_genes <- rbind.data.frame(cyp2b6_1, cyp2c9_1, cyp2c19_1, cyp3a4_1, cyp3a5_1, nat2_1, dpyd_1)
aldy_genes$ancestry <- as.factor(aldy_genes$ancestry)
#aldy_genes$phenotype2 <- as.factor(aldy_genes$phenotype2)

# Generate scatterplot for 15X vs 30X
ggplot(aldy_genes, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("Carrier freq (15X samples)") + ylab("Carrier freq (30X samples)")

cor.test(aldy_genes$s_15x, aldy_genes$s_30x)
cor.test(aldy_genes$s_15x, aldy_genes$s_30x)$p.value

# ======================================================================================================================================
# Section 3.0 : PGX genes using VCF-Force call - CACNA1S, CFTR, CYP4F2, F5, IFNL3, IFNL4, NUDT15, RYR1, SLCO1B1, TPMT, UGT1A1, VKORC1
# ======================================================================================================================================

# Create individual-level summary by read depth
indv <- read.table("r5.3_9105_UnrelatedSamples_Demographics_20211117.txt", header = TRUE, sep = "\t")
dem_15x <- indv %>% filter(target_depth == "15x")
dem_30x <- indv %>% filter(target_depth == "30x")

# =========================================================
# Section 3.1 : PGX genes using VCF-Force call - CACNA1S
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_CACNA1S_v2.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(df1$rs772226819 == 1 ~ "MHS",
                                            df1$rs772226819 == 0 ~ "WT",
                                            TRUE ~ "Untyped"))
df3 <- df2 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "MHS"] <- "CACNA1S_MHS"
df3$phenotype[df3$phenotype == "WT"] <- "CACNA1S_WT"
cacna1s_15x <- left_join(dem_15x, df3)
cacna1s_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cacna1s_15x_1 <- cacna1s_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
cacna1s_15x_1 <- as.data.frame(cacna1s_15x_1)
cacna1s_15x_2 <- cacna1s_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3802, IND = I/1667, MY = M/1220) %>% slice(1:(n() - 1))
cacna1s_15x_3 <- cacna1s_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cacna1s_30x_1 <- cacna1s_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
cacna1s_30x_1 <- as.data.frame(cacna1s_30x_1)
cacna1s_30x_2 <- cacna1s_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
cacna1s_30x_3 <- cacna1s_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cacna1s_1 <- full_join(cacna1s_15x_3, cacna1s_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cacna1s_1$ancestry <- as.factor(cacna1s_1$ancestry)
ggplot(cacna1s_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CACNA1S CF (15X samples)") + ylab("CACNA1S CF (30X samples)")

cor.test(cacna1s_1$s_15x, cacna1s_1$s_30x)
cor.test(cacna1s_1$s_15x, cacna1s_1$s_30x)$p.value

# =========================================================
# Section 3.2 : PGX genes using VCF-Force call - CFTR
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_CFTR_v2.csv", header = TRUE)
df2 <- c()
for(i in 1:nrow(df1))
{
  myrow <- df1[i,]
  nasum <- length(which(is.na(myrow[,2:5])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:5]) > 0))
  {
    temp_phenotype <- "Favourable"
  }else if(nasum == 0 & (sum(myrow[,2:5]) == 0)){
    temp_phenotype <- "WT"
  }
  df2 <- rbind(df2, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2) <- c('npm_research_id', 'rs115545701', 'rs202179988', 'rs78769542', 'rs74503330', 'NAcount','phenotype')
df2 <- data.frame(df2)
df3 <- df2[, c(1,7)]
df3$phenotype[df3$phenotype == "Favourable"] <- "CFTR_sens"
df3$phenotype[df3$phenotype == "WT"] <- "CFTR_WT"
cftr_15x <- left_join(dem_15x, df3)
cftr_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cftr_15x_1 <- cftr_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
cftr_15x_1 <- as.data.frame(cftr_15x_1)
cftr_15x_2 <- cftr_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3667, IND = I/1574, MY = M/1159) %>% slice(1:(n() - 1))
cftr_15x_3 <- cftr_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cftr_30x_1 <- cftr_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
cftr_30x_1 <- as.data.frame(cftr_30x_1)
cftr_30x_2 <- cftr_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1572, IND = I/222, MY = M/350) %>% slice(1:(n() - 1))
cftr_30x_3 <- cftr_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cftr_1 <- full_join(cftr_15x_3, cftr_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cftr_1$ancestry <- as.factor(cftr_1$ancestry)
ggplot(cftr_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CFTR CF (15X samples)") + ylab("CFTR CF (30X samples)")

cor.test(cftr_1$s_15x, cftr_1$s_30x)
cor.test(cftr_1$s_15x, cftr_1$s_30x)$p.value

# =========================================================
# Section 3.3 : PGX genes using VCF-Force call - CYP4F2
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_CYP4F2_v1.csv", header = TRUE)
df2 <-
  df1 %>% mutate(phenotype = case_when((
    df1$rs2108622 == 0) ~ "WT",
    (df1$rs2108622 == 1 | df1$rs2108622 == 2) ~ "Increased_dose_req",
    TRUE ~ "Untyped"
  )) 
df3 <- df2 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Increased_dose_req"] <- "CYP4F2_sens"
df3$phenotype[df3$phenotype == "WT"] <- "CYP4F2_WT"
cyp4f2_15x <- left_join(dem_15x, df3)
cyp4f2_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
cyp4f2_15x_1 <- cyp4f2_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
cyp4f2_15x_1 <- as.data.frame(cyp4f2_15x_1)
cyp4f2_15x_2 <- cyp4f2_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3714, IND = I/1615, MY = M/1197) %>% slice(1:(n() - 1))
cyp4f2_15x_3 <- cyp4f2_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
cyp4f2_30x_1 <- cyp4f2_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
cyp4f2_30x_1 <- as.data.frame(cyp4f2_30x_1)
cyp4f2_30x_2 <- cyp4f2_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1569, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
cyp4f2_30x_3 <- cyp4f2_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

cyp4f2_1 <- full_join(cyp4f2_15x_3, cyp4f2_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
cyp4f2_1$ancestry <- as.factor(cyp4f2_1$ancestry)
ggplot(cyp4f2_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("CYP4F2 CF (15X samples)") + ylab("CYP4F2 CF (30X samples)")

cor.test(cyp4f2_1$s_15x, cyp4f2_1$s_30x)
cor.test(cyp4f2_1$s_15x, cyp4f2_1$s_30x)$p.value

# =========================================================
# Section 3.4 : PGX genes using VCF-Force call - F5
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_F5_v2.csv", header = TRUE)
df2 <-
  df1 %>% mutate(phenotype = case_when(df1$rs6025 == 1 ~ "Increased_risk",
                                       df1$rs6025 == 0 ~ "WT",
                                       TRUE ~ "Untyped")) 
df3 <- df2 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Increased_risk"] <- "F5_risk"
df3$phenotype[df3$phenotype == "WT"] <- "F5_WT"
f5_15x <- left_join(dem_15x, df3)
f5_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
f5_15x_1 <- f5_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
f5_15x_1 <- as.data.frame(f5_15x_1)
f5_15x_2 <- f5_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3803, IND = I/1648, MY = M/1214) %>% slice(1:(n() - 1))
f5_15x_3 <- f5_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
f5_30x_1 <- f5_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
f5_30x_1 <- as.data.frame(f5_30x_1)
f5_30x_2 <- f5_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
f5_30x_3 <- f5_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

f5_1 <- full_join(f5_15x_3, f5_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
f5_1$ancestry <- as.factor(f5_1$ancestry)
ggplot(f5_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("F5 CF (15X samples)") + ylab("F5 CF (30X samples)")

cor.test(f5_1$s_15x, f5_1$s_30x)
cor.test(f5_1$s_15x, f5_1$s_30x)$p.value

# =========================================================
# Section 3.5 : PGX genes using VCF-Force call - IFNL3
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_IFNL3_v2.csv", header = TRUE)
df2 <- c()
for(i in 1:nrow(df1))
{
  myrow <- df1[i,]
  nasum <- length(which(is.na(myrow[,2:3])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:3]) > 0))
  {
    temp_phenotype <- "Decreased_response"
  } else if(nasum == 0 & (sum(myrow[,2:3]) == 0))
  {
    temp_phenotype <- "WT"
  }
  df2 <- rbind(df2, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2) <- c('npm_research_id', 'rs12979860', 'rs8099917', 'NAcount','phenotype')
df2 <- data.frame(df2)
df3 <- df2[, c(1,5)]
df3$phenotype[df3$phenotype == "Decreased_response"] <- "IFNL3_sens"
df3$phenotype[df3$phenotype == "WT"] <- "IFNL3_WT"
ifnl3_15x <- left_join(dem_15x, df3)
ifnl3_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
ifnl3_15x_1 <- ifnl3_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
ifnl3_15x_1 <- as.data.frame(ifnl3_15x_1)
ifnl3_15x_2 <- ifnl3_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3595, IND = I/1577, MY = M/1170) %>% slice(1:(n() - 1))
ifnl3_15x_3 <- ifnl3_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
ifnl3_30x_1 <- ifnl3_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
ifnl3_30x_1 <- as.data.frame(ifnl3_30x_1)
ifnl3_30x_2 <- ifnl3_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1573, IND = I/224, MY = M/350) %>% slice(1:(n() - 1))
ifnl3_30x_3 <- ifnl3_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

ifnl3_1 <- full_join(ifnl3_15x_3, ifnl3_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
ifnl3_1$ancestry <- as.factor(ifnl3_1$ancestry)
ggplot(ifnl3_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("IFNL3 CF (15X samples)") + ylab("IFNL3 CF (30X samples)")

cor.test(ifnl3_1$s_15x, ifnl3_1$s_30x)
cor.test(ifnl3_1$s_15x, ifnl3_1$s_30x)$p.value

# =========================================================
# Section 3.6 : PGX genes using VCF-Force call - IFNL4
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_IFNL4_v2.csv", header = TRUE)
df2 <- c()
for(i in 1:nrow(df1))
{
  myrow <- df1[i,]
  nasum <- length(which(is.na(myrow[,2:3])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:3]) > 0))
  {
    temp_phenotype <- "Decreased_response"
  } else if(nasum == 0 & (sum(myrow[,2:3]) == 0))
  {
    temp_phenotype <- "WT"
  }
  df2 <- rbind(df2, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2) <- c('npm_research_id', 'rs12979860', 'rs11322783CT.C', 'NAcount','phenotype')
df2 <- data.frame(df2)
df3 <- df2[, c(1,5)]
df3$phenotype[df3$phenotype == "Decreased_response"] <- "IFNL4_sens"
df3$phenotype[df3$phenotype == "WT"] <- "IFNL4_WT"
ifnl4_15x <- left_join(dem_15x, df3)
ifnl4_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
ifnl4_15x_1 <- ifnl4_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
ifnl4_15x_1 <- as.data.frame(ifnl4_15x_1)
ifnl4_15x_2 <- ifnl4_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3512, IND = I/1548, MY = M/1171) %>% slice(1:(n() - 1))
ifnl4_15x_3 <- ifnl4_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
ifnl4_30x_1 <- ifnl4_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
ifnl4_30x_1 <- as.data.frame(ifnl4_30x_1)
ifnl4_30x_2 <- ifnl4_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1567, IND = I/220, MY = M/350) %>% slice(1:(n() - 1))
ifnl4_30x_3 <- ifnl4_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

ifnl4_1 <- full_join(ifnl4_15x_3, ifnl4_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
ifnl4_1$ancestry <- as.factor(ifnl4_1$ancestry)
ggplot(ifnl4_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("IFNL4 CF (15X samples)") + ylab("IFNL4 CF (30X samples)")

cor.test(ifnl4_1$s_15x, ifnl4_1$s_30x)
cor.test(ifnl4_1$s_15x, ifnl4_1$s_30x)$p.value

# =========================================================
# Section 3.7 : PGX genes using VCF-Force call - NUDT15
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_NUDT15_v2.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(df1$rs116855232 == 0 ~ "Normal_metabolizer",
                                            df1$rs116855232 == 1 ~ "Intermediate_metabolizer",
                                            df1$rs116855232 == 2 ~ "Poor_metabolizer",
                                            TRUE ~ "Untyped"))
df3 <- df2 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Normal_metabolizer"] <- "NUDT15_NM"
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "NUDT15_IM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "NUDT15_PM"
nudt15_15x <- left_join(dem_15x, df3)
nudt15_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
nudt15_15x_1 <- nudt15_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
nudt15_15x_1 <- as.data.frame(nudt15_15x_1)
nudt15_15x_2 <- nudt15_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3791, IND = I/1659, MY = M/1209) %>% slice(1:(n() - 1))
nudt15_15x_3 <- nudt15_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
nudt15_30x_1 <- nudt15_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
nudt15_30x_1 <- as.data.frame(nudt15_30x_1)
nudt15_30x_2 <- nudt15_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1571, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
nudt15_30x_3 <- nudt15_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

nudt15_1 <- full_join(nudt15_15x_3, nudt15_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
nudt15_1$ancestry <- as.factor(nudt15_1$ancestry)
ggplot(nudt15_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("NUDT15 CF (15X samples)") + ylab("NUDT15 CF (30X samples)")

cor.test(nudt15_1$s_15x, nudt15_1$s_30x)
cor.test(nudt15_1$s_15x, nudt15_1$s_30x)$p.value

# =========================================================
# Section 3.8 : PGX genes using VCF-Force call - RYR1
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_RYR1_v2.csv", header = TRUE)
df2 <- c()
for(i in 1:nrow(df1))
{
  myrow <- df1[i,]
  nasum <- length(which(is.na(myrow[,2:4])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:4]) > 0))
  {
    temp_phenotype <- "MHS"
  } else if(nasum == 0 & (sum(myrow[,2:4]) == 0))
  {
    temp_phenotype <- "WT"
  }
  df2 <- rbind(df2, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2) <- c('npm_research_id', 'rs112563513', 'rs121918593', 'rs118192168', 'NAcount','phenotype')
df2 <- data.frame(df2)
df3 <- df2[, c(1,6)]
df3$phenotype[df3$phenotype == "MHS"] <- "RYR1_MHS"
df3$phenotype[df3$phenotype == "WT"] <- "RYR1_WT"
ryr1_15x <- left_join(dem_15x, df3)
ryr1_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
ryr1_15x_1 <- ryr1_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
ryr1_15x_1 <- as.data.frame(ryr1_15x_1)
ryr1_15x_2 <- ryr1_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3374, IND = I/1503, MY = M/1135) %>% slice(1:(n() - 1))
ryr1_15x_3 <- ryr1_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
ryr1_30x_1 <- ryr1_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
ryr1_30x_1 <- as.data.frame(ryr1_30x_1)
ryr1_30x_2 <- ryr1_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1573, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
ryr1_30x_3 <- ryr1_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

ryr1_1 <- full_join(ryr1_15x_3, ryr1_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
ryr1_1$ancestry <- as.factor(ryr1_1$ancestry)
ggplot(ryr1_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("RYR1 CF (15X samples)") + ylab("RYR1 CF (30X samples)")

cor.test(ryr1_1$s_15x, ryr1_1$s_30x)
cor.test(ryr1_1$s_15x, ryr1_1$s_30x)$p.value

# =========================================================
# Section 3.9 : PGX genes using VCF-Force call - TPMT
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_TPMT_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_TPMT_v1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_TPMT_v1_phenotypetable.csv", header = TRUE)
df2 <- df1 %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename(Func_alle1 = Allele_Func.x, Func_alle2 = Allele_Func.y)
df2b <- c()
for(i in 1:nrow(df2))
{
  myrow <- df2[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df2b <- rbind(df2b, unlist(c(myrow, phenotype)))
}
colnames(df2b) <- c('npm_research_id', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df2b <- data.frame(df2b)
df3 <- df2b %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Normal_metabolizer"] <- "TPMT_NM"
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "TPMT_IM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "TPMT_PM"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "TPMT_ID"
tpmt_15x <- left_join(dem_15x, df3)
tpmt_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
tpmt_15x_1 <- tpmt_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
tpmt_15x_1 <- as.data.frame(tpmt_15x_1)
tpmt_15x_2 <- tpmt_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/2732, IND = I/1024, MY = M/859) %>% slice(1:(n() - 1))
tpmt_15x_3 <- tpmt_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
tpmt_30x_1 <- tpmt_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
tpmt_30x_1 <- as.data.frame(tpmt_30x_1)
tpmt_30x_2 <- tpmt_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1565, IND = I/221, MY = M/347) %>% slice(1:(n() - 1))
tpmt_30x_3 <- tpmt_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

tpmt_1 <- full_join(tpmt_15x_3, tpmt_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
tpmt_1$ancestry <- as.factor(tpmt_1$ancestry)
ggplot(tpmt_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("TPMT CF (15X samples)") + ylab("TPMT CF (30X samples)")

cor.test(tpmt_1$s_15x, tpmt_1$s_30x)
cor.test(tpmt_1$s_15x, tpmt_1$s_30x)$p.value

# =========================================================
# Section 3.10 : PGX genes using VCF-Force call - SLCO1B1
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_SLCO1B1_v2.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(df1$rs4149056 == 0 ~ "WT",
                                            df1$rs4149056 == 1 ~ "Intermediate_function",
                                            df1$rs4149056 == 2 ~ "Low_function", 
                                            TRUE ~"Untyped"))
df3 <- df2 %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Intermediate_function"] <- "SLCO1B1_Int"
df3$phenotype[df3$phenotype == "Low_function"] <- "SLCO1B1_Low"
df3$phenotype[df3$phenotype == "WT"] <- "SLCO1B1_Wt"
slco1b1_15x <- left_join(dem_15x, df3)
slco1b1_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
slco1b1_15x_1 <- slco1b1_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
slco1b1_15x_1 <- as.data.frame(slco1b1_15x_1)
slco1b1_15x_2 <- slco1b1_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3794, IND = I/1672, MY = M/1224) %>% slice(1:(n() - 1))
slco1b1_15x_3 <- slco1b1_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
slco1b1_30x_1 <- slco1b1_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
slco1b1_30x_1 <- as.data.frame(slco1b1_30x_1)
slco1b1_30x_2 <- slco1b1_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1572, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
slco1b1_30x_3 <- slco1b1_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

slco1b1_1 <- full_join(slco1b1_15x_3, slco1b1_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
slco1b1_1$ancestry <- as.factor(slco1b1_1$ancestry)
ggplot(slco1b1_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("SLCO1B1 CF (15X samples)") + ylab("SLCO1B1 CF (30X samples)")

cor.test(slco1b1_1$s_15x, slco1b1_1$s_30x)
cor.test(slco1b1_1$s_15x, slco1b1_1$s_30x)$p.value

# =========================================================
# Section 3.11 : PGX genes using VCF-Force call - UGT1A1
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_UGT1A1_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_UGT1A1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_UGT1A1_phenotypetable.csv", header = TRUE)
df2 <- df1 %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename(Func_alle1 = Allele_Func.x, Func_alle2 = Allele_Func.y)
df2b <- c()
for(i in 1:nrow(df2))
{
  myrow <- df2[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df2b <- rbind(df2b, unlist(c(myrow, phenotype)))
}
colnames(df2b) <- c('npm_research_id', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df2b <- data.frame(df2b)
df3 <- df2b %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Normal_metabolizer"] <- "UGT1A1_NM"
df3$phenotype[df3$phenotype == "Intermediate_metabolizer"] <- "UGT1A1_IM"
df3$phenotype[df3$phenotype == "Poor_metabolizer"] <- "UGT1A1_PM"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "UGT1A1_ID"
ugt1a1_15x <- left_join(dem_15x, df3)
ugt1a1_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
ugt1a1_15x_1 <- ugt1a1_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
ugt1a1_15x_1 <- as.data.frame(ugt1a1_15x_1)
ugt1a1_15x_2 <- ugt1a1_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/2901, IND = I/1280, MY = M/875) %>% slice(1:(n() - 1))
ugt1a1_15x_3 <- ugt1a1_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
ugt1a1_30x_1 <- ugt1a1_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
ugt1a1_30x_1 <- as.data.frame(ugt1a1_30x_1)
ugt1a1_30x_2 <- ugt1a1_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1440, IND = I/212, MY = M/326) %>% slice(1:(n() - 1))
ugt1a1_30x_3 <- ugt1a1_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

ugt1a1_1 <- full_join(ugt1a1_15x_3, ugt1a1_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
ugt1a1_1$ancestry <- as.factor(ugt1a1_1$ancestry)
ggplot(ugt1a1_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("UGT1A1 CF (15X samples)") + ylab("UGT1A1 CF (30X samples)")

cor.test(ugt1a1_1$s_15x, ugt1a1_1$s_30x)
cor.test(ugt1a1_1$s_15x, ugt1a1_1$s_30x)$p.value

# =========================================================
# Section 3.11 : PGX genes using VCF-Force call - VKORC1
# =========================================================

# Annotate by phenotype
df1 <- read.csv("R_VKORC1_v2.csv", header = TRUE)
df2 <- c()
for(i in 1:nrow(df1))
{
  myrow <- df1[i,]
  nasum <- length(which(is.na(myrow[,2:5])))
  df2 <- rbind(df2, unlist(c(myrow, nasum)))
}
colnames(df2) <- c('npm_research_id', 'rs7294', 'rs9934438', 'rs9923231', 'rs2359612', 'NAcount')
df2 <- data.frame(df2)
df2b <- df2 %>% mutate(rs7294_sensitive = case_when(df2$rs7294 == 0 ~ 1, TRUE ~ 0), 
                       rs9934438_sensitive = case_when(df2$rs9934438 != 0 ~ 1, TRUE ~ 0),
                       rs9923231_sensitive = case_when(df2$rs9923231 != 0 ~ 1, TRUE ~ 0),
                       rs2359612_sensitive = case_when(df2$rs2359612 !=2 ~ 1, TRUE ~ 0)) %>% mutate(sensitive_score = rowSums(.[7:10])) 

df2c <- df2b %>% mutate(phenotype = case_when((NAcount == 0 & sensitive_score > 0) ~ "Sensitive", (NAcount == 0 & sensitive_score == 0) ~ "WT", TRUE ~ "Untyped"))
df3 <- df2c %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Sensitive"] <- "VKORC1_sens"
df3$phenotype[df3$phenotype == "WT"] <- "VKORC1_wt"
vkorc1_15x <- left_join(dem_15x, df3)
vkorc1_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
vkorc1_15x_1 <- vkorc1_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
vkorc1_15x_1 <- as.data.frame(vkorc1_15x_1)
vkorc1_15x_2 <- vkorc1_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3062, IND = I/1284, MY = M/1000) %>% slice(1:(n() - 1))
vkorc1_15x_3 <- vkorc1_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
vkorc1_30x_1 <- vkorc1_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
vkorc1_30x_1 <- as.data.frame(vkorc1_30x_1)
vkorc1_30x_2 <- vkorc1_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1568, IND = I/223, MY = M/348) %>% slice(1:(n() - 1))
vkorc1_30x_3 <- vkorc1_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

vkorc1_1 <- full_join(vkorc1_15x_3, vkorc1_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
vkorc1_1$ancestry <- as.factor(vkorc1_1$ancestry)
ggplot(vkorc1_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("VKORC1 CF (15X samples)") + ylab("VKORC1 CF (30X samples)")

cor.test(vkorc1_1$s_15x, vkorc1_1$s_30x)
cor.test(vkorc1_1$s_15x, vkorc1_1$s_30x)$p.value

# =========================================================================
# Section 3.12.1 : Combine all PGX genes using VCF-Force call (15X vs 30X)
# =========================================================================

# Combine all forcecall genes
forcecall_genes <- rbind.data.frame(cacna1s_1, cftr_1, cyp4f2_1, f5_1, ifnl3_1, ifnl4_1, nudt15_1, ryr1_1, slco1b1_1, tpmt_1, ugt1a1_1, vkorc1_1)
forcecall_genes$ancestry <- as.factor(forcecall_genes$ancestry)

# Generate scatterplot
ggplot(forcecall_genes, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("Carrier freq (15X samples)") + ylab("Carrier freq (30X samples)")

cor.test(forcecall_genes$s_15x, forcecall_genes$s_30x)
cor.test(forcecall_genes$s_15x, forcecall_genes$s_30x)$p.value

# ==============================================================================================================
# Section 3.12.1 : Combine all PGX genes using VCF-Force call (15X vs 30X) - excluding genes with low counts
# ==============================================================================================================

# Exclude genes with very low counts of actionable carriers: CACNA1S, CFTR, F5, RYR1
forcecall_genes.2 <- rbind.data.frame(cyp4f2_1, ifnl3_1, ifnl4_1, nudt15_1, slco1b1_1, tpmt_1, ugt1a1_1, vkorc1_1)
forcecall_genes.2$ancestry <- as.factor(forcecall_genes.2$ancestry)

# Generate scatterplot
ggplot(forcecall_genes.2, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("Carrier freq (15X samples)") + ylab("Carrier freq (30X samples)")

cor.test(forcecall_genes.2$s_15x, forcecall_genes.2$s_30x)
cor.test(forcecall_genes.2$s_15x, forcecall_genes.2$s_30x)$p.value

# =====================================================================
# Section 4.0 : PGX genes using VCF-force call : G6PD, HLA-A, HLA-B
# =====================================================================

# ==========================================================
# Section 4.1 : PGX genes using VCF-force call : G6PD
# ==========================================================

# Annotate by phenotype
df1 <- read.csv("R_G6PD_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_G6PD_v1_alleletable.csv", header = TRUE)
tabl_alle_b <- tabl_alle %>% select(Allele_SG, WHO_class)
tabl_phe <- read.csv("R_G6PD_v1_phenotypetable.csv", header = TRUE)
tabl_phe$Allele_2[is.na(tabl_phe$Allele_2)] <- 0
df2 <- df1 %>% left_join(tabl_alle_b, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle_b, c("Allele_2" = "Allele_SG")) %>% rename("WHO_alle1" = "WHO_class.x", "WHO_alle2" = "WHO_class.y")
df2b <- c()
for(i in 1:nrow(df2))
{
  myrow <- df2[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$WHO_alle1 & tabl_phe$Allele_2 == myrow$WHO_alle2),3]  
  df2b <- rbind(df2b, unlist(c(myrow, phenotype)))
}
colnames(df2b) <- c('npm_research_id', 'Diplotype', 'Allele_1', 'Allele_2', 'WHO_alle1', 'WHO_alle2', 'phenotype')
df2b <- data.frame(df2b)
df3 <- df2b %>% select(npm_research_id, phenotype)
df3$phenotype[df3$phenotype == "Normal"] <- "G6PD_Norm"
df3$phenotype[df3$phenotype == "Deficient"] <- "G6PD_Def"
df3$phenotype[df3$phenotype == "Variable"] <- "G6PD_Var"
df3$phenotype[df3$phenotype == "Indeterminate"] <- "G6PD_ID"
g6pd_15x <- left_join(dem_15x, df3)
g6pd_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
g6pd_15x_1 <- g6pd_15x %>% group_by(phenotype, genetic_sex, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
g6pd_15x_1 <- as.data.frame(g6pd_15x_1)
g6pd_15x_2.F <- g6pd_15x_1 %>% filter(genetic_sex == "F") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/2373, IND = I/1009, MY = M/747) %>% slice(1:(n() - 1))
g6pd_15x_3.F <- g6pd_15x_2.F %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
g6pd_15x_2.M <- g6pd_15x_1 %>% filter(genetic_sex == "M") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1554, IND = I/708, MY = M/510) %>% slice(1:(n() - 1))
g6pd_15x_3.M <- g6pd_15x_2.M %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
g6pd_30x_1 <- g6pd_30x %>% group_by(phenotype, genetic_sex, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
g6pd_30x_1 <- as.data.frame(g6pd_30x_1)
g6pd_30x_2.F <- g6pd_30x_1 %>% filter(genetic_sex == "F") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/788, IND = I/97, MY = M/168) %>% slice(1:(n() - 1))
g6pd_30x_3.F <- g6pd_30x_2.F %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
g6pd_30x_2.M <- g6pd_30x_1 %>% filter(genetic_sex == "M") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/787, IND = I/127, MY = M/183) %>% slice(1:(n() - 1))
g6pd_30x_3.M <- g6pd_30x_2.M %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

g6pd_15x_4 <- full_join(g6pd_15x_3.F, g6pd_15x_3.M) %>% replace(is.na(.),0)
g6pd_30x_4 <- full_join(g6pd_30x_3.F, g6pd_30x_3.M) %>% replace(is.na(.),0)
g6pd_1 <- full_join(g6pd_15x_4, g6pd_30x_4) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots - overall
g6pd_1$ancestry <- as.factor(g6pd_1$ancestry)
ggplot(g6pd_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("G6PD CF (15X samples)") + ylab("G6PD CF (30X samples)")

cor.test(g6pd_1$s_15x, g6pd_1$s_30x)
cor.test(g6pd_1$s_15x, g6pd_1$s_30x)$p.value

# Generate scatterplot with ancestry indicated dots - by sex
# data subset - females
g6pd_2.F <- full_join(g6pd_15x_3.F, g6pd_30x_3.F) %>% replace(is.na(.),0)
g6pd_3.F <- g6pd_2.F %>% filter(phenotype != "G6PD_Norm")
g6pd_3.F$ancestry <- as.factor(g6pd_3.F$ancestry)
# data subset - males
g6pd_2.M <- full_join(g6pd_15x_3.M, g6pd_30x_3.M) %>% replace(is.na(.),0)
g6pd_3.M <- g6pd_2.M %>% filter(phenotype != "G6PD_Norm")
g6pd_3.M$ancestry <- as.factor(g6pd_3.M$ancestry)

# scatterplot - females
ggplot(g6pd_3.F, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("G6PD CF (15X samples, F)") + ylab("G6PD CF (30X samples, F)")

cor.test(g6pd_3.F$s_15x, g6pd_3.F$s_30x)

# scatterplot - males
ggplot(g6pd_3.M, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("G6PD CF (15X samples, M)") + ylab("G6PD CF (30X samples, M)")

cor.test(g6pd_3.M$s_15x, g6pd_3.M$s_30x)

# Generate scatterplot with ancestry indicated dots - excluding Normal phenotype carriers, by sex

g6pd_15x_b1 <- g6pd_15x %>% group_by(phenotype, genetic_sex, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype %in% c("G6PD_Def", "G6PD_Var", "G6PD_ID")) %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
g6pd_15x_b1 <- as.data.frame(g6pd_15x_b1)
g6pd_15x_b2.F <- g6pd_15x_b1 %>% filter(genetic_sex == "F") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/131, IND = I/72, MY = M/36) %>% slice(1:(n() - 1))
g6pd_15x_b3.F <- g6pd_15x_b2.F %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
g6pd_15x_b2.M <- g6pd_15x_b1 %>% filter(genetic_sex == "M") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/23, IND = I/25, MY = M/6) %>% slice(1:(n() - 1))
g6pd_15x_b3.M <- g6pd_15x_b2.M %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

g6pd_30x_b1 <- g6pd_30x %>% group_by(phenotype, genetic_sex, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype %in% c("G6PD_Def", "G6PD_Var", "G6PD_ID")) %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
g6pd_30x_b1 <- as.data.frame(g6pd_30x_b1)
g6pd_30x_b2.F <- g6pd_30x_b1 %>% filter(genetic_sex == "F") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/50, IND = I/8, MY = M/15) %>% slice(1:(n() - 1))
g6pd_30x_b3.F <- g6pd_30x_b2.F %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
g6pd_30x_b2.M <- g6pd_30x_b1 %>% filter(genetic_sex == "M") %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/30, IND = I/8, MY = M/12) %>% slice(1:(n() - 1))
g6pd_30x_b3.M <- g6pd_30x_b2.M %>% select(phenotype, genetic_sex, CH, IND, MY) %>% pivot_longer(cols = 3:5, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

g6pd_b.F <- full_join(g6pd_15x_b3.F, g6pd_30x_b3.F) %>% replace(is.na(.),0)
g6pd_b.F$ancestry <- as.factor(g6pd_b.F$ancestry)
g6pd_b.M <- full_join(g6pd_15x_b3.M, g6pd_30x_b3.M) %>% replace(is.na(.),0)
g6pd_b.M$ancestry <- as.factor(g6pd_b.M$ancestry)

# scatterplot - excl Normal, females
ggplot(g6pd_b.F, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("G6PD CF (15X samples, F)") + ylab("G6PD CF (30X samples, F)")

cor.test(g6pd_b.F$s_15x, g6pd_b.F$s_30x)

# scatterplot - excl Normal, males
ggplot(g6pd_b.M, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("G6PD CF (15X samples, M)") + ylab("G6PD CF (30X samples, M)")

cor.test(g6pd_b.M$s_15x, g6pd_b.M$s_30x)

# ==========================================================
# Section 4.2 : PGX genes using VCF-force call : HLA-A
# ==========================================================

# Annotate by phenotype
df1 <- read.csv("R_HLA-A_v1.csv", header = TRUE)
strip_HLA <- function(x){
  paste(unlist(strsplit(x, ":"))[1], unlist(strsplit(x, ":"))[2], sep = ":")
}
newdf1 <- df1
for(i in 1:nrow(newdf1))
{
  newdf1$Allele_1[i] <- strip_HLA(newdf1$Allele_1[i])
  newdf1$Allele_2[i] <- strip_HLA(newdf1$Allele_2[i])
}
df2 <-
  newdf1 %>% mutate(
    Def_alle1_3101 = case_when(Allele_1 == "A*31:01" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_3101 = case_when(Allele_2 == "A*31:01" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle_3101_score = Def_alle1_3101 + Def_alle2_3101,
    phenotype = case_when((Def_alle_3101_score == 1 | Def_alle_3101_score == 2) ~ "Increased_SCAR",
                          Def_alle_3101_score == 0 ~ "WT",
                          TRUE ~ "Untyped")
  )
df3 <- df2 %>% select(npm_research_id, phenotype) 
df3$phenotype[df3$phenotype == "Increased_SCAR"] <- "HLA-A_SCAR"
df3$phenotype[df3$phenotype == "WT"] <- "HLA-A_WT"
hlaa_15x <- left_join(dem_15x, df3)
hlaa_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
hlaa_15x_1 <- hlaa_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
hlaa_15x_1 <- as.data.frame(hlaa_15x_1)
hlaa_15x_2 <- hlaa_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3893, IND = I/1695, MY = M/1245) %>% slice(1:(n() - 1))
hlaa_15x_3 <- hlaa_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
hlaa_30x_1 <- hlaa_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
hlaa_30x_1 <- as.data.frame(hlaa_30x_1)
hlaa_30x_2 <- hlaa_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
hlaa_30x_3 <- hlaa_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

hlaa_1 <- full_join(hlaa_15x_3, hlaa_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
hlaa_1$ancestry <- as.factor(hlaa_1$ancestry)
ggplot(hlaa_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("HLA-A CF (15X samples)") + ylab("HLA-A CF (30X samples)")

cor.test(hlaa_1$s_15x, hlaa_1$s_30x)
cor.test(hlaa_1$s_15x, hlaa_1$s_30x)$p.value

# ==========================================================
# Section 4.3 : PGX genes using VCF-force call : HLA-B
# ==========================================================

# Annotate by phenotype
df1 <- read.csv("R_HLA-B_v1.csv", header = TRUE)
strip_HLB <- function(x){
  paste(unlist(strsplit(x, ":"))[1], unlist(strsplit(x, ":"))[2], sep = ":")
}
newdf1 <- df1
for(i in 1:nrow(newdf1))
{
  newdf1$Allele_1[i] <- strip_HLB(newdf1$Allele_1[i])
  newdf1$Allele_2[i] <- strip_HLB(newdf1$Allele_2[i])
}
df2 <-
  newdf1 %>% mutate(
    Def_alle1_1502 = case_when(Allele_1 == "B*15:02" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_1502 = case_when(Allele_2 == "B*15:02" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle1_5801 = case_when(Allele_1 == "B*58:01" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_5801 = case_when(Allele_2 == "B*58:01" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle1_5701 = case_when(Allele_1 == "B*57:01" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_5701 = case_when(Allele_2 == "B*57:01" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_Balle_score = rowSums(across(where(is.numeric))),
    phenotype = case_when(Def_Balle_score > 0 ~ "Increased_SCAR",
                          Def_Balle_score == 0 ~ "WT",
                          Def_Balle_score < 0 ~ "Untyped")
  )
df3 <- df2 %>% select(npm_research_id, phenotype) 
df3$phenotype[df3$phenotype == "Increased_SCAR"] <- "HLA-A_SCAR"
df3$phenotype[df3$phenotype == "WT"] <- "HLA-A_WT"
hlab_15x <- left_join(dem_15x, df3)
hlab_30x <- left_join(dem_30x, df3)

# Annotate carrier frequency by ancestry for 15X and 30X samples
hlab_15x_1 <- hlab_15x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_15x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_15x) %>% replace(is.na(.),0)
hlab_15x_1 <- as.data.frame(hlab_15x_1)
hlab_15x_2 <- hlab_15x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/3893, IND = I/1694, MY = M/1245) %>% slice(1:(n() - 1))
hlab_15x_3 <- hlab_15x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_15x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)
hlab_30x_1 <- hlab_30x %>% group_by(phenotype, genetic_ethnicity) %>% count(phenotype) %>% rename(s_30x = n) %>% filter(phenotype != "Untyped") %>% pivot_wider(names_from = genetic_ethnicity, values_from = s_30x) %>% replace(is.na(.),0)
hlab_30x_1 <- as.data.frame(hlab_30x_1)
hlab_30x_2 <- hlab_30x_1 %>% bind_rows(summarise(., across(where(is.numeric), sum), across(where(is.character), ~"Total"))) %>% mutate(CH = C/1575, IND = I/224, MY = M/351) %>% slice(1:(n() - 1))
hlab_30x_3 <- hlab_30x_2 %>% select(phenotype, CH, IND, MY) %>% pivot_longer(cols = 2:4, names_to = "ancestry", values_to = "s_30x") %>% unite("phenotype2", phenotype:ancestry, sep = "_", remove = FALSE)

hlab_1 <- full_join(hlab_15x_3, hlab_30x_3) %>% replace(is.na(.),0)

# Generate scatterplot with ancestry indicated dots
hlab_1$ancestry <- as.factor(hlab_1$ancestry)
ggplot(hlab_1, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("HLA-B CF (15X samples)") + ylab("HLA-B CF (30X samples)")

cor.test(hlab_1$s_15x, hlab_1$s_30x)
cor.test(hlab_1$s_15x, hlab_1$s_30x)$p.value

# ===============================================================
# Section 4.4 : PGX genes using VCF-force call : HLA-A & HLA-B
# ===============================================================

# Generate scatterplot of both HLA genes together
hlagenes <- rbind.data.frame(hlaa_1, hlab_1)
hlagenes$ancestry <- as.factor(hlagenes$ancestry)
ggplot(hlagenes, aes(x = s_15x, y = s_30x)) + geom_point(aes(color = ancestry), size = 2) + scale_color_manual(values = c('#d73027', '#4575b4', '#fdae61')) + geom_smooth(method = lm, color = "#000000") + theme(
  panel.background = element_blank(),
  panel.border = element_rect(
    fill = NA,
    linetype = "solid",
    size = 0.5
  ),
  legend.background = element_blank(),
  legend.key = element_blank()
)  + xlab("HLA genes CF (15X samples)") + ylab("HLA genes CF (30X samples)")

cor.test(hlagenes$s_15x, hlagenes$s_30x)
cor.test(hlagenes$s_15x, hlagenes$s_30x)$p.value

