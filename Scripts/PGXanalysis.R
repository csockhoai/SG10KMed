# This code is to analyze the fraction of individuals with a risk allele in one of the 23 pharmacogenes
# and the fractions of individuals with concurrent P/LP variants in a CDC Tier1 genetic condition and a PGX risk allele relevant to their condition (Figure 3).
# This script is independent of PGXfreq.R
# This script has three sections:
# 1) Analysis of individuals with risk alleles in pharmacogenes
# 2) Consolidating carrier frequencies by pharmacophenotypes with therapeutic recommendations (Supplementary Data 7)
# 3) Intersection of individuals with CDC Tier1 PLP variants and relevant PGX alleles
# 4) Descriptive analysis of novel PGX variants

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)
library(janitor)
library(data.table)
library(patchwork)

# ================================
# Pre-02 : Load main data files
# ================================
# individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
dem_sub <- dem %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
# ACMG SF version 3 (v.3) genes list
acmg73 <- read.csv("Z0_ACMG73list.csv", header = TRUE)
# exclusion list of APOB loss-of-function variants: APOB LOF variants are not associated with autosomal dominant familial hypercholesterolemia
apoblof <- read.table("ExclusionList_01_APOB_LOFv.txt", header = TRUE, sep = "\t")

# ===========================================================================
# Section 1.0 : Analysis of individuals with risk alleles in pharmacogenes
# ===========================================================================
# First, to generate a master list of consolidated PGX risk allele carriers by pharmacogene
# For brevity, variables for each pharmacogene in code block is represented by alphabets (e.g. df1A, df2A is coded for CACNA1S; df1B, df2B is coded for CFTR)
# The 23 pharmacogenes are listed in alphabetical order (i.e. A=CACNA1S, B=CFTR, ... V=VKORC1). Refer to Supplementary Data 10 for list of 23 pharmacogenes

df1A <- read.csv("R_CACNA1S_v2.csv", header = TRUE)
df2A <- df1A %>% mutate(phenotype = case_when(df1A$rs772226819 == 1 ~ "MHS",
                                              df1A$rs772226819 == 0 ~ "WT",
                                              TRUE ~ "Untyped"))
df3A <- df2A %>% mutate(CACNA1S = case_when(df2A$phenotype == "MHS" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, CACNA1S)
df4A <- left_join(dem_sub, df2A, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5A_all <- df4A %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CACNA1S = n)

df1B <- read.csv("R_CFTR_v2.csv", header = TRUE)
df2B <- c()
for(i in 1:nrow(df1B))
{
  myrow <- df1B[i,]
  nasum <- length(which(is.na(myrow[,2:5])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:5]) > 0))
  {
    temp_phenotype <- "Favourable"
  }else if(nasum == 0 & (sum(myrow[,2:5]) == 0)){
    temp_phenotype <- "WT"
  }
  df2B <- rbind(df2B, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2B) <- c('npm_research_id', 'rs115545701', 'rs202179988', 'rs78769542', 'rs74503330', 'NAcount','phenotype')
df2B <- data.frame(df2B)
df3B <- df2B %>% mutate(CFTR = case_when(df2B$phenotype == "Favourable" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, CFTR)
df4B <- left_join(dem_sub, df2B, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5B_all <- df4B %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CFTR = n)

df1C <- read.csv("R_CYP2B6_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_CYP2B6_v1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_CYP2B6_v1_phenotypetable.csv", header = TRUE)
df2C <- df1C %>% left_join(tabl_alle, c("Allele_1" = "Allele")) %>% left_join(tabl_alle, c("Allele_2" = "Allele")) %>% rename(Func_alle1 = Allele_function.x, Func_alle2 = Allele_function.y)
df3C <- c()
for(i in 1:nrow(df2C))
{
  myrow <- df2C[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df3C <- rbind(df3C, unlist(c(myrow, phenotype)))
}
colnames(df3C) <- c('npm_research_id', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df3C <- data.frame(df3C)
df4C <-
  df3C %>% mutate(CYP2B6 = case_when((
    phenotype == "Intermediate_metabolizer" |
      phenotype == "Poor_metabolizer"
  ) ~ 1,
  TRUE ~ 0
  )) %>% select(npm_research_id, CYP2B6)
df5C <- left_join(dem_sub, df3C, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df6C_all <- df5C %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CYP2B6 = n)

df1D <- read.csv("R_CYP2C9_v1.csv", header = TRUE)
df1D$phenotype <- trimws(df1D$phenotype, which = c("right"))
df2D <- df1D %>% mutate(CYP2C9 = case_when((phenotype == "Intermediate_metabolizer" | phenotype == "Poor_metabolizer") ~ 1, TRUE ~ 0))
df3D <- df2D %>% select(npm_research_id, CYP2C9)
df4D <- left_join(dem_sub, df1D, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5D_all <- df4D %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CYP2C9 = n)

df1E <- read.csv("R_CYP2C19_v1.csv", header = TRUE)
tabl_phe <- read.csv("R_CYP2C19_v1_phenotypetable.csv", header = TRUE)
df2E <- df1E %>% left_join(tabl_phe, c("Diplotype_updated" = "Diplotype_CPIC"))
df3E <-
  df2E %>% mutate(CYP2C19 = case_when((
    phenotype == "Ultrarapid_metabolizer " |
      phenotype == "Rapid_metabolizer " |
      phenotype == "Intermediate_metabolizer" |
      phenotype == "Poor_metabolizer"
  ) ~ 1,
  TRUE ~ 0
  )) %>% select(npm_research_id, CYP2C19)
df4E <- left_join(dem_sub, df2E, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5E_all <- df4E %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CYP2C19 = n)

df1F <- read.csv("R_CYP2D6_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_CYP2D6_SGallele_v3_final.csv", header = TRUE)
tabl_alle_b <- tabl_alle %>% select(Allele_SG, Activity_score)
df2F <- df1F %>% left_join(tabl_alle_b, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle_b, c("Allele_2" = "Allele_SG")) %>% rename(Activity_alle1 = Activity_score.x, Activity_alle2 = Activity_score.y) %>% mutate(phenotype_score = Activity_alle1 + Activity_alle2)
df3F <- df2F %>% mutate(phenotype = case_when(phenotype_score > 2 ~ "Ultrarapid_metabolizer",
                                              (phenotype_score > 1 & phenotype_score <= 2) ~ "Normal_metabolizer",
                                              (phenotype_score >= 0.5 & phenotype_score <= 1) ~ "Intermediate_metabolizer",
                                              (phenotype_score >= 0 & phenotype_score < 0.5) ~ "Poor_metabolizer",
                                              phenotype_score < 0 ~ "Uncallable",
                                              TRUE ~ "Indeterminate"))
df4F <-
  df3F %>% mutate(CYP2D6 = case_when((
    phenotype == "Ultrarapid_metabolizer" |
      phenotype == "Intermediate_metabolizer" |
      phenotype == "Poor_metabolizer"
  ) ~ 1,
  TRUE ~ 0
  )) %>% select(npm_research_id, CYP2D6)
df5F <- left_join(dem_sub, df3F, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df6F_all <- df5F %>% filter(phenotype != "Uncallable") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CYP2D6 = n)

df1G <- read.csv("R_CYP3A5_v1.csv", header = TRUE)
df1G[c("Allele_1", "Allele_2")][is.na(df1G[c("Allele_1", "Allele_2")])] <- 0
tabl_alle <- read.csv("R_CYP3A5_v1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_CYP3A5_v1_phenotypetable.csv", header = TRUE)
df2G <- df1G %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename("Func_alle1" = "Allele_func.x", "Func_alle2" = "Allele_func.y")
df3G <- c()
for(i in 1:nrow(df2G))
{
  myrow <- df2G[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df3G <- rbind(df3G, unlist(c(myrow, phenotype)))
}
colnames(df3G) <- c('npm_research_id', 'Diplotype', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df3G <- data.frame(df3G)
df4G <- df3G %>% mutate(CYP3A5 = case_when((phenotype == "Normal_metabolizer" | phenotype == "Intermediate_metabolizer") ~ 1, TRUE ~ 0)) %>% select(npm_research_id, CYP3A5)
df5G <- left_join(dem_sub, df3G, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df6G_all <- df5G %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CYP3A5 = n)

df1H <- read.csv("R_CYP4F2_v1.csv", header = TRUE)
df2H <-
  df1H %>% mutate(phenotype = case_when((
    df1H$rs2108622 == 0) ~ "WT",
    (df1H$rs2108622 == 1 | df1H$rs2108622 == 2) ~ "Increased_dose_req",
    TRUE ~ "Untyped"
  )) 
df3H <- df2H %>% mutate(CYP4F2 = case_when(phenotype == "Increased_dose_req" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, CYP4F2)
df4H <- left_join(dem_sub, df2H, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5H_all <- df4H %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_CYP4F2 = n)

df1I <- read.csv("R_DPYD_v1.csv", header = TRUE)
df2I <-
  df1I %>% mutate(DPYD = case_when((
    phenotype == "Intermediate_metabolizer" |
      phenotype == "Poor_metabolizer"
  ) ~ 1,
  TRUE ~ 0
  )) %>% select(npm_research_id, DPYD)
df3I <- left_join(dem_sub, df1I, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df4I_all <- df3I %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_DPYD = n)

df1J <- read.csv("R_F5_v2.csv", header = TRUE)
df2J <-
  df1J %>% mutate(phenotype = case_when(df1J$rs6025 == 1 ~ "Increased_risk",
                                        df1J$rs6025 == 0 ~ "WT",
                                        TRUE ~ "Untyped")) 
df3J <- df2J %>% mutate(F5 = case_when(phenotype == "Increased_risk" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, F5)
df4J <- left_join(dem_sub, df2J, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5J_all <- df4J %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_F5 = n)

df1K <- read.csv("R_G6PD_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_G6PD_v1_alleletable.csv", header = TRUE)
tabl_alle_b <- tabl_alle %>% select(Allele_SG, WHO_class)
tabl_phe <- read.csv("R_G6PD_v1_phenotypetable.csv", header = TRUE)
tabl_phe$Allele_2[is.na(tabl_phe$Allele_2)] <- 0
df2K <- df1K %>% left_join(tabl_alle_b, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle_b, c("Allele_2" = "Allele_SG")) %>% rename("WHO_alle1" = "WHO_class.x", "WHO_alle2" = "WHO_class.y")
df3K <- c()
for(i in 1:nrow(df2K))
{
  myrow <- df2K[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$WHO_alle1 & tabl_phe$Allele_2 == myrow$WHO_alle2),3]  
  df3K <- rbind(df3K, unlist(c(myrow, phenotype)))
}
colnames(df3K) <- c('npm_research_id', 'Diplotype', 'Allele_1', 'Allele_2', 'WHO_alle1', 'WHO_alle2', 'phenotype')
df3K <- data.frame(df3K)
df4K <- df3K %>% mutate(G6PD = case_when((phenotype == "Deficient" | phenotype == "Variable") ~ 1, TRUE ~ 0)) %>% select(npm_research_id, G6PD)
df5K <- left_join(dem_sub, df3K, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df6K_all <- df5K %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_G6PD = n)

df1L <- read.csv("R_HLA-A_v1.csv", header = TRUE)
strip_HLA <- function(x){
  paste(unlist(strsplit(x, ":"))[1], unlist(strsplit(x, ":"))[2], sep = ":")
}
newdf1L <- df1L
for(i in 1:nrow(newdf1L))
{
  newdf1L$Allele_1[i] <- strip_HLA(newdf1L$Allele_1[i])
  newdf1L$Allele_2[i] <- strip_HLA(newdf1L$Allele_2[i])
}
df2L <-
  newdf1L %>% mutate(
    Def_alle1_3101 = case_when(Allele_1 == "A*31:01" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_3101 = case_when(Allele_2 == "A*31:01" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle_3101_score = Def_alle1_3101 + Def_alle2_3101,
    phenotype = case_when((Def_alle_3101_score == 1 | Def_alle_3101_score == 2) ~ "Increased_SCAR",
                          Def_alle_3101_score == 0 ~ "WT",
                          TRUE ~ "Untyped")
  )
df3L <- df2L %>% mutate(HLAA = case_when(phenotype == "Increased_SCAR" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, HLAA)
df4L <- left_join(dem_sub, df2L, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5L_all <- df4L %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_HLAA = n)

df1M <- read.csv("R_HLA-B_v1.csv", header = TRUE)
strip_HLB <- function(x){
  paste(unlist(strsplit(x, ":"))[1], unlist(strsplit(x, ":"))[2], sep = ":")
}
newdf1M <- df1M
for(i in 1:nrow(newdf1M))
{
  newdf1M$Allele_1[i] <- strip_HLB(newdf1M$Allele_1[i])
  newdf1M$Allele_2[i] <- strip_HLB(newdf1M$Allele_2[i])
}
df2M <-
  newdf1M %>% mutate(
    Def_alle1_1502 = case_when(Allele_1 == "B*15:02" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_1502 = case_when(Allele_2 == "B*15:02" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle1_5801 = case_when(Allele_1 == "B*58:01" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_5801 = case_when(Allele_2 == "B*58:01" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle1_5701 = case_when(Allele_1 == "B*57:01" ~ 1, Allele_1 == "Not_typed:NA" ~ -99, TRUE ~ 0),
    Def_alle2_5701 = case_when(Allele_2 == "B*57:01" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0))
df2.2M <- df2M %>% rowwise() %>% mutate(Def_Balle_score_SCAR = sum(Def_alle1_1502, Def_alle2_1502, Def_alle1_5801, Def_alle2_5801),
                                      Def_Balle_score_AHS = sum(Def_alle1_5701, Def_alle2_5701),
                                      phenotype = case_when(Def_Balle_score_SCAR > 0 & Def_Balle_score_AHS > 0 ~ "Increased_SCAR_AHS",
                                                            Def_Balle_score_SCAR > 0 & Def_Balle_score_AHS == 0 ~ "Increased_SCAR",
                                                            Def_Balle_score_SCAR > 0 & Def_Balle_score_AHS < 0 ~ "Increased_SCAR",
                                                            Def_Balle_score_SCAR == 0 & Def_Balle_score_AHS > 0 ~ "Increased_AHS",
                                                            Def_Balle_score_SCAR == 0 & Def_Balle_score_AHS == 0 ~ "WT",
                                                            Def_Balle_score_SCAR == 0 & Def_Balle_score_AHS < 0 ~ "WT",
                                                            Def_Balle_score_SCAR < 0 & Def_Balle_score_AHS > 0 ~ "Increased_AHS",
                                                            Def_Balle_score_SCAR < 0 & Def_Balle_score_AHS == 0 ~ "WT",
                                                            Def_Balle_score_SCAR < 0 & Def_Balle_score_AHS < 0 ~ "Untyped"))
df3M <- df2.2M %>% mutate(HLAB = case_when((phenotype == "Increased_SCAR" | phenotype == "Increased_SCAR_AHS" | phenotype == "Increased_AHS" ~ 1), TRUE ~ 0)) %>% select(npm_research_id, HLAB)
df4M <- left_join(dem_sub, df2.2M, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5M_all <- df4M %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_HLAB = n)

df1N <- read.csv("R_IFNL3_v2.csv", header = TRUE)
df2N <- c()
for(i in 1:nrow(df1N))
{
  myrow <- df1N[i,]
  nasum <- length(which(is.na(myrow[,2:3])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:3]) > 0))
  {
    temp_phenotype <- "Decreased_response"
  } else if(nasum == 0 & (sum(myrow[,2:3]) == 0))
  {
    temp_phenotype <- "WT"
  }
  df2N <- rbind(df2N, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2N) <- c('npm_research_id', 'rs12979860', 'rs8099917', 'NAcount','phenotype')
df2N <- data.frame(df2N)
df3N <- df2N %>% mutate(IFNL3 = case_when(phenotype == "Decreased_response" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, IFNL3)
df4N <- left_join(dem_sub, df2N, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5N_all <- df4N %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_IFNL3 = n)

df1O <- read.csv("R_IFNL4_v2.csv", header = TRUE)
df2O <- c()
for(i in 1:nrow(df1O))
{
  myrow <- df1O[i,]
  nasum <- length(which(is.na(myrow[,2:3])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:3]) > 0))
  {
    temp_phenotype <- "Decreased_response"
  } else if(nasum == 0 & (sum(myrow[,2:3]) == 0))
  {
    temp_phenotype <- "WT"
  }
  df2O <- rbind(df2O, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2O) <- c('npm_research_id', 'rs12979860', 'rs11322783CT.C', 'NAcount','phenotype')
df2O <- data.frame(df2O)
df3O <- df2O %>% mutate(IFNL4 = case_when(phenotype == "Decreased_response" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, IFNL4)
df4O <- left_join(dem_sub, df2O, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5O_all <- df4O %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_IFNL4 = n)

df1P <- read.csv("R_NAT2_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_NAT2_alleletable.csv", header = TRUE)
df2P <- df1P %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename("Func_alle1" = "Allele_func.x", "Func_alle2" = "Allele_func.y")
df3P <- df2P %>% mutate(phenotype = case_when((Func_alle1 == "rapid" & Func_alle2 == "rapid") ~ "Rapid_acetylator",
                                              (Func_alle1 == "slow" & Func_alle2 == "slow") ~ "Slow_acetylator",
                                              (Func_alle1 == "rapid" & Func_alle2 == "slow") ~ "Intermediate_acetylator",
                                              (Func_alle1 == "slow" & Func_alle2 == "rapid") ~ "Intermediate_acetylator",
                                              (Func_alle1 == "rapid" & Func_alle2 == "unknown") ~ "Indeterminate",
                                              (Func_alle1 == "unknown" & Func_alle2 == "rapid") ~ "Indeterminate",
                                              (Func_alle1 == "slow" & Func_alle2 == "unknown") ~ "Indeterminate",
                                              (Func_alle1 == "unknown" & Func_alle2 == "slow") ~ "Indeterminate",
                                              (Func_alle1 == "unknown" & Func_alle2 == "unknown") ~ "Indeterminate",
                                              TRUE ~ "Untyped"))
df4P <- df3P %>% mutate(NAT2 = case_when(phenotype == "Slow_acetylator" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, NAT2)
df5P <- left_join(dem_sub, df3P, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df6P_all <- df5P %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_NAT2 = n)

df1Q <- read.csv("R_NUDT15_v2.csv", header = TRUE)
df2Q <- df1Q %>% mutate(phenotype = case_when(df1Q$rs116855232 == 0 ~ "Normal_metabolizer",
                                              df1Q$rs116855232 == 1 ~ "Intermediate_metabolizer",
                                              df1Q$rs116855232 == 2 ~ "Poor_metabolizer",
                                              TRUE ~ "Untyped"))
df3Q <- df2Q %>% mutate(NUDT15 = case_when((phenotype == "Intermediate_metabolizer" | phenotype == "Poor_metabolizer") ~ 1, TRUE ~ 0)) %>% select(npm_research_id, NUDT15)
df4Q <- left_join(dem_sub, df2Q, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5Q_all <- df4Q %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_NUDT15 = n)

df1R <- read.csv("R_RYR1_v2.csv", header = TRUE)
df2R <- c()
for(i in 1:nrow(df1R))
{
  myrow <- df1R[i,]
  nasum <- length(which(is.na(myrow[,2:4])))
  temp_phenotype <- "Untyped"
  if(nasum == 0 & (sum(myrow[,2:4]) > 0))
  {
    temp_phenotype <- "MHS"
  } else if(nasum == 0 & (sum(myrow[,2:4]) == 0))
  {
    temp_phenotype <- "WT"
  }
  df2R <- rbind(df2R, unlist(c(myrow, nasum, temp_phenotype)))
}
colnames(df2R) <- c('npm_research_id', 'rs112563513', 'rs121918593', 'rs118192168', 'NAcount','phenotype')
df2R <- data.frame(df2R)
df3R <- df2R %>% mutate(RYR1 = case_when(phenotype == "MHS" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, RYR1)
df4R <- left_join(dem_sub, df2R, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5R_all <- df4R %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_RYR1 = n)

df1S <- read.csv("R_SLCO1B1_v2.csv", header = TRUE)
df2S <- df1S %>% mutate(phenotype = case_when(df1S$rs4149056 == 0 ~ "WT",
                                              df1S$rs4149056 == 1 ~ "Intermediate_function",
                                              df1S$rs4149056 == 2 ~ "Low_function", 
                                              TRUE ~"Untyped"))
df3S <- df2S %>% mutate(SLCO1B1 = case_when((phenotype == "Intermediate_function" | phenotype == "Low_function") ~ 1, TRUE ~ 0)) %>% select(npm_research_id, SLCO1B1)
df4S <- left_join(dem_sub, df2S, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df5S_all <- df4S %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_SLCO1B1 = n)

df1T <- read.csv("R_TPMT_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_TPMT_v1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_TPMT_v1_phenotypetable.csv", header = TRUE)
df2T <- df1T %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename(Func_alle1 = Allele_Func.x, Func_alle2 = Allele_Func.y)
df3T <- c()
for(i in 1:nrow(df2T))
{
  myrow <- df2T[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df3T <- rbind(df3T, unlist(c(myrow, phenotype)))
}
colnames(df3T) <- c('npm_research_id', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df3T <- data.frame(df3T)
df4T <- df3T %>% mutate(TPMT = case_when((phenotype == "Intermediate_metabolizer" | phenotype == "Poor_metabolizer") ~ 1, TRUE ~ 0)) %>% select(npm_research_id, TPMT)
df5T <- left_join(dem_sub, df3T, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df6T_all <- df5T %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_TPMT = n)

df1U <- read.csv("R_UGT1A1_v1.csv", header = TRUE)
tabl_alle <- read.csv("R_UGT1A1_alleletable.csv", header = TRUE)
tabl_phe <- read.csv("R_UGT1A1_phenotypetable.csv", header = TRUE)
df2U <- df1U %>% left_join(tabl_alle, c("Allele_1" = "Allele_SG")) %>% left_join(tabl_alle, c("Allele_2" = "Allele_SG")) %>% rename(Func_alle1 = Allele_Func.x, Func_alle2 = Allele_Func.y)
df3U <- c()
for(i in 1:nrow(df2U))
{
  myrow <- df2U[i,]
  phenotype <- tabl_phe[which(tabl_phe$Allele_1 == myrow$Func_alle1 & tabl_phe$Allele_2 == myrow$Func_alle2),3]  
  df3U <- rbind(df3U, unlist(c(myrow, phenotype)))
}
colnames(df3U) <- c('npm_research_id', 'Allele_1', 'Allele_2', 'Func_alle1', 'Func_alle2', 'phenotype')
df3U <- data.frame(df3U)
df4U <- df3U %>% mutate(UGT1A1 = case_when((phenotype == "Intermediate_metabolizer" | phenotype == "Poor_metabolizer") ~ 1, TRUE ~ 0)) %>% select(npm_research_id, UGT1A1)
df5U <- left_join(dem_sub, df3U, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df6U_all <- df5U %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_UGT1A1 = n)

df1V <- read.csv("R_VKORC1_v2.csv", header = TRUE)
df2V <- c()
for(i in 1:nrow(df1V))
{
  myrow <- df1V[i,]
  nasum <- length(which(is.na(myrow[,2:5])))
  df2V <- rbind(df2V, unlist(c(myrow, nasum)))
}
colnames(df2V) <- c('npm_research_id', 'rs7294', 'rs9934438', 'rs9923231', 'rs2359612', 'NAcount')
df2V <- data.frame(df2V)
df3V <- df2V %>% mutate(rs7294_sensitive = case_when(df2V$rs7294 == 0 ~ 1, TRUE ~ 0), 
                        rs9934438_sensitive = case_when(df2V$rs9934438 != 0 ~ 1, TRUE ~ 0),
                        rs9923231_sensitive = case_when(df2V$rs9923231 != 0 ~ 1, TRUE ~ 0),
                        rs2359612_sensitive = case_when(df2V$rs2359612 !=2 ~ 1, TRUE ~ 0)) %>% mutate(sensitive_score = rowSums(.[7:10])) 
df4V <- df3V %>% mutate(phenotype = case_when((NAcount == 0 & sensitive_score > 0) ~ "Sensitive", (NAcount == 0 & sensitive_score == 0) ~ "WT", TRUE ~ "Untyped"))
df5V <- df4V %>% mutate(VKORC1 = case_when(phenotype == "Sensitive" ~ 1, TRUE ~ 0)) %>% select(npm_research_id, VKORC1)
df6V <- left_join(dem_sub, df4V, by = "npm_research_id") %>% select(npm_research_id, genetic_ethnicity, phenotype)
df7V_all <- df6V %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% count(genetic_ethnicity) %>% rename(n.all_VKORC1 = n)

# ===============================================================================
# Section 1.1 : Analysis of carrier distribution for actionable PGX phenotypes
# ===============================================================================
# Consolidate the PGX risk allele carrier master file with basic demographics
indv_dem <- dem %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
templist <- list(df3A, df3B, df4C, df3D, df3E, df4F, df4G, df3H, df2I, df3J, df4K, df3L, df3M, df3N, df3O, df4P, df3Q, df3R, df3S, df4T, df4U, df5V)
masterlist.1 <- Reduce(function(d1,d2) merge(d1,d2, by = "npm_research_id", all.x = TRUE, all.y = FALSE), templist)
masterlist.1a <- left_join(indv_dem, masterlist.1, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M"))
masterlist.1a[is.na(masterlist.1a)] <- 0

# Summarize number of individuals carrying at least one PGX risk allele
masterlist.1b <- masterlist.1a %>% rowwise() %>% mutate(PGX_present = rowSums(across(where(is.numeric))))
# number of individuals with no PGX risk allele:
length(which(masterlist.1b$PGX_present == 0))
# number of individuals with at least one PGX risk allele:
length(which(masterlist.1b$PGX_present > 0))
# median of PGX risk alleles per individual
masterlist.1b$PGX_present %>% median()

# Summarize number of individuals with PGX risk alleles by class of effect (1 = Efficacy, 2 = Efficacy & Toxicity, 3 = Toxicity)
masterlist.1c <-
  masterlist.1b %>% mutate(
    Toxicity = case_when((
      CACNA1S == 1 |
        DPYD == 1 |
        F5 == 1 |
        G6PD == 1 |
        HLAA == 1 |
        HLAB == 1 | NAT2 == 1 | RYR1 == 1 | SLCO1B1 == 1
    ) ~ 3,
    TRUE ~ 0
    ),
    Efficacy = case_when((CFTR == 1 |
                            CYP3A5 == 1 |
                            CYP4F2 == 1 |
                            IFNL3 == 1 | IFNL4 == 1 | VKORC1 == 1) ~ 1,
                         TRUE ~ 0
    ),
    Toxicity_Efficacy = case_when((
      CYP2B6 == 1 |
        CYP2C9 == 1 |
        CYP2C19 == 1 |
        CYP2D6 == 1 |
        NUDT15 == 1 | TPMT == 1 | UGT1A1 == 1
    ) ~ 2,
    TRUE ~ 0
    )
  )
masterlist.1d <- masterlist.1c %>% rowwise() %>% mutate(PGx_type_max = max(c(Toxicity, Efficacy, Toxicity_Efficacy)))
table(masterlist.1d$PGx_type_max)

# Summarize number of individuals with life-threatening toxicity phenotypes (SJS/TEN, DPD toxicity, MHS)
# First, regenerate masterlist to exclude individuals with HLA-B AHS phenotype only
df3M_b <- df2.2M %>% mutate(HLAB = case_when((phenotype == "Increased_SCAR" | phenotype == "Increased_SCAR_AHS" ~ 1), TRUE ~ 0)) %>% select(npm_research_id, HLAB)

templist <- list(df3A, df3B, df4C, df3D, df3E, df4F, df4G, df3H, df2I, df3J, df4K, df3L, df3M_b, df3N, df3O, df4P, df3Q, df3R, df3S, df4T, df4U, df5V)
masterlist.2 <- Reduce(function(d1,d2) merge(d1,d2, by = "npm_research_id", all.x = TRUE, all.y = FALSE), templist)
masterlist.2a <- left_join(indv_dem, masterlist.2, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M"))
masterlist.2a[is.na(masterlist.2a)] <- 0
# Second, collapse genotypes from multiple genes associated with MHS and SJS/TEN, respectively
fatal_tox.2a <-
  masterlist.2a %>% mutate(MHS_positive = case_when((CACNA1S == 1 |
                                                       RYR1 == 1) ~ 1, TRUE ~ 0),
                           SCAR_positive = case_when((HLAA == 1 |
                                                        HLAB == 1) ~ 1, TRUE ~ 0)) 

# number of individuals with risk alleles for MHS
length(which(fatal_tox.2a$MHS_positive > 0))
length(which(fatal_tox.2a$MHS_positive > 0))/9051
# number of individuals with risk alleles HLA-A/HLA-B SJS/TEN
length(which(fatal_tox.2a$SCAR_positive > 0))
length(which(fatal_tox.2a$SCAR_positive > 0))/9051
# number of individuals with risk allele for DPYD
length(fatal_tox.2a$npm_research_id[fatal_tox.2a$DPYD == 1])
length(fatal_tox.2a$npm_research_id[fatal_tox.2a$DPYD == 1])/9051

# Combined list of all individuals with life-threatening toxicity risk alleles
fatal_tox.2b <- masterlist.2a %>% mutate(SAE = case_when((CACNA1S == 1 | RYR1 == 1 | HLAA == 1 | HLAB == 1 | DPYD == 1) ~ 1, TRUE ~ 0))

# Summarize individuals with genotype associated with life-threatening toxicities (SJS/TEN, DPD toxicity, MHS)
# number of individuals with risk alleles for life-threatening toxicities
length(which(fatal_tox.2b$SAE > 0))
# fraction of individuals with risk alleles for life-threatening toxicities
length(which(fatal_tox.2b$SAE > 0))/9051

# ===============================================================================================================================
# Section 2 : Consolidating carrier frequencies by pharmacophenotypes with therapeutic recommendations (Supplementary Data 7)
# ===============================================================================================================================
# Regenerate masterlist with HLA-B SCAR and AHS phenotypes in separate columns
df3M_c <- df2.2M %>% mutate(HLAB_SCAR = case_when((phenotype == "Increased_SCAR" | phenotype == "Increased_SCAR_AHS" ~ 1), TRUE ~ 0),
                            HLAB_AHS = case_when((phenotype == "Increased_AHS" | phenotype == "Increased_SCAR_AHS" ~ 1), TRUE ~ 0)) %>% select(npm_research_id, HLAB_SCAR, HLAB_AHS)

templist <- list(df3A, df3B, df4C, df3D, df3E, df4F, df4G, df3H, df2I, df3J, df4K, df3L, df3M_c, df3N, df3O, df4P, df3Q, df3R, df3S, df4T, df4U, df5V)
masterlist.3 <- Reduce(function(d1,d2) merge(d1,d2, by = "npm_research_id", all.x = TRUE, all.y = FALSE), templist)
masterlist.3a <- left_join(indv_dem, masterlist.3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M"))
masterlist.3a[is.na(masterlist.3a)] <- 0

# Collapse by pharmacophenotype (for multiple genes)
masterlist.byTx_1 <- masterlist.3a %>% mutate(MHS = case_when((CACNA1S == 1 | RYR1 == 1) ~ 1, TRUE ~ 0),
                                              CFTR_ivacaftor_plus = case_when(CFTR == 1 ~1, TRUE ~ 0),
                                              CYP2B6_metab = case_when(CYP2B6 == 1 ~ 1, TRUE ~ 0),
                                              CYP2C9_metab = case_when(CYP2C9 == 1 ~ 1, TRUE ~ 0),
                                              CYP2C19_metab = case_when(CYP2C19 == 1 ~ 1, TRUE ~ 0),
                                              CYP2D6_metab = case_when(CYP2D6 == 1 ~ 1, TRUE ~ 0),
                                              CYP3A5_metab = case_when(CYP3A5 == 1 ~ 1, TRUE ~ 0),
                                              CYP4F2_warfarin_highdose = case_when(CYP4F2 == 1 ~ 1, TRUE ~ 0),
                                              VKORC1_warfarin_sensitivity = case_when(VKORC1 == 1 ~ 1, TRUE ~ 0),
                                              DPD_tox = case_when(DPYD == 1 ~ 1, TRUE ~ 0),
                                              F5_VTE = case_when(F5 == 1 ~ 1, TRUE ~ 0),
                                              G6PD_metab = case_when(G6PD == 1 ~ 1, TRUE ~ 0),
                                              SCAR = case_when((HLAA == 1 | HLAB_SCAR == 1) ~ 1, TRUE ~ 0),
                                              AHS = case_when(HLAB_AHS == 1 ~ 1, TRUE ~ 0),
                                              PEGifn_less = case_when((IFNL3 == 1 | IFNL4 == 1) ~ 1, TRUE ~ 0),
                                              NAT2_isoniazidtox = case_when(NAT2 == 1 ~ 1, TRUE ~ 0),
                                              Thiopurine_less = case_when((NUDT15 == 1 | TPMT == 1) ~ 1, TRUE ~ 0),
                                              SLCO1B1_statintox = case_when(SLCO1B1 == 1 ~ 1, TRUE ~ 0),
                                              UGT1A1_metab = case_when(UGT1A1 == 1 ~ 1, TRUE ~ 0))

masterlist.byTx_2 <- masterlist.byTx_1[, c(3,27:45)]
masterlist.byTx_3 <- setDT(masterlist.byTx_2)[,lapply(.SD, sum), genetic_ethnicity]
masterlist.byTx_4 <- data.frame(t(masterlist.byTx_3))
masterlist.byTx_4 <- masterlist.byTx_4[-1,]
colnames(masterlist.byTx_4) <- c("CH", "MY", "IND")
masterlist.byTx_4 <- tibble::rownames_to_column(masterlist.byTx_4, "Phenotype")
masterlist.byTx_4[, 2:4] <- lapply(masterlist.byTx_4[, 2:4], as.integer)
masterlist.byTx_5 <- masterlist.byTx_4 %>% rowwise() %>% mutate(Overall = (CH+IND+MY))

# Consolidate masterlist of genotyped individuals for all genes
templist_gt.indv <- list(df5A_all, df5B_all, df6C_all, df5D_all, df5E_all, df6F_all, df6G_all, df5H_all, df4I_all, df5J_all, df6K_all, df5L_all, df5M_all, df5N_all, df5O_all, df6P_all, df5Q_all, df5R_all, df5S_all, df6T_all, df6U_all, df7V_all)
masterlist_gt.indv_1 <- Reduce(function(d1,d2) merge(d1,d2, by = "genetic_ethnicity", all.x = TRUE, all.y = FALSE), templist_gt.indv)
masterlist_gt.indv_2 <-
  masterlist_gt.indv_1 %>% rowwise() %>% mutate(
    MHS = max(n.all_CACNA1S, n.all_RYR1),
    SCAR = max(n.all_HLAA, n.all_HLAB),
    AHS = n.all_HLAB,
    PEGifn_less = max(n.all_IFNL3, n.all_IFNL4),
    Thiopurine_less = max(n.all_NUDT15, n.all_TPMT)
  ) %>% select(
    genetic_ethnicity,
    MHS,
    n.all_CFTR,
    n.all_CYP2B6,
    n.all_CYP2C9,
    n.all_CYP2C19,
    n.all_CYP2D6,
    n.all_CYP3A5,
    n.all_CYP4F2,
    n.all_VKORC1,
    n.all_DPYD,
    n.all_F5,
    n.all_G6PD,
    SCAR,
    AHS,
    PEGifn_less,
    n.all_NAT2,
    Thiopurine_less,
    n.all_SLCO1B1,
    n.all_UGT1A1
  ) %>% rename(
    CFTR_ivacaftor_plus = n.all_CFTR,
    CYP2B6_metab = n.all_CYP2B6,
    CYP2C9_metab = n.all_CYP2C9,
    CYP2C19_metab = n.all_CYP2C19,
    CYP2D6_metab = n.all_CYP2D6,
    CYP3A5_metab = n.all_CYP3A5,
    CYP4F2_warfarin_highdose = n.all_CYP4F2,
    VKORC1_warfarin_sensitivity = n.all_VKORC1,
    DPD_tox = n.all_DPYD,
    F5_VTE = n.all_F5,
    G6PD_metab = n.all_G6PD,
    NAT2_isoniazidtox = n.all_NAT2,
    SLCO1B1_statintox = n.all_SLCO1B1,
    UGT1A1_metab = n.all_UGT1A1
  ) %>% adorn_totals()

masterlist_gt.indv_3 <- data.frame(t(masterlist_gt.indv_2[-1]))
colnames(masterlist_gt.indv_3) <- c("ntyped_CH", "ntyped_IND", "ntyped_MY", "ntyped_Overall")
masterlist_gt.indv_3 <- tibble::rownames_to_column(masterlist_gt.indv_3, "Phenotype")

# Generate masterlist of carrier frequency by collapsed pharmacophenotype (Supplementary Data 7)
masterlist_freq.byTx <-
  left_join(masterlist.byTx_5, masterlist_gt.indv_3, by = "Phenotype") %>% rowwise() %>% mutate(
    Freq_CH = (CH / ntyped_CH),
    Freq_IND = (IND / ntyped_IND),
    Freq_MY = (MY / ntyped_MY),
    Freq_Overall = (Overall / ntyped_Overall)
  ) %>% mutate(
    Freq_CH_percent = Freq_CH * 100,
    Freq_IND_percent = Freq_IND * 100,
    Freq_MY_percent = Freq_MY * 100,
    Freq_Overall_percent = Freq_Overall * 100
  )
#write.table(masterlist_freq.byTx, "output_PGX_carrierfreq_byTx.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# =========================================================================================================
# Section 3: Intersection of individuals with CDC Tier1 PLP variants and relevant PGX alleles (Figure 3)
# =========================================================================================================
# To identify individuals with genetic predisposition to ACMG SF v3.0 AD condition and their corresponding PGX risk alleles
# Consolidate individual demographic and variant data
indv_t1 <- left_join(dem, PLP_indv)

# Consolidate individual-level summary of PLP variant carriers by gene for autosomal dominant (AD) genes in ACMG SF v3
indv_acmg73 <- indv_t1 %>% filter(gene_symbol %in% acmg73$Gene[acmg73$Inheritance == "AD"] & !(hgvs_c %in% apoblof$hgvs_c)) %>% group_by(npm_research_id, genetic_sex, genetic_ethnicity, gene_symbol, genotype_code) %>% summarise(Count = n())
indv_masterlist.2 <- indv_acmg73[,c(1:4,6)]
indv_masterlist.2a <- indv_masterlist.2 %>% pivot_wider(names_from = gene_symbol, values_from = Count, values_fill = 0)

# consolidate individual-level ACMG AD genes variant and PGX variant carrier status into one master list
indv_masterlist.3 <- left_join(indv_masterlist.2a, masterlist.1a, by = c("npm_research_id", "genetic_sex", "genetic_ethnicity"))

# identify individuals with ACMG + PGX variant for CDC tier1 conditions (HBOC, Lynch, FH)
cdc_pgx1 <-
  indv_masterlist.3 %>% mutate(
    HBOC = case_when((BRCA1 != 0 | BRCA2 != 0 | PALB2 != 0) ~ 1, TRUE ~ 0),
    LS = case_when((MLH1 != 0 | MSH6 != 0 | PMS2 != 0) ~ 1, TRUE ~ 0),
    FH = case_when((LDLR != 0 | APOB != 0 | PCSK9 != 0) ~ 1, TRUE ~ 0))

cdc_pgx2 <-
  cdc_pgx1 %>% mutate(
    HBOC_PGX = case_when(((BRCA1 != 0 & CYP2D6 != 0) | 
                            (BRCA2 != 0 & CYP2D6 != 0) |
                            (PALB2 != 0 & CYP2D6 != 0)) ~ 1, 
                         TRUE ~ 0),
    LS_PGX = case_when(((MLH1 != 0 & UGT1A1 != 0) |
                          (MSH6 != 0 & UGT1A1 != 0) |
                          (PMS2 != 0 & UGT1A1 != 0)) ~ 1,
                       TRUE ~ 0),
    FH_PGX = case_when(((LDLR != 0 & SLCO1B1 != 0) |
                          (APOB != 0 & SLCO1B1 != 0) |
                          (PCSK9 != 0 & SLCO1B1 != 0)) ~ 1, 
                       TRUE ~ 0))

cdc_pgx3 <- cdc_pgx2 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, BRCA1, BRCA2, PALB2, MLH1, MSH6, PMS2, LDLR, APOB, PCSK9, HBOC, LS, FH, HBOC_PGX, LS_PGX, FH_PGX)

df3F_b <- df3F %>% select(npm_research_id, phenotype) %>% rename(npm_research_id = npm_research_id, phe_CYP2D6 = phenotype)
df3U_b <- df3U %>% select(npm_research_id, phenotype) %>% rename(npm_research_id = npm_research_id, phe_UGT1A1 = phenotype)
df3S_b <- df2S %>% select(npm_research_id, phenotype) %>% rename(npm_research_id = npm_research_id, phe_SLCO1B1 = phenotype)
templist_2 <- list(df3F_b, df3U_b, df3S_b)
cdc_pgx4 <- Reduce(function(d1,d2) merge(d1,d2, by = "npm_research_id", all.x = TRUE, all.y = FALSE), templist_2)
cdc_pgx5 <- left_join(cdc_pgx3, cdc_pgx4, by = "npm_research_id")

# summary of HBOC-predisposed individuals with their CYP2D6 phenotypes
cdc_pgx5_hboc <- cdc_pgx5 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, HBOC, HBOC_PGX, phe_CYP2D6) %>% filter(HBOC == 1) %>% rename(Ancestry = genetic_ethnicity) %>% group_by(Ancestry) %>% count(phe_CYP2D6) %>% pivot_wider(names_from = phe_CYP2D6, values_from = n)
cdc_pgx5_hboc[is.na(cdc_pgx5_hboc)] <- 0
cdc_pgx5_hboc2 <- cdc_pgx5_hboc %>% mutate(total = rowSums(across(where(is.numeric))))
cdc_pgx5_hboc2$Uncallable <- as.numeric(cdc_pgx5_hboc2$Uncallable)
cdc_pgx5_hboc2$total <- as.numeric(cdc_pgx5_hboc2$total)
cdc_pgx5_hboc2$total_typed <- (cdc_pgx5_hboc2$total - cdc_pgx5_hboc2$Uncallable)
cdc_pgx5_hboc2$GT_NA <- (cdc_pgx5_hboc2$Uncallable + cdc_pgx5_hboc2$Indeterminate)
cdc_pgx5_hboc2$Disorder <- "HBOC"
cdc_pgx5_hboc3 <- cdc_pgx5_hboc2 %>% rename(IM = Intermediate_metabolizer, NM = Normal_metabolizer, PM = Poor_metabolizer) %>% select(Disorder, Ancestry, NM, IM, PM, GT_NA) %>% pivot_longer(cols = 3:6, names_to = "PGX_phenotype", values_to = "value")
cdc_pgx5_hboc3

# summary of Lynch-predisposed individuals with their UGT1A1 phenotypes
cdc_pgx5_ls <- cdc_pgx5 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, LS, LS_PGX, phe_UGT1A1) %>% filter(LS == 1) %>% rename(Ancestry = genetic_ethnicity) %>% group_by(Ancestry) %>% count(phe_UGT1A1) %>% pivot_wider(names_from = phe_UGT1A1, values_from = n)
cdc_pgx5_ls[is.na(cdc_pgx5_ls)] <- 0
cdc_pgx5_ls2 <- cdc_pgx5_ls %>% mutate(total = rowSums(across(where(is.numeric))))
cdc_pgx5_ls2$Untyped <- as.numeric(cdc_pgx5_ls2$Untyped)
cdc_pgx5_ls2$total <- as.numeric(cdc_pgx5_ls2$total)
cdc_pgx5_ls2$total_typed <- (cdc_pgx5_ls2$total - cdc_pgx5_ls2$Untyped)
cdc_pgx5_ls2$Disorder <- "LS"
cdc_pgx5_ls3 <- cdc_pgx5_ls2 %>% rename(IM = Intermediate_metabolizer, NM = Normal_metabolizer, PM = Poor_metabolizer) %>% rename(GT_NA = Untyped) %>% select(Disorder, Ancestry, NM, IM, PM, GT_NA) %>% pivot_longer(cols = 3:6, names_to = "PGX_phenotype", values_to = "value")
cdc_pgx5_ls3

# summary of FH-predisposed individuals with their SLCO1B1 phenotypes
cdc_pgx5_fh <- cdc_pgx5 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, FH, FH_PGX, phe_SLCO1B1) %>% filter(FH == 1) %>% rename(Ancestry = genetic_ethnicity) %>% group_by(Ancestry) %>% count(phe_SLCO1B1) %>% pivot_wider(names_from = phe_SLCO1B1, values_from = n)
cdc_pgx5_fh[is.na(cdc_pgx5_fh)] <- 0
cdc_pgx5_fh2 <- cdc_pgx5_fh %>% mutate(total = rowSums(across(where(is.numeric))))
cdc_pgx5_fh2$Untyped <- as.numeric(cdc_pgx5_fh2$Untyped)
cdc_pgx5_fh2$total <- as.numeric(cdc_pgx5_fh2$total)
cdc_pgx5_fh2$total_typed <- (cdc_pgx5_fh2$total - cdc_pgx5_fh2$Untyped)
cdc_pgx5_fh2$Disorder <- "FH"
cdc_pgx5_fh3 <- cdc_pgx5_fh2 %>% rename(IM = Intermediate_function, GT_NA = Untyped, NM = WT) %>% select(Disorder, Ancestry, NM, IM, GT_NA) %>% pivot_longer(cols = 3:5, names_to = "PGX_phenotype", values_to = "value")
cdc_pgx5_fh3

# calculate fraction of individuals the relevant intersect (dual carriers of HBOC-CYP2D6, LS-UGT1A1, FH-SLCO1B1)
cdc_pgx.all.1 <- rbind(cdc_pgx5_hboc3, cdc_pgx5_ls3, cdc_pgx5_fh3)
cdc_pgx.all.2 <- cdc_pgx.all.1 %>% filter(value != 0)
cdc_pgx.all.2$value <- as.numeric(cdc_pgx.all.2$value)
cdc_pgx.all.2 <- as.data.frame(cdc_pgx.all.2)

hboc_all <- sum(cdc_pgx.all.2[which(cdc_pgx.all.2$Disorder == "HBOC"), 4])
hboc_all
hboc_pgx <-
  sum(cdc_pgx.all.2[which((
    cdc_pgx.all.2$Disorder == "HBOC" &
      cdc_pgx.all.2$PGX_phenotype == "IM"
  ) |
    (
      cdc_pgx.all.2$Disorder == "HBOC" &
        cdc_pgx.all.2$PGX_phenotype == "PM"
    )
  ), 4])
hboc_pgx
# fraction of individuals predisposed to HBOC and are CYP2D6 IM/PM:
hboc_pgx/hboc_all

ls_all <- sum(cdc_pgx.all.2[which(cdc_pgx.all.2$Disorder == "LS"), 4])
ls_all
ls_pgx <-
  sum(cdc_pgx.all.2[which((
    cdc_pgx.all.2$Disorder == "LS" &
      cdc_pgx.all.2$PGX_phenotype == "IM"
  ) |
    (
      cdc_pgx.all.2$Disorder == "LS" &
        cdc_pgx.all.2$PGX_phenotype == "PM"
    )
  ), 4])
ls_pgx
# fraction of individuals predisposed to Lynch and are UGT1A1 IM/PM:
ls_pgx/ls_all

fh_all <- sum(cdc_pgx.all.2[which(cdc_pgx.all.2$Disorder == "FH"), 4])
fh_all
fh_pgx <- sum(cdc_pgx.all.2[which(cdc_pgx.all.2$Disorder == "FH" & cdc_pgx.all.2$PGX_phenotype == "IM"), 4])
fh_pgx
# fraction of individuals predisposed to FH and are SLCO1B1 IM/PM:
fh_pgx/fh_all

cdct1_all <- sum(cdc_pgx.all.2$value)
cdct1_all
pgx_all <- sum(hboc_pgx, ls_pgx, fh_pgx)
pgx_all
# fraction of all individuals predisposed to CDC Tier1 genetic conditions and are carriers of PGX risk alleles:
pgx_all/cdct1_all


# To compare the frequency of pharmacophenotype carriers between CDCT1 at-risk vs non-at-risk individuals
#summarise PGX phenotype carrier counts for HBOC & non-HBOC individuals

#subset individual-level PGX status for CYP2D6, UGT1A1, SLCO1B1
pgxphe_a.1 <- masterlist.1a %>% select(npm_research_id, genetic_sex, genetic_ethnicity, CYP2D6, UGT1A1, SLCO1B1)
pgxphe_a.2 <- left_join(pgxphe_a.1, cdc_pgx4)
#subset individual-level PLP carrier status for HBOC, LS, FH genes
cdct1_plp_carriers <- indv_masterlist.2a %>% select(npm_research_id, genetic_sex, genetic_ethnicity, BRCA1, BRCA2, PALB2, MLH1, MSH6, PMS2, APOB, LDLR, PCSK9)
cdct1_plp_carriers.dem <- left_join(dem_sub, cdct1_plp_carriers) %>% replace(is.na(.), 0)
cdct1_plp_carriers.dem.b <- cdct1_plp_carriers.dem %>% mutate(HBOC = case_when((BRCA1 != 0 | BRCA2 != 0 | PALB2 != 0) ~ 1, TRUE ~ 0), 
                                                              LS = case_when((MLH1 != 0 | MSH6 != 0 | PMS2 != 0) ~ 1, TRUE ~ 0), 
                                                              FH = case_when((LDLR != 0 | APOB != 0 | PCSK9 != 0) ~ 1, TRUE ~ 0))
#merge cdct1 status and pgx status dataframes
pgxphe_b.1 <- left_join(cdct1_plp_carriers.dem.b, pgxphe_a.2)

#select hboc
pgxphe_b_hboc.1 <- pgxphe_b.1 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, HBOC, CYP2D6, phe_CYP2D6)
pgxphe_b_hboc.2 <- pgxphe_b_hboc.1 %>% group_by(HBOC, phe_CYP2D6) %>% count() %>% pivot_wider(names_from = phe_CYP2D6, values_from = n) %>% replace(is.na(.), 0)
pgxphe_b_hboc.2b <- pgxphe_b_hboc.2 %>% select(-Uncallable)
pgxphe_b_hboc.2b$HBOC[pgxphe_b_hboc.2b$HBOC == "0"] <- "non_HBOC"
pgxphe_b_hboc.2b$HBOC[pgxphe_b_hboc.2b$HBOC == "1"] <- "HBOC"
pgxphe_b_hboc.2c <- pgxphe_b_hboc.2b %>% mutate(across(.cols = where(is.integer), .fns = as.numeric)) %>% mutate(GTyped = rowSums(across(where(is.numeric))))
pgxphe_b_hboc.2d <- pgxphe_b_hboc.2c %>% mutate(IMPM = (Intermediate_metabolizer + Poor_metabolizer),
                                                UM_frac = (Ultrarapid_metabolizer/GTyped),
                                                NM_frac = (Normal_metabolizer/GTyped),
                                                IM_frac = (Intermediate_metabolizer/GTyped),
                                                PM_frac = (Poor_metabolizer/GTyped),
                                                Indeterminate_frac = (Indeterminate/GTyped),
                                                IMPM_frac = (Intermediate_metabolizer + Poor_metabolizer)/GTyped,
                                                non_IMPM = (GTyped - IMPM))
pgxphe_b_hboc.3 <- pgxphe_b_hboc.2d %>% select(HBOC, IMPM, non_IMPM)
pgxphe_b_hboc.3_data <- pgxphe_b_hboc.3[,-1]
rownames(pgxphe_b_hboc.3_data) <- pgxphe_b_hboc.3$HBOC
fisher.test(pgxphe_b_hboc.3_data)

#select LS
pgxphe_b_ls.1 <- pgxphe_b.1 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, LS, UGT1A1, phe_UGT1A1)
pgxphe_b_ls.2 <- pgxphe_b_ls.1 %>% group_by(LS, phe_UGT1A1) %>% count() %>% pivot_wider(names_from = phe_UGT1A1, values_from = n) %>% replace(is.na(.), 0)
pgxphe_b_ls.2b <- pgxphe_b_ls.2 %>% select(-Untyped)
pgxphe_b_ls.2b$LS[pgxphe_b_ls.2b$LS == "0"] <- "non_LS"
pgxphe_b_ls.2b$LS[pgxphe_b_ls.2b$LS == "1"] <- "LS"
pgxphe_b_ls.2c <- pgxphe_b_ls.2b %>% mutate(across(.cols = where(is.integer), .fns = as.numeric)) %>% mutate(GTyped = rowSums(across(where(is.numeric))))
pgxphe_b_ls.2d <- pgxphe_b_ls.2c %>% mutate(IMPM = (Intermediate_metabolizer + Poor_metabolizer),
                                            NM_frac = (Normal_metabolizer/GTyped),
                                            IM_frac = (Intermediate_metabolizer/GTyped),
                                            PM_frac = (Poor_metabolizer/GTyped),
                                            Indeterminate_frac = (Indeterminate/GTyped),
                                            IMPM_frac = (Intermediate_metabolizer + Poor_metabolizer)/GTyped,
                                            non_IMPM = (GTyped - IMPM))
pgxphe_b_ls.3 <- pgxphe_b_ls.2d %>% select(LS, IMPM, non_IMPM)
pgxphe_b_ls.3_data <- pgxphe_b_ls.3[,-1]
rownames(pgxphe_b_ls.3_data) <- pgxphe_b_ls.3$LS
fisher.test(pgxphe_b_ls.3_data)

#select FH
pgxphe_b_fh.1 <- pgxphe_b.1 %>% select(npm_research_id, genetic_sex, genetic_ethnicity, FH, SLCO1B1, phe_SLCO1B1)
pgxphe_b_fh.2 <- pgxphe_b_fh.1 %>% group_by(FH, phe_SLCO1B1) %>% count() %>% pivot_wider(names_from = phe_SLCO1B1, values_from = n) %>% replace(is.na(.), 0)
pgxphe_b_fh.2b <- pgxphe_b_fh.2 %>% select(-Untyped)
pgxphe_b_fh.2b$FH[pgxphe_b_fh.2b$FH == "0"] <- "non_FH"
pgxphe_b_fh.2b$FH[pgxphe_b_fh.2b$FH == "1"] <- "FH"
pgxphe_b_fh.2c <- pgxphe_b_fh.2b %>% mutate(across(.cols = where(is.integer), .fns = as.numeric)) %>% mutate(GTyped = rowSums(across(where(is.numeric))))
pgxphe_b_fh.2d <- pgxphe_b_fh.2c %>% mutate(IMPM = (Intermediate_function + Low_function),
                                            NM_frac = (WT/GTyped),
                                            IM_frac = (Intermediate_function/GTyped),
                                            PM_frac = (Low_function/GTyped),
                                            IMPM_frac = (Intermediate_function + Low_function)/GTyped,
                                            non_IMPM = (GTyped - IMPM))
pgxphe_b_fh.3 <- pgxphe_b_fh.2d %>% select(FH, IMPM, non_IMPM)
pgxphe_b_fh.3_data <- pgxphe_b_fh.3[,-1]
rownames(pgxphe_b_fh.3_data) <- pgxphe_b_fh.3$FH
fisher.test(pgxphe_b_fh.3_data)

# To generate multiple stacked bar plot (Figure 3) of showing the proportions of pharmacophenotype of individuals predisposed to CDC Tier1 genetic conditions by ancestry
# Figure legends and annotations were added using Adobe Illustrator
#rename PGX_phenotype values so that can reorder
fig3_df <- cdc_pgx.all.2
fig3_hboc <- fig3_df %>% filter(Disorder == "HBOC")
fig3_ls <- fig3_df %>% filter(Disorder == "LS")
fig3_fh <- fig3_df %>% filter(Disorder == "FH")

ancestrycol <- c('#fdae61', '#4575b4', '#d73027')

#HBOC
fig3_hboc.1 <- fig3_hboc %>% mutate(PGX_phenotype = recode(PGX_phenotype, "GT_NA" = "1", "NM" = "2", "IM" = "3", "PM" = "4"))
fig3_hboc.1$Ancestry <- as.character(fig3_hboc.1$Ancestry)
fig3_hboc.1$PGX_phenotype <- as.character(fig3_hboc.1$PGX_phenotype)
fig3_hboc.1b <- fig3_hboc.1 %>% arrange(desc(Ancestry), desc(PGX_phenotype))
fig3_hboc.1b$phenotype <- str_c(fig3_hboc.1b$Ancestry, '_', fig3_hboc.1b$PGX_phenotype)
fig3_hboc.1b$Disorder <- as.factor(fig3_hboc.1b$Disorder)
fig3_hboc.1b$phenotype <- factor(fig3_hboc.1b$phenotype, levels = unique(fig3_hboc.1b$phenotype))
fig3_hboc.1b$Ancestry <- factor(fig3_hboc.1b$Ancestry, levels = unique(fig3_hboc.1b$Ancestry))
barcol_hboc <- c("#0c2c84", "#1d91c0", "#bdbdbd", "#dedede", "#1d91c0", "#bdbdbd", "#dedede", "#0c2c84", "#1d91c0", "#bdbdbd", "#dedede")
hboc_bar1 <- ggplot(fig3_hboc.1b, aes(x = Disorder, y = value, fill = Ancestry)) + geom_bar(position = "stack", stat = "identity") + scale_fill_manual(values = ancestrycol) + theme(panel.background = element_rect(fill = "transparent"))
hboc_bar1
hboc_bar2 <- ggplot(fig3_hboc.1b, aes(x = Disorder, y = value, fill = phenotype)) + geom_bar(position = "stack", stat = "identity", colour = "black", size = 0.05) + scale_fill_manual(values = barcol_hboc) + theme(panel.background = element_rect(fill = "transparent"))
hboc_bar2
hboc_bar1 + hboc_bar2

#Lynch
fig3_ls.1 <- fig3_ls %>% mutate(PGX_phenotype = recode(PGX_phenotype, "GT_NA" = "1", "NM" = "2", "IM" = "3", "PM" = "4"))
fig3_ls.1$Ancestry <- as.character(fig3_ls.1$Ancestry)
fig3_ls.1$PGX_phenotype <- as.character(fig3_ls.1$PGX_phenotype)
fig3_ls.1b <- fig3_ls.1 %>% arrange(desc(Ancestry), desc(PGX_phenotype))
fig3_ls.1b$phenotype <- str_c(fig3_ls.1b$Ancestry, '_', fig3_ls.1b$PGX_phenotype)
fig3_ls.1b$Disorder <- as.factor(fig3_ls.1b$Disorder)
fig3_ls.1b$phenotype <- factor(fig3_ls.1b$phenotype, levels = unique(fig3_ls.1b$phenotype))
fig3_ls.1b$Ancestry <- factor(fig3_ls.1b$Ancestry, levels = unique(fig3_ls.1b$Ancestry))
barcol_ls <- c("#1d91c0", "#bdbdbd", "#1d91c0", "#0c2c84", "#1d91c0", "#bdbdbd", "#dedede")
ls_bar1 <- ggplot(fig3_ls.1b, aes(x = Disorder, y = value, fill = Ancestry)) + geom_bar(position = "stack", stat = "identity") + scale_fill_manual(values = ancestrycol) + theme(panel.background = element_rect(fill = "transparent"))
ls_bar1
ls_bar2 <- ggplot(fig3_ls.1b, aes(x = Disorder, y = value, fill = phenotype)) + geom_bar(position = "stack", stat = "identity", colour = "black", size = 0.05) + scale_fill_manual(values = barcol_ls) + theme(panel.background = element_rect(fill = "transparent"))
ls_bar2
ls_bar1 + ls_bar2

#FH
fig3_fh.1 <- fig3_fh %>% mutate(PGX_phenotype = recode(PGX_phenotype, "GT_NA" = "1", "NM" = "2", "IM" = "3"))
fig3_fh.1$Ancestry <- as.character(fig3_fh.1$Ancestry)
fig3_fh.1$PGX_phenotype <- as.character(fig3_fh.1$PGX_phenotype)
fig3_fh.1b <- fig3_fh.1 %>% arrange(desc(Ancestry), desc(PGX_phenotype))
fig3_fh.1b$phenotype <- str_c(fig3_fh.1b$Ancestry, '_', fig3_fh.1b$PGX_phenotype)
fig3_fh.1b$Disorder <- as.factor(fig3_fh.1b$Disorder)
fig3_fh.1b$phenotype <- factor(fig3_fh.1b$phenotype, levels = unique(fig3_fh.1b$phenotype))
fig3_fh.1b$Ancestry <- factor(fig3_fh.1b$Ancestry, levels = unique(fig3_fh.1b$Ancestry))
barcol_fh <- c("#bdbdbd", "#bdbdbd", "#1d91c0", "#bdbdbd", "#dedede")
fh_bar1 <- ggplot(fig3_fh.1b, aes(x = Disorder, y = value, fill = Ancestry)) + geom_bar(position = "stack", stat = "identity") + scale_fill_manual(values = ancestrycol) + theme(panel.background = element_rect(fill = "transparent"))
fh_bar1
fh_bar2 <- ggplot(fig3_fh.1b, aes(x = Disorder, y = value, fill = phenotype)) + geom_bar(position = "stack", stat = "identity", colour = "black", size = 0.05) + scale_fill_manual(values = barcol_fh) + theme(panel.background = element_rect(fill = "transparent"))
fh_bar2
fh_bar1 + fh_bar2

# =======================================================
# Section 4: Descriptive analysis of novel PGX variants
# =======================================================

#read Tier1B LOF file for PGX genes
pgx_t1b_1a <- read.table("PGX_Tier1BLOF_variants_curated.txt", sep = "\t", header = TRUE)
#subset LOF variants that meet criteria (1. on MANE transcript, and 2. autoPVS1 criterion = NF1/SS1, and 3. not present in CPIC/PharmVar, and 4. on gene for which LOF is mechanism/pgx phenotype outcome is due to LOF)
pgx_t1b_1b <- pgx_t1b_1a %>% filter(Variant_remarks == "met_criteria")
#remove SLCO1B1 NM_006446.5:c.1738C>T due to poor read coverage
pgx_t1b_1c <- pgx_t1b_1b %>% filter(HGVSc != "NM_006446.5:c.1738C>T")
pgx_t1b_1c_n <- length(pgx_t1b_1c$HGVSc)
pgx_t1b_1c_n

#count number of variants by gene
pgx_t1b_1d <- pgx_t1b_1c %>% group_by(SYMBOL) %>% count() %>% arrange(desc(n))
#count number of individuals with variants in CYP2C9, CYP2C19, CYP2D6 genes
pgx_t1b_1c_cyp2c <- pgx_t1b_1c %>% filter(SYMBOL %in% c("CYP2C9", "CYP2C19", "CYP2D6"))
length(pgx_t1b_1c_cyp2c$HGVSc)
pgx_t1b_1c_cyp2c_indv.1 <- pgx_t1b_1c_cyp2c %>% select(HetSamples) %>% separate_rows(HetSamples, sep = ",") %>% distinct(HetSamples)
length(pgx_t1b_1c_cyp2c_indv.1$HetSamples)

#number of variants with SG10K_Health minor allele frequency (MAF) < 1%
pgx_t1b_maf.less1 <- pgx_t1b_1c %>% filter(SG10K_AF < 0.01)
pgx_t1b_maf.less1_n <- length(pgx_t1b_maf.less1$HGVSc)
pgx_t1b_maf.less1_n
#overall percentage
(pgx_t1b_maf.less1_n / pgx_t1b_1c_n)*100
#number of subset variants that are singleton or doubletons
pgx_t1b_maf.less1.b <- pgx_t1b_maf.less1 %>% filter(SG10K_AC <= 2)
pgx_t1b_maf.less1.b_n <- length(pgx_t1b_maf.less1.b$HGVSc)
pgx_t1b_maf.less1.b_n
#percentage over variants with MAF < 1%
(pgx_t1b_maf.less1.b_n / pgx_t1b_maf.less1_n)*100

# identify number of samples for all singleton and doubleton variants
pgx_t1b_maf.less1.b_indv <- pgx_t1b_maf.less1.b %>% select(HetSamples) %>% separate_rows(HetSamples, sep = ",") %>% distinct(HetSamples)
length(pgx_t1b_maf.less1.b_indv$HetSamples)

