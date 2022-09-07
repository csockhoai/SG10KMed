# This codes calculates the carrier frequency of pharmacophenotypes for 23 pharmacogenes with high confidence gene-drug interactions, allele frequencies and diplotype frequencies for genes with star allele nomenclature
# Generated outcomes are displayed in :
# a) Table 2 - frequency of pharmacophenotypes
# b) Supplementary Data 6 - allele frequencies
# c) Supplementary Data 11 - diplotype frequencies for genes with star nomenclature
# Code block for each gene is repeated, please note that variables may be overwritten.

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)
library(RVAideMemoire)

# ================================
# Pre-02 : Load main data files
# ================================
# individual-level basic demographic data
indv <- read.table("r5.3_9105_UnrelatedSamples_Demographics_20211117.txt", header = TRUE, sep = "\t")

# =========================================================
# PGX Gene 1: CACNA1S
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_CACNA1S_v2.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(df1$rs772226819 == 1 ~ "MHS",
                                            df1$rs772226819 == 0 ~ "WT",
                                            TRUE ~ "Untyped"))
df3 <- df2 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_a <- left_join(indv, df2, by = "npm_research_id") %>% select(genetic_ethnicity, rs772226819)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# To tabulate number of carriers by zygosity
df4_c <- left_join(indv, df2, by = "npm_research_id") %>% group_by(genetic_ethnicity, rs772226819, phenotype) %>% count() 
df4_d <- df4_c %>% filter(phenotype != "Untyped")
df4_d

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# ================================================================================================
# PGX Gene 2: CFTR
# allele type: rsID, VCF-derived
# individuals with at least 1 untyped genotype in any 1 of the 4 loci is considered "Untyped"
# ================================================================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
# This block of code is looped for each allele with actionable phenotype
df4_a <- left_join(indv, df2, by = "npm_research_id") %>% select(6,9:12)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 3: CYP2B6
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
tabl_alle_phe <- tabl_alle %>% filter(Allele_function %in% c('Decreased_function', 'No_function', 'Increased_function'))
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
df4_a <- df2b %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id")
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(5:ncol(df4_b)) %>% filter(genetic_ethnicity %in% c('C', 'I', 'M'))
# This block of code is looped for each allele with actionable phenotype: Increased/Decreased/No function (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
df5_a <- df2b %>% unite("Diplotype", Allele_1:Allele_2, sep = "/", remove = FALSE) %>% select(npm_research_id, Diplotype, phenotype)
diplo_1 <- left_join(indv, df5_a, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Ultrarapid_metabolizer", "Rapid_metabolizer", "Intermediate_metabolizer", "Poor_metabolizer"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:8,10)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 4: CYP2C9
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_CYP2C9_v1.csv", header = TRUE)
df3 <- df1 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
tabl_alle <- read.table("R_CYP2C9_alleletable.txt", header = TRUE, sep = "\t")
tabl_alle_phe <- tabl_alle %>% filter(Allele_function %in% c('Decreased function', 'No function'))
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
df2 <- df1 %>% separate(Diplotype, into = c('Allele_1', 'Allele_2'), sep = "/")
df4_a <- df2 %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id")
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(5:ncol(df4_b)) %>% filter(genetic_ethnicity %in% c('C', 'I', 'M'))
# This block of code is looped for each allele with actionable phenotype: Decreased/No function (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo_1 <- left_join(indv, df1, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Intermediate_metabolizer ", "Poor_metabolizer "))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:6,8)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 5: CYP2C19
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_CYP2C19_v1.csv", header = TRUE)
tabl_phe <- read.csv("R_CYP2C19_v1_phenotypetable.csv", header = TRUE)
df2 <- df1 %>% left_join(tabl_phe, c("Diplotype_updated" = "Diplotype_CPIC"))
df3 <- df2 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
tabl_alle <- read.csv("R_CYP2C19_v1_alleletable.csv", header = TRUE)
tabl_alle_phe <- tabl_alle %>% filter(Allele_function %in% c('Decreased function', 'No function', 'Increased function'))
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
df2b <- df2 %>% filter(phenotype != "Untyped") %>% separate(Diplotype_updated, into = c('Allele_1', 'Allele_2'), sep = "/")
df4_a <- df2b %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))
# This block of code is looped for each allele with actionable phenotype: Decreased/No function (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo_1 <- left_join(indv, df2, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype_updated, phenotype, genetic_ethnicity) %>% count(Diplotype_updated)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Intermediate_metabolizer", "Poor_metabolizer"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype_updated) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:8,10)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 6: CYP2D6
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

#To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
tabl_alle_phe <- tabl_alle %>% filter(V1_ClinicalAnnotation %in% c('Decreased function', 'No function', 'Increased function'))
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
df4_a <- df2b %>% filter(phenotype != "Uncallable") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))
# This block of code is looped for each allele with actionable phenotype: Increased/Decreased/No function (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo_1 <- left_join(indv, df2b, by = "npm_research_id") %>% unite("Diplotype", Allele_1:Allele_2, sep = "/", remove = FALSE) %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Ultrarapid_metabolizer", "Intermediate_metabolizer", "Poor_metabolizer"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Uncallable") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Uncallable <- as.numeric(t1a_b2$Uncallable)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Uncallable)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:7,9)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# statistical test of CYP2D6 IM+PM across ancestry groups
impm_df1 <- t1_c.T[c(1,3,5,8),c(1:3)]
colnames(impm_df1) <- impm_df1[1,]
impm_df1 <- impm_df1[-1,]
impm_df1 <- as.data.frame(impm_df1)
impm_df1[,1:3] <- lapply(impm_df1[,1:3], as.numeric)
impm_df1_impmsum <- apply(impm_df1[1:2,],2, sum)
impm_df1_residual <- impm_df1[3,] - impm_df1_impmsum
impm_df1_matrix <- rbind.data.frame(impm_df1_impmsum, impm_df1_residual)
impm_df1_matrix <- as.matrix(impm_df1_matrix)
fisher.test(impm_df1_matrix)
fisher.multcomp(impm_df1_matrix, p.method = "BH")

# =========================================================
# PGX Gene 7: CYP3A4
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_CYP3A4_v1.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(Diplotype == "*1/*1" ~ "Normal_metabolism",
                                            (Diplotype == "*1/*22" | Diplotype == "*22/*22") ~ "Decreased_metabolism",
                                            Diplotype == "Untyped" ~ "Untyped",
                                            TRUE ~ "Unknown"))
df3 <- df2 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
tabl_alle_phe <- read.csv("R_CYP3A4_v1_alleletable.csv", header = TRUE)
df4_a <- df2 %>% filter(phenotype != "Untyped") %>% separate(Diplotype, into = c('Allele_1', 'Allele_2'), sep = "/") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))
# This block of code is looped for each allele with actionable phenotype: Decreased function (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo_1 <- left_join(indv, df2, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype == "Decreased_metabolism")
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:5,7)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 8: CYP3A5
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
tabl_alle_phe <- tabl_alle %>% filter(Allele_func == 'No_function')
df4_a <- df2b %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))
# This block of code is looped for each allele with actionable phenotype: No function
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo_1 <- left_join(indv, df2b, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Normal_metabolizer", "Intermediate_metabolizer", "Poor_metabolizer"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:6,8)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 9: CYP4F2
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_CYP4F2_v1.csv", header = TRUE)
df2 <-
  df1 %>% mutate(phenotype = case_when((
    df1$rs2108622 == 0) ~ "WT",
    (df1$rs2108622 == 1 | df1$rs2108622 == 2) ~ "Increased_dose_req",
    TRUE ~ "Untyped"
  )) 
df3 <- df2 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_a <- left_join(indv, df2, by = "npm_research_id") %>% select(genetic_ethnicity, rs2108622)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# tabulate number of carriers by zygosity
df4_c <- left_join(indv, df2, by = "npm_research_id") %>% group_by(genetic_ethnicity, rs2108622, phenotype) %>% count() 
df4_d <- df4_c %>% filter(phenotype != "Untyped")
df4_d

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 10: DPYD
# allele type: rsID/star allele (customized Aldy calling)
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_DPYD_v1.csv", header = TRUE)
df3 <- df1 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
tabl_alle_phe <- read.csv("R_DPYD_deficientalleles.csv", header = TRUE)
tabl_alle_phe_SG <- read.table("R_DPYD_deficientalleles_SG.txt", header = TRUE, sep = "\t")
df4_a <- df1 %>% filter(phenotype != "Untyped") %>% separate(Diplotype, into = c('Allele_1', 'Allele_2'), sep = "/", remove = FALSE) %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))

df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  if(nrow(temp_AC) <1){
    next
  }
  temp_df <- cbind.data.frame(allelename, temp_AC)
  df4_d <- rbind.data.frame(df4_d, temp_df)
}

alle_split <- str_split_fixed(df4_d$allelename, pattern = "\\+", n = 3)
colnames(alle_split) <- c("Alle_1", "Alle_2", "Alle_3")
df4_e <- cbind.data.frame(df4_d[, 2:3], alle_split)
df4_f <- df4_e %>% pivot_longer(cols = 3:5, names_to = "serial", values_to = "Allele") %>% filter(Allele != "") %>% select(genetic_ethnicity, total_AC, Allele)
df4_g <- as.data.frame(sapply(df4_f, gsub, pattern = "\\*rs", replacement = "rs"))
df4_g[df4_g == "*HapB3" | df4_g == "rs75017182"] <- "rs75017182/HapB3"
df4_g[df4_g == "*2" | df4_g == "rs3918290"] <- "*2/rs3918290"
df4_g$total_AC <- as.numeric(df4_g$total_AC)
df4_h <- df4_g %>% select(Allele, genetic_ethnicity, total_AC) %>% group_by(Allele, genetic_ethnicity) %>% summarise(AC = sum(total_AC))
alle_number_1 <- df4_b %>% group_by(genetic_ethnicity) %>% count()
alle_number_2 <- alle_number_1 %>% summarise(GTyped = sum(n)) %>% rowwise() %>% mutate(AN = (GTyped*2))
df5_a <- df4_h %>% filter(Allele %in% tabl_alle_phe_SG$Allele)
df5_b <- left_join(df5_a, alle_number_2, by = "genetic_ethnicity") %>% mutate(AF = AC/AN)
df5_b

# To summarize frequency by diplotype (Supplementary Data 11)
diplo_1 <- left_join(indv, df1, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Intermediate_metabolizer", "Poor_metabolizer"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% mutate(phenotype = ifelse(is.na(phenotype), "Untyped", phenotype)) %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:5,7)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 11: F5
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_F5_v2.csv", header = TRUE)
df2 <-
  df1 %>% mutate(phenotype = case_when(df1$rs6025 == 1 ~ "Increased_risk",
                                       df1$rs6025 == 0 ~ "WT",
                                       TRUE ~ "Untyped")) 
df3 <- df2 %>% select(npm_research_id, phenotype)

#To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_a <- left_join(indv, df2, by = "npm_research_id") %>% select(genetic_ethnicity, rs6025)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# tabulate number of carriers by zygosity
df4_c <- left_join(indv, df2, by = "npm_research_id") %>% group_by(genetic_ethnicity, rs6025, phenotype) %>% count() 
df4_d <- df4_c %>% filter(phenotype != "Untyped")
df4_d

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 12: G6PD
# allele type: rsID, VCF-derived, WHO nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
tabl_alle_phe <- tabl_alle_b %>% filter(WHO_class %in% c("II", "III", "III-IV"))
df4_a <- df2b %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Allele_1, Allele_2)
df4_a[is.na(df4_a)] <- "none"
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_sex, genetic_ethnicity, Allele_1, Allele_2)

df4_m <- df4_b %>% filter(genetic_sex == "M") %>% select(genetic_ethnicity, Allele_1)
df4_f <- df4_b %>% filter(genetic_sex == "F")

for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_f))
  {
    #    print(j)
    temprow <- df4_f[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_f <- cbind(df4_f, temp_alle_list)
  colnames(df4_f)[ncol(df4_f)] <- current_alle_phe
}

for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_m))
  {
    #    print(j)
    temprow <- df4_m[j,]
    temp_alle1 <- temprow$Allele_1
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_m <- cbind(df4_m, temp_alle_list)
  colnames(df4_m)[ncol(df4_m)] <- current_alle_phe
}
df4_f_2 <- df4_f %>% select(2, 5:ncol(df4_f))
df4_m_2 <- df4_m %>% select(1, 3:ncol(df4_m))

alle_num_f <- df4_f_2 %>% group_by(genetic_ethnicity) %>% count() %>% summarise(Carrier_f = sum(n)) %>% rowwise() %>% mutate(AN_f = (Carrier_f*2))
alle_num_m <- df4_m_2 %>% group_by(genetic_ethnicity) %>% count() %>% summarise(Carrier_m = sum(n)) %>% rowwise() %>% mutate(AN_m = (Carrier_m*1))
alle_num <- left_join(alle_num_f, alle_num_m) %>% mutate(AN = sum(AN_f, AN_m)) %>% select(genetic_ethnicity, AN)
# This block of code is looped for each allele with actionable phenotype: WHO II, III, III-IV alleles
df5_a <- rbind.data.frame(df4_f_2, df4_m_2)
df5_b <- data.frame()
for(i in 2:ncol(df5_a))
{
  allelename <- colnames(df5_a)[i]
  output_df <- df5_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_AN <- alle_num
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df5_b <- rbind.data.frame(df5_b, temp_AF)
}
df5_b

# To summarize frequency by diplotype
diplo_1 <- left_join(indv, df2b, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Deficient", "Variable"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b2))
{
  myrow <- t1a_b2[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b2)[1], colnames(t1a_b2[2:ncol(t1a_b2)]), colnames(t1a_b2[2:ncol(t1a_b2)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 13: HLA-A
# alleles called using HLA-HD
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
df4_a <- df2 %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Def_alle_3101_score, phenotype)
df4_b <- left_join(df4_a, indv, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Def_alle_3101_score)
df4_c <- data.frame()
for(i in 2:ncol(df4_b))
{
  allelename <- colnames(df4_b)[i]
  output_df <- df4_b[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_c <- rbind.data.frame(df4_c, temp_AF)
}
df4_c

# To generate phenotype summary count
t1a <- left_join(indv,df3, by = "npm_research_id") %>% mutate(phenotype = ifelse(is.na(phenotype), "Untyped", phenotype)) %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 14: HLA-B
# alleles called using HLA-HD
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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
    Def_alle2_5701 = case_when(Allele_2 == "B*57:01" ~ 1, Allele_2 == "Not_typed:NA" ~ -99, TRUE ~ 0))
df2.2 <- df2 %>% rowwise() %>% mutate(Def_Balle_score_SCAR = sum(Def_alle1_1502, Def_alle2_1502, Def_alle1_5801, Def_alle2_5801),
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
df3 <- df2.2 %>% select(npm_research_id, phenotype) 

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
df2b <- df2.2 %>% mutate(Def_1502 = (Def_alle1_1502 + Def_alle2_1502),
                       Def_5801 = (Def_alle1_5801 + Def_alle2_5801),
                       Def_5701 = (Def_alle1_5701 + Def_alle2_5701))
df4_a <- df2b %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Def_1502, Def_5801, Def_5701, phenotype)
df4_b <- left_join(df4_a, indv, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Def_1502, Def_5801, Def_5701)
df4_c <- data.frame()
for(i in 2:ncol(df4_b))
{
  allelename <- colnames(df4_b)[i]
  output_df <- df4_b[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_c <- rbind.data.frame(df4_c, temp_AF)
}
df4_c

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% mutate(phenotype = ifelse(is.na(phenotype), "Untyped", phenotype)) %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
#c alculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:6,8)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 15: IFNL3
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
df4_a <- df2 %>% select(npm_research_id, rs12979860, rs8099917, phenotype)
df4_b <- left_join(df4_a, indv, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, rs12979860, rs8099917)
df4_c <- data.frame()
for(i in 2:ncol(df4_b))
{
  allelename <- colnames(df4_b)[i]
  output_df <- df4_b[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_c <- rbind.data.frame(df4_c, temp_AF)
}
df4_c

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 16: IFNL4
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
df4_a <- df2 %>% select(npm_research_id, rs12979860, rs11322783CT.C, phenotype)
df4_b <- left_join(df4_a, indv, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, rs12979860, rs11322783CT.C)
df4_c <- data.frame()
for(i in 2:ncol(df4_b))
{
  allelename <- colnames(df4_b)[i]
  output_df <- df4_b[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_c <- rbind.data.frame(df4_c, temp_AF)
}
df4_c

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 17: NAT2
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
tabl_alle_phe <- tabl_alle %>% filter(Allele_func %in% c("rapid", "slow"))
df4_a <- df2b %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))
# This block of code is looped for each allele with actionable phenotype: rapid/slow acetylator
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo_1 <- left_join(indv, df2b, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype == "Slow_acetylator")
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:6,8)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 18: NUDT15
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_NUDT15_v2.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(df1$rs116855232 == 0 ~ "Normal_metabolizer",
                                            df1$rs116855232 == 1 ~ "Intermediate_metabolizer",
                                            df1$rs116855232 == 2 ~ "Poor_metabolizer",
                                            TRUE ~ "Untyped"))
df3 <- df2 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_a <- left_join(indv, df2, by = "npm_research_id") %>% select(genetic_ethnicity, rs116855232)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# tabulate number of carriers by zygosity
df4_c <- left_join(indv, df2, by = "npm_research_id") %>% group_by(genetic_ethnicity, rs116855232, phenotype) %>% count() 
df4_d <- df4_c %>% filter(phenotype != "Untyped")
df4_d

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:5,7)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 19: RYR1
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
# This block of code is looped for each allele with actionable phenotype
df4_a <- left_join(indv, df2, by = "npm_research_id") %>% select(6,9:11)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 20: SLCO1B1
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
df1 <- read.csv("R_SLCO1B1_v2.csv", header = TRUE)
df2 <- df1 %>% mutate(phenotype = case_when(df1$rs4149056 == 0 ~ "WT",
                                            df1$rs4149056 == 1 ~ "Intermediate_function",
                                            df1$rs4149056 == 2 ~ "Low_function", 
                                            TRUE ~"Untyped"))
df3 <- df2 %>% select(npm_research_id, phenotype)

# To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
df4_a <- left_join(indv, df2, by = "npm_research_id") %>% select(genetic_ethnicity, rs4149056)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# tabulate number of carriers by zygosity
df4_c <- left_join(indv, df2, by = "npm_research_id") %>% group_by(genetic_ethnicity, rs4149056, phenotype) %>% count() 
df4_d <- df4_c %>% filter(phenotype != "Untyped")
df4_d

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:5,7)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 21: TPMT
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
tabl_alle_phe <- tabl_alle %>% filter(Allele_Func == "No_function")
df4_a <- df2b %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))
# This block of code is looped for each allele with actionable phenotype: No function
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo <- df2b %>% unite("Diplotype", Allele_1:Allele_2, sep = "/", remove = FALSE)
diplo_1 <- left_join(indv, diplo, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Intermediate_metabolizer", "Poor_metabolizer"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:6,8)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

# =========================================================
# PGX Gene 22: UGT1A1
# allele type: star nomenclature
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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

# To calculate allele count, allele frequency and number of carriers (Supplementary Data 6)
dem <- indv %>% select(npm_research_id, genetic_sex, genetic_ethnicity)
tabl_alle_phe <- tabl_alle %>% filter(Allele_Func %in% c("Decreased_function", "Increased_function"))
df4_a <- df2b %>% filter(phenotype != "Untyped") %>% select(npm_research_id, Allele_1, Allele_2)
df4_b <- left_join(df4_a, dem, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% select(genetic_ethnicity, Allele_1, Allele_2)
for(i in 1:nrow(tabl_alle_phe))
{
  current_alle_phe <- tabl_alle_phe$Allele_SG[i]
  temp_alle_list <- c()
  for(j in 1:nrow(df4_b))
  {
    #    print(j)
    temprow <- df4_b[j,]
    temp_alle1 <- temprow$Allele_1
    temp_alle2 <- temprow$Allele_2
    temp_code <- 0
    
    if(temp_alle1 == current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 2
    }else if(temp_alle1 == current_alle_phe & temp_alle2 != current_alle_phe){
      temp_code <- 1
    }else if(temp_alle1 != current_alle_phe & temp_alle2 == current_alle_phe){
      temp_code <- 1
    }
    temp_alle_list <- c(temp_alle_list, temp_code)
  }
  df4_b <- cbind(df4_b, temp_alle_list)
  colnames(df4_b)[ncol(df4_b)] <- current_alle_phe
}
df4_c <- df4_b %>% select(1,4:ncol(df4_b))
# This block of code is looped for each allele with actionable phenotype: No function
df4_d <- data.frame()
for(i in 2:ncol(df4_c))
{
  allelename <- colnames(df4_c)[i]
  output_df <- df4_c[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  temp1 <- output_df %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_d <- rbind.data.frame(df4_d, temp_AF)
}
df4_d

# To summarize frequency by diplotype (Supplementary Data 11)
diplo <- df2b %>% unite("Diplotype", Allele_1:Allele_2, sep = "/", remove = FALSE)
diplo_1 <- left_join(indv, diplo, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(Diplotype, phenotype, genetic_ethnicity) %>% count(Diplotype)
diplo_2 <- diplo_1 %>% filter(phenotype %in% c("Intermediate_metabolizer", "Poor_metabolizer"))
carriers_1 <- diplo_1 %>% filter(phenotype != "Untyped") %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n))
diplo_df <- data.frame()
for(i in 1:nrow(diplo_2))
{
  myrow <- diplo_2[i,]
  myrow$CarrierCount <- carriers_1$Carrier[which(carriers_1$genetic_ethnicity == myrow$genetic_ethnicity)]
  diplo_df <- rbind.data.frame(diplo_df, myrow)
}
diplo_3 <- diplo_df %>% group_by(Diplotype) %>% mutate(CarrierFreq = (n/CarrierCount))
diplo_3

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:6,8)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

#statistical test of UGT1A1 IM+PM across ancestry groups
impm_df1 <- t1_c.T[c(1,3,5,7),c(1:3)]
colnames(impm_df1) <- impm_df1[1,]
impm_df1 <- impm_df1[-1,]
impm_df1 <- as.data.frame(impm_df1)
impm_df1[,1:3] <- lapply(impm_df1[,1:3], as.numeric)
impm_df1_impmsum <- apply(impm_df1[1:2,],2, sum)
impm_df1_residual <- impm_df1[3,] - impm_df1_impmsum
impm_df1_matrix <- rbind.data.frame(impm_df1_impmsum, impm_df1_residual)
impm_df1_matrix <- as.matrix(impm_df1_matrix)
fisher.test(impm_df1_matrix)
fisher.multcomp(impm_df1_matrix, p.method = "BH")

#statistical test of UGT1A1 PM across ancestry groups
pm_df1 <- t1_c.T[c(1,5,7),c(1:3)]
colnames(pm_df1) <- pm_df1[1,]
pm_df1 <- pm_df1[-1,]
pm_df1 <- as.data.frame(pm_df1)
pm_df1[,1:3] <- lapply(pm_df1[,1:3], as.numeric)
pm_df1_pm <- pm_df1[1,]
pm_df1_residual <- pm_df1[2,] - pm_df1[1,]
pm_df1_matrix <- rbind.data.frame(pm_df1_pm, pm_df1_residual)
pm_df1_matrix <- as.matrix(pm_df1_matrix)
fisher.multcomp(pm_df1_matrix, p.method = "BH")

# =========================================================
# PGX Gene 23: VKORC1
# allele type: rsID, VCF-derived
# =========================================================
# To consolidate genotypes/phenotypes for all individuals
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
df3b <- left_join(indv, df2c, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M"))

# To calculate allele count, allele frequency and number of carriers (genotype legend: 0=wildtype, 1=heterozygous, 2=homozygous)
# This block of code is looped for each allele with actionable phenotype (except rs7294, rs2359612)
df4_a <- left_join(indv, df2c, by = "npm_research_id") %>% select(genetic_ethnicity, rs9934438, rs9923231)
df4_b <- data.frame()
for(i in 2:ncol(df4_a))
{
  allelename <- colnames(df4_a)[i]
  output_df <- df4_a[,c(1,i)]
  colnames(output_df)[2] <- "Allele"
  alle_typed <- output_df %>% filter(!(is.na(Allele)))
  temp1 <- alle_typed %>% group_by(genetic_ethnicity) %>% count(Allele)
  temp2 <- temp1 %>% filter(Allele != 0) %>% group_by(genetic_ethnicity) %>% rowwise() %>% mutate(AC = case_when((Allele == 1) ~ (n*1),
                                                                                                                 (Allele == 2) ~ (n*2)))
  temp_AC <- temp2 %>% group_by(genetic_ethnicity) %>% summarise(total_AC = sum(AC))
  temp_alle_num <- temp1 %>% group_by(genetic_ethnicity) %>% summarise(Carrier = sum(n)) %>% rowwise() %>% mutate(AN = (Carrier*2))
  temp_AN <- temp_alle_num %>% select(genetic_ethnicity, AN)
  temp_df1 <- left_join(temp_AN, temp_AC, by = "genetic_ethnicity")
  temp_df1[is.na(temp_df1)] <- 0
  temp_df2 <- rbind.data.frame(temp_df1, c("Overall", apply(temp_df1[,2:ncol(temp_df1)], 2, sum)))
  temp_df2[,c(2,3)] <- sapply(temp_df2[,c(2,3)], as.numeric)
  temp_AF <- temp_df2 %>% rowwise() %>% mutate(AF = total_AC/AN)
  temp_AF <- cbind.data.frame(rep(allelename, nrow(temp_AF)), temp_AF)
  colnames(temp_AF)[1] <- 'Allele'
  df4_b <- rbind.data.frame(df4_b, temp_AF)
}
df4_b

# for allele rs7294:
df4_c <- left_join(indv, df2c, by = "npm_research_id") %>% select(genetic_ethnicity, rs7294)
df4_d <- df4_c %>% mutate(AC = case_when(rs7294 == 0 ~ 2,
                                         rs7294 == 1 ~ 1,
                                         rs7294 == 2 ~ 0))
df4_e <- df4_d %>% group_by(genetic_ethnicity) %>% filter(!(is.na(AC))) %>% summarise(total_AC = sum(AC))
alle_num <- df4_c %>% filter(!(is.na(rs7294))) %>% group_by(genetic_ethnicity) %>% count() %>% mutate(AN = n*2) %>% select(genetic_ethnicity, AN)
df4_f <- left_join(df4_e, alle_num, by = "genetic_ethnicity") %>% mutate(AF = total_AC/AN)
df4_f

# for allele rs2359612
df4_g <- left_join(indv, df2c, by = "npm_research_id") %>% select(genetic_ethnicity, rs2359612)
df4_h <- df4_g %>% mutate(AC = case_when(rs2359612 == 0 ~ 2,
                                         rs2359612 == 1 ~ 1,
                                         rs2359612 == 2 ~ 0))
df4_i <- df4_h %>% group_by(genetic_ethnicity) %>% filter(!(is.na(AC))) %>% summarise(total_AC = sum(AC))
alle_num <- df4_g %>% filter(!(is.na(rs2359612))) %>% group_by(genetic_ethnicity) %>% count() %>% mutate(AN = n*2) %>% select(genetic_ethnicity, AN)
df4_j <- left_join(df4_i, alle_num, by = "genetic_ethnicity") %>% mutate(AF = total_AC/AN)
df4_j

# To generate phenotype summary count
t1a <- left_join(indv, df3, by = "npm_research_id") %>% filter(genetic_ethnicity %in% c("C", "I", "M")) %>% group_by(genetic_ethnicity) %>% count(phenotype) %>% pivot_wider(names_from = phenotype, values_from = n)
t1a[is.na(t1a)] <- 0
t1a2 <- t1a %>% mutate(total = rowSums(across(where(is.numeric))))
# add summary row
t1a_b2 <- rbind.data.frame(t1a2, c("Total", apply(t1a2[,2:ncol(t1a2)], 2, sum)))
# calculate total genotyped individuals
t1a_b2$Untyped <- as.numeric(t1a_b2$Untyped)
t1a_b2$total <- as.numeric(t1a_b2$total)
t1a_b2$total_typed <- (t1a_b2$total - t1a_b2$Untyped)
t1a_b2.T <- t(t1a_b2)
t1a_b <- t1a_b2[, c(1:4,6)]
# calculate fraction of each phenotype and combine into final table
t1_c <- c()
for(i in 1:nrow(t1a_b))
{
  myrow <- t1a_b[i,]
  tempout <- tempout2 <- c()
  for(j in 2:length(myrow))
  {
    tempout <- c(tempout, as.numeric(myrow[j]) / as.numeric(myrow[length(myrow)]))
    tempout2 <- c(tempout2, as.numeric(myrow[j]))
  }
  outrow <- c(as.character(myrow[1,1]), tempout2, tempout)
  t1_c <- rbind(t1_c, outrow)
}
colnames(t1_c) <- c(colnames(t1a_b)[1], colnames(t1a_b[2:ncol(t1a_b)]), colnames(t1a_b[2:ncol(t1a_b)]))
t1_c.T <- t(t1_c)
t1_c.T

