library(tidyverse)

#To read dataset file: individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")

#To calculate basic demographics (sex, ancestry) for Supplementary Table 2
dem %>% group_by(genetic_ethnicity, genetic_sex) %>% count()

#To calculate age distribution (median, range[min-max]) : overall cohort
median(dem$age, na.rm = TRUE)
min(dem$age, na.rm = TRUE)
max(dem$age, na.rm = TRUE)

#To calculate age distribution (median, range[min-max]) : for each ancestry group
#The following block of code (Block1) calculates for one ancestry group
#Please modify 'genetic_ethnicity' variable (Chinese: 'C', Indian: 'I', Malay: 'M') to calculate age distribution for each ancestry group
###start_of_Block1###
dem_ancestry <- dem %>% group_by(genetic_ethnicity, age) %>% filter(genetic_ethnicity == 'C')
median(dem_ancestry$age, na.rm = TRUE)
min(dem_ancestry$age, na.rm = TRUE)
max(dem_ancestry$age, na.rm = TRUE)
###end_of_Block1###

#To tabulate the counts of individuals in the cohort by age bracket in Table S2
dem_age.band <- dem %>% mutate(age_band = case_when((age < 10) ~ "birth-9y",
                                                      (age >= 10 & age < 20) ~ "10-19y",
                                                      (age >= 20 & age < 30) ~ "20-29y",
                                                      (age >= 30 & age < 40) ~ "30-39y",
                                                      (age >= 40 & age < 50) ~ "40-49y",
                                                      (age >= 50 & age < 60) ~ "50-59y",
                                                      (age >= 60 & age < 70) ~ "60-69y",
                                                      (age >= 70 & age < 80) ~ "70-79y",
                                                      (age >= 80 & age < 90) ~ "80-89y",
                                                      (age >= 90 & age < 100) ~ "90-99y")) %>% group_by(genetic_ethnicity, age_band) %>% count()

#To tabulate the number of individuals with a pathogenic/likely pathogenic (PLP) variant in Supplementary Table 2
#First, read dataset file: individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")
PLP_indv1 <- PLP_indv %>% group_by(npm_research_id) %>% count()
count(distinct(PLP_indv1))

#To tabulate number of individuals with PLP variants by ancestry group in Supplementary Table 2
dem_base <- dem %>% select(npm_research_id,genetic_sex,genetic_ethnicity,age,study,target_depth)
PLP_dem <- left_join(dem, PLP_indv1)
PLP_dem_ancestry <- PLP_dem %>% group_by(genetic_ethnicity) %>% filter(n >0) %>% count()
View(PLP_dem_ancestry)
