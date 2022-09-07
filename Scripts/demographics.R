# This code is to tabulate basic demographics and overall pathogenic (P/LP) variant carriers in Supplementary Table 2

# Data files required:
# Full_Sample_Info_9051_edit_v2.txt
# All_Tier1_Pathogenic_Variants_r5.3_20211117.txt

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)

# ================================
# Pre-02 : Load data files
# ================================
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
# individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")

# =================================================================
# Section 01: Cohort distribution - genetic sex, genetic ancestry
#==================================================================
# summarizing cohort-level distribution by genetic ancestry and genetic sex
dem %>% group_by(genetic_ethnicity, genetic_sex) %>% count()

# =================================================================
# Section 02: Cohort distribution - age
#==================================================================
# summarizing age distribution (median, range[min-max]) of overall cohort
dem_overall <- dem %>% replace(is.na(.), 0) %>% summarise(
  Age_median = median(age),
  Age_min = min(age),
  Age_max = max(age)
)
dem_overall

# summarizing age distribution (median, range[min-max]) for each ancestry group
dem_ancestry <- dem %>% group_by(genetic_ethnicity) %>% filter(!(is.na(age))) %>% summarise(
  Age_median = median(age),
  Age_min = min(age),
  Age_max = max(age)
)
dem_ancestry

# summarizing cohort distribution by age bracket
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
dem_age.band

# displaying distribution of cohort age by genetic ancestry and testing for difference in age distribution
#subset demographic data by genetic ancestry and age, and exclude age=NA individuals
dem_ancestry_1 <- dem %>% select(genetic_ethnicity, age) %>% filter(!is.na(age))
#statistical test using nonparametric one-way ANOVA - Kruskal-Wallis test
age_ch <- dem_ancestry_1 %>% filter(genetic_ethnicity == "C") %>% select(age)
age_ind <- dem_ancestry_1 %>% filter(genetic_ethnicity == "I") %>% select(age)
age_my <- dem_ancestry_1 %>% filter(genetic_ethnicity == "M") %>% select(age)
kruskal.test(list(age_ch$age, age_ind$age, age_my$age))
#statistical test using ANOVA followed by Tukey test
table(dem_ancestry_1$genetic_ethnicity)
myaov <- aov(dem_ancestry_1$age~dem_ancestry_1$genetic_ethnicity)
summary(myaov)
TukeyHSD(myaov)
plot(density(dem_ancestry_1$age[which(dem_ancestry_1$genetic_ethnicity == "C")]), col = "red")
lines(density(dem_ancestry_1$age[which(dem_ancestry_1$genetic_ethnicity == "I")]), col = "green", add = T)
lines(density(dem_ancestry_1$age[which(dem_ancestry_1$genetic_ethnicity == "M")]), col = "blue", add = T)

# =================================================================
# Section 03: P/LP variant distribution by genetic ancestry 
#==================================================================
# list of individuals with pathogenic/likely pathogenic (P/LP) variants identified
PLP_indv.count <- PLP_indv %>% group_by(npm_research_id) %>% count()
count(distinct(PLP_indv.count))

# number of individuals with pathogenic/likely pathogenic (P/LP) variants identified
length(PLP_indv.count$npm_research_id)

# tabulate P/LP carriers by ancestry group
PLP_dem <- left_join(dem, PLP_indv.count)
PLP_dem.ancestry <- PLP_dem %>% group_by(genetic_ethnicity) %>% filter(n >0) %>% count()
PLP_dem.ancestry
