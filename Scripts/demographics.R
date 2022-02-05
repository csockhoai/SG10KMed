# This code is to tabulate basic demographics and overall pathogenic (PLP) variant carriers in Supplementary Table 2

# Data files required:
# Full_Sample_Info_9051_edit_v2.txt
# All_Tier1_Pathogenic_Variants_r5.3_20211117.txt

library(tidyverse)

# load data files:
# individual-level cohort demographic data
dem <- read.table("Full_Sample_Info_9051_edit_v2.txt", header = TRUE, sep = "\t")
# individual-level PLP variant carrier status data
PLP_indv <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117.txt", header = TRUE, sep = "\t")

dem %>% group_by(genetic_ethnicity, genetic_sex) %>% count()
#age distribution (median, range[min-max]) of overall cohort
dem_overall <- dem %>% replace(is.na(.), 0) %>% summarise(
  Age_median = median(age),
  Age_min = min(age),
  Age_max = max(age)
)
dem_overall

#age distribution (median, range[min-max]) for each ancestry group
dem_ancestry <- dem %>% group_by(genetic_ethnicity) %>% replace(is.na(.), 0) %>% summarise(
  Age_median = median(age),
  Age_min = min(age),
  Age_max = max(age)
)
dem_ancestry

# cohort distribution by age bracket
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

# list of individuals with pathogenic/likely pathogenic (PLP) variants identified
PLP_indv.count <- PLP_indv %>% group_by(npm_research_id) %>% count()
count(distinct(PLP_indv.count))

# tabulate PLP carriers by ancestry group
PLP_dem <- left_join(dem, PLP_indv.count)
PLP_dem.ancestry <- PLP_dem %>% group_by(genetic_ethnicity) %>% filter(n >0) %>% count()
PLP_dem.ancestry
