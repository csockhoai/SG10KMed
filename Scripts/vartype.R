# This code is to tabulate the number of pathogenic (P/LP) variants by consequence type

# Data files required:
# All_Tier1_Pathogenic_Variants_r5.3_20211117_varonly.txt
# Variant_Level_r5.3_20211117.txt
# Z0_var_consequence_class.txt

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)

# ================================
# Pre-02 : Load data files
# ================================
PLP_var <- read.table("All_Tier1_Pathogenic_Variants_r5.3_20211117_varonly.txt", sep = "\t", header = TRUE, quote = "") #disable quotes due to occurring quote marks in clinvar_clndn field
var_class <- read.table("Z0_var_consequence_class.txt", sep = "\t", header = TRUE)

# ==============================================================
# Section 1.0 : Tabulate variant counts by consequence type
# ==============================================================

# re-label variant consequence by defined categories
PLP_var_2 <- left_join(PLP_var, var_class, by = "consequence")
# count number of variants in each variant consequence category
PLP_var_class <- PLP_var_2 %>% group_by(consequence_class) %>% count() %>% rename(var_count = n) %>% arrange(desc(var_count))

# collapsing variant consequences into protein-length altering changes, other SNVs and non-coding groups
PLP_var_class_1 <- PLP_var_class %>% mutate(class = case_when((consequence_class == "Frameshift_indels" | consequence_class == "Nonsense" | consequence_class == "Essential_splice") ~ "PTV", (consequence_class == "Missense" | consequence_class == "Inframe_indels" | consequence_class == "Start-Stop_lost" | consequence_class == "Synonymous") ~ "Oth_SNV", TRUE ~ "Noncoding"))

# calculate fraction of protein-length altering variants
PLP_var_class_2 <- PLP_var_class_1 %>% relocate(class, .after = consequence_class)
PLP_var_class_2b <- PLP_var_class_2[,c(2:3)]
PLP_var_class_2b$var_count <- as.numeric(PLP_var_class_2b$var_count)
PLP_var_class_2c <- PLP_var_class_2b %>% group_by(class) %>% summarise(sumcount = sum(var_count))
totalvar <- sum(PLP_var_class_2c$sumcount)
PLP_var_class_2d <- PLP_var_class_2c %>% rowwise() %>% mutate(var_frac = (sumcount/totalvar))

# barplot of distribution of variant consequence types
PLP_var_class_1$consequence_class <- as.character(PLP_var_class_1$consequence_class)
PLP_var_class_1$consequence_class <- factor(PLP_var_class_1$consequence_class, levels = rev(unique(PLP_var_class_1$consequence_class)))
PLP_var_class_1$class <- as.factor(PLP_var_class_1$class)
mycol <- c('#9DABAB', '#D29B95', '#c25b5f')
names(mycol) <- unique(levels(PLP_var_class_1$class))
F03 <- ggplot(PLP_var_class_1, aes(x=consequence_class, y=var_count, fill=class)) + geom_bar(stat = "identity") + scale_fill_manual(values = mycol) + theme(panel.background = element_rect(fill = "transparent"), axis.line = element_line(linetype = "solid", size = 0.5), axis.text.x = element_text(angle = 45, hjust = 0.95)) + coord_flip()
F03

# Consolidate list of P/LP variants with ancestry level carrier frequency details for Supplementary Data 2
PLP_var_car.freq <- read.table("Variant_Level_r5.3_20211117.txt", sep = "\t", header = TRUE)
PLP_var_2_car.freq <- left_join(PLP_var_2, PLP_var_car.freq)
