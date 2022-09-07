# This code is to generate summary figure of CNV deletions identified by size of deletion and carrier frequency

# Data files required:
# CNVDEL_20211016_plot.txt
# CNVDEL_summary_20211117_final.txt
# SMA_CarrierStatus.txt

# ===================================
# Pre-01 : Load required packages
# ===================================
library(tidyverse)

# ================================
# Pre-02 : Load data files
# ================================
# SMA deletion carrier summary (only processed for 30X WGS samples)
sma_del <- read.table("SMA_CarrierStatus.txt", sep = "\t", header = TRUE)
# CNV deletion summary with ancestry group breakdown
cnv_del <- read.table("CNVDEL_summary_20211117_final.txt", sep = "\t", header = TRUE)
# CNV deletion summary to plot (note: AffectedGeneCount = 5 denotes the number of affected genes is greater than 5)
cnv_plot <- read.table("CNVDEL_20211016_plot.txt", header = TRUE, sep = "\t")

# ==================================================
# Section 1 : Summary of SMA deletion carriers
# ==================================================

sma_alltested <- length(sma_del$Sample)
sma_alltested
sma_carrier <- length(sma_del$Sample[sma_del$isCarrier == "TRUE"])
sma_carrier
sma_carrier_frac <- sma_carrier/sma_alltested
sma_carrier_frac

# ====================================================================================
# Section 2 : Carrier distribution of CNV deletions by size of deletion event
# ====================================================================================
# Consolidate CNV deletion events with overall carrier frequency > 0.1%
cnvdel_1 <- cnv_del %>% mutate(CF_Total = (Total/9051), max_frac = pmax(CF_C, CF_I, CF_M))
cnvdel_1b <- cnvdel_1 %>% filter(CF_Total >= 0.001) %>% arrange(desc(CF_Total))
cnvdel_1b_c <- cnvdel_1 %>% filter(CF_C >= 0.001) %>% arrange(desc(CF_C))
cnvdel_1b_i <- cnvdel_1 %>% filter(CF_I >= 0.001) %>% arrange(desc(CF_I))
cnvdel_1b_m <- cnvdel_1 %>% filter(CF_M >= 0.001) %>% arrange(desc(CF_M))

# Plot distribution of CNV deletions by size and number of carriers
cnv_plot$CarrierCount_range <- as.character(cnv_plot$CarrierCount_range)
cnv_plot$CarrierCount_range <- factor(cnv_plot$CarrierCount_range, levels = rev(unique(cnv_plot$CarrierCount_range)))
cnv_plot$AffectedGeneCount <- as.factor(cnv_plot$AffectedGeneCount)

F7 <- ggplot(cnv_plot, aes(x = CarrierCount_range, y = Size_log10, size = AffectedGeneCount, color = AffectedGeneCount)) + geom_point(alpha = 0.5, shape = 16)
F7 + scale_color_manual(values = c("#F01111", "#a06fb9", "#30c030", "#002cff")) + theme(panel.background = element_rect(fill = "transparent"), axis.line = element_line(linetype = "solid", size = 0.5), legend.key = element_blank()) + ylab("CNV size in kb (log10)") + xlab("Number of carriers")
