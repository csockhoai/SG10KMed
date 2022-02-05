# This code is to generate summary figure of CNV deletions identified by size of deletion and carrier frequency

# Data files required:
# CNVDEL_20211016_plot.txt

library(tidyverse)

# load data files:
# CNV deletion summary (note: AffectedGeneCount = 5 denotes the number of affected genes is greater than 5)
cnv <- read.table("CNVDEL_20211016_plot.txt", header = TRUE, sep = "\t")

cnv$CarrierCount_range <- as.character(cnv$CarrierCount_range)
cnv$CarrierCount_range <- factor(cnv$CarrierCount_range, levels = rev(unique(cnv$CarrierCount_range)))
cnv$AffectedGeneCount <- as.factor(cnv$AffectedGeneCount)

F7 <- ggplot(cnv, aes(x = CarrierCount_range, y = Size_log10, size = AffectedGeneCount, color = AffectedGeneCount)) + geom_point(alpha = 0.5, shape = 16)
F7 + scale_color_manual(values = c("#F01111", "#a06fb9", "#30c030", "#002cff")) + theme(panel.background = element_rect(fill = "transparent"), axis.line = element_line(linetype = "solid", size = 0.5), legend.key = element_blank()) + ylab("CNV size in kb (log10)") + xlab("Number of carriers")
