
me=system("whoami", intern = TRUE)
setwd(paste0("/Users/", me, "/Dropbox/tim3/"))

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(Seurat)

## Set global ggplot2-themes
theme_set(theme_classic(base_size = 17))

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
