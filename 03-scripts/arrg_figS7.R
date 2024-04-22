# arrange Fig. 4 for export
# the sourced scripts should be fast so long as their respective prep_*
# scripts aren't sourced downstream

library(tidyverse)
library(here)
library(gridExtra)

# all grobs referenced here come from these scripts
# if one is not found, check to make sure the names match!
#source(here("03-scripts", "plot_lipidmaps.R"))
source(here("03-scripts", "plot_curvature_cteno.R"))
source(here("03-scripts", "pgls_curvature.R")) # need to split this script

arrangeGrob(
  grobs = list(
    panel_s7c,
    panel_s7abcd
  ),
  #widths = c(1.8, 1.2),
  heights = c(0.5,1),
  #layout_matrix = rbind(
  #  c(1, 2)
  #  #c(NA, 2, NA)
  #)
) %>% 
  ggsave(here("04-suppfigs", "FigS5_barc0", "FigS5_arranged_20240407a.pdf"), plot = ., width = 120, height = 180, units = "mm")

