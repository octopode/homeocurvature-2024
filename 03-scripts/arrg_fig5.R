# arrange Fig. 5 for export
# the sourced scripts should be fast so long as their respective prep_*
# scripts aren't sourced downstream

library(tidyverse)
library(here)
library(gridExtra)

# all grobs referenced here come from these scripts
# if one is not found, check to make sure the names match!
source(here("03-scripts", "plot_od600_ecoli.R"))
source(here("03-scripts", "plot_saxs_ecoli.R")) # need to split this script
# so there's a plotting one that's quick to source!

# placeholder
#panel_3b = ggplot() + theme_tiny()

arrangeGrob(
  grobs = list(
    panel_5b,
    panel_5c
  ),
  #widths = c(1,2,3,3),
  heights = c(0.5,1,1,1,1),
  #layout_matrix = rbind(
  #  c(NA,NA,NA,1,1),
  #  c(2,2,NA,NA,NA)
  #)
  layout_matrix = cbind(
    c(NA,NA,2,2,2),
    c(1,1,1,NA,NA)
  )
) %>% 
ggsave(here("04-mainfigs", "Fig5_ecolimodel", "Fig5_arranged_20230802b.pdf"), plot = ., width = 121, height = 70, units = "mm")

arrangeGrob(
  grobs = list(
    panel_5b,
    panel_5c
  ),
  widths = c(1, 1),
  heights = c(0.4, 0.4, 1),
  layout_matrix = cbind(
    c(NA,2,2),
    c(1,1,NA)
  )
) %>% 
  ggsave(here("04-mainfigs", "Fig5_ecolimodel", "Fig5_arranged_20230802d.pdf"), plot = ., width = 121, height = 70, units = "mm")

# changing the arrangement
arrangeGrob(
  grobs = list(
    panel_5b,
    panel_5c
  ),
  widths = c(1, 1),
  heights = c(1,2),
  layout_matrix = rbind(
    c(NA,1),
    c(NA,2)
  )
) %>% 
  ggsave(here("04-mainfigs", "Fig5_ecolimodel", "Fig5_arranged_20230802e.pdf"), plot = ., width = 121, height = 60, units = "mm")

# single column!
# changing the arrangement
arrangeGrob(
  grobs = list(
    panel_5b,
    panel_5d
  ),
  heights = c(0.6,0.6,0.85,0.45),
  layout_matrix = cbind(
    c(NA,1,NA,2)
  )
) %>% 
  ggsave(here("04-mainfigs", "Fig5_ecolimodel", "Fig5_arranged_20230803b.pdf"), plot = ., width = 57, height = 110, units = "mm")
