# arrange Fig. 3 for export
# the sourced scripts should be fast so long as their respective prep_*
# scripts aren't sourced downstream

library(tidyverse)
library(here)
library(gridExtra)

# all grobs referenced here come from these scripts
# if one is not found, check to make sure the names match!
#source(here("03-scripts", "plot_lipidmaps.R"))
source(here("03-scripts", "plot_pgls_lipidmaps.R")) # need to split this script
source(here("03-scripts", "plot_md.R"))
# so there's a plotting one that's quick to source!

# placeholder
#panel_3b = ggplot() + theme_tiny()

arrangeGrob(
  grobs = list(
    panel_3a_bubble,
    panel_3b,
    panel_3cd,
    panel_3efgh
  ),
  #widths = c(1,2,3,3),
  heights = c(
    0.45, 
    0.05, # buffer
    0.15, # buffer
    0.75
  ),
  layout_matrix = rbind(
    c(rep(1, 4), rep(2, 6)),
    #c(rep(NA, 4), rep(NA, 6)), # make panel A shorter
    c(rep(1, 4), rep(NA, 6)),
    c(rep(NA,6), rep(4, 4)),
    c(rep(3, 6), rep(4, 4))#,
    #c(4,4,4,4,4),
    #c(NA,NA,5,5,5)
  )
) %>% 
ggsave(here("04-mainfigs", "Fig3_ctenolipidomics", "Fig3_arranged_20231002c.pdf"), plot = ., width = 121, height = 90, units = "mm")

