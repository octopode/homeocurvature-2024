## Plot lipidomes' mean curvature under multiple estimation schemes

library(tidyverse)
library(here)
library(ggpubr)

source(here("03-scripts", "prep_curvature.R"))
source(here("03-scripts", "plot_helpers.R"))

# same as above but linreg only, formatted as panel for Fig. S7
panel_s9bd = plcidata_ecoli %>% 
  filter(scheme == "linreg") %>% 
  # consolidate headgroup blocks
  group_by(scheme, sid, class) %>%
  summarize(plci = sum(plci)) %>% 
  gg_plcurv(
    #rows = scheme,
    cols = sid,
    y = plci
  ) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x = unit(6, "mm"),
    #axis.text.y = element_text(angle = 90, hjust = 0.5)
  ) +
  labs(
    x = "Strain"
  )
panel_s9bd
ggsave(here("04-suppfigs", "FigS9_Ecoli", "panel_s9bd_20231008a.pdf"), width = 120, height = 60, unit = "mm")
