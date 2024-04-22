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
source(here("03-scripts", "pgls_curvature.R")) # need to split this script? Slow bc jackknife.
source(here("03-scripts", "plot_curvature_soppe.R")) # need to split this script
# so there's a plotting one that's quick to source!

arrangeGrob(
  grobs = list(
    panel_4b,
    panel_4clegend,
    panel_4c
  ),
  widths = c(1.0, 0.7, 0.1, 0.5, 0.9, 1.3),
  #heights = c(1,1),
  layout_matrix = rbind(
    c(NA, 1, NA, 2, 3, 3)#,
    #c(NA, 2, NA)
  )
) %>% 
  ggsave(here("04-mainfigs", "Fig4_homeocurvature", "Fig4_arranged_20240405a.pdf"), plot = ., width = 180, height = 40, units = "mm")

# a little synthesis where we compare the slopes in panels B and C
bind_rows(
  # slope of c0 vs. P?
  curv_so %>% 
    group_by(class) %>% 
    reframe(lm(c0~press, cur_data()) %>% tidy()) %>% 
    filter(term == "press") %>% 
    mutate(id = class),
  # slope of c0 vs. depth
  allmods_plci %>% 
    filter(jknife == "nothing") %>% 
    filter(predvar == "depth_col") %>% 
    filter(scheme == "linreg") %>% 
    filter(subset == "<= 10 deg C") %>% 
    filter(term == "pred") %>% 
    mutate(estimate = -10*estimate) %>%  # fix sign and 10 m/bar
    mutate(id = model)
) %>%
  ggplot(
    aes(
      x = id,
      y = estimate
    )
  ) +
    geom_col() +
    geom_errorbar(
      width = 0.2,
      aes(
        ymin = estimate - std.error,
        ymax = estimate + std.error
      )
    ) +
  theme_tiny() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    x = "Regression",
    y = "Slope (Å^-1 * bar^-1)"
  )
# They are about 5-fold off (with "evolutionary compensation" appearing stronger)
ggsave(here("04-mainfigs", "Fig4_homeocurvature", "SlopeComparison_20230801.pdf"), width = 60, height = 60, units = "mm")
