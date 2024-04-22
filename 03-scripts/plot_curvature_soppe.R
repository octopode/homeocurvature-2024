## Plot SOPE and SOPPE curvature values determined as part of this study
## They are read in from the same GSheet as literature values.

library(tidyverse)
library(here)

source(here("03-scripts", "plot_helpers.R"))
source(here("03-scripts", "c0_meta_analysis.R"))

curv_so = data_curv %>% 
  filter(class %in% c("PE", "PPE")) %>% 
  filter(
    (carbsn1 == 18) &
      (dbonsn1 == 0) &
      (carbsn2 == 18) &
      (dbonsn2 == 1)
  )

panel_4b = curv_so %>% 
  ggplot(
    aes(
      x = press,
      y = c0,
      ymin = c0 - tol,
      ymax = c0 + tol,
      fill = class,
      color = class,
      group = class
    )
  ) +
  geom_smooth(
    method = "lm", 
    se = FALSE,
    size = 1/.pt
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey25",
    width = 50
  ) +
  geom_point(
    size = 1.5, # ends up covering error bars
    shape = 21,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  scale_color_manual(values = chroma_cl) +
  theme_tiny() +
  theme(
    legend.position = c(0.25, 0.9),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
  ) +
  scale_x_continuous(breaks = seq(0, 2500, 500)) +
  guides(fill = guide_legend(nrow = 2, keyheight = 0.4)) +
  labs(
    x = "Pressure (bar)",
    y = "c0 (Å^-1)"
  )
panel_4b
ggsave(here("04-mainfigs", "FIg4_homeocurvature", "panel_4b_20240407a.pdf"), width = 35, height = 40, unit = "mm")

## SCRATCH ##

