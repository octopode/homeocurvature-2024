## Plot results of MD simulations

library(tidyverse)
library(here)
library(httr)

source(here("03-scripts", "plot_helpers.R"))

file_md = here("01-rawdata", "md_results.csv")
# this is a Google Sheets key for the curvature info file, to facilitate a nice Excel-like interface
gsht_md = "1MvKmItEAIrCP9TnqtK-WCYajjjIdhiRAS1ZKh_fWVbE"

# refresh c0 spreadsheet from online interface
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_md}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_md)

results_md = file_md %>% read_csv()

# dF/dR + dAPL plots for bottom row of Fig. 3
slide_mdpars = results_md %>% 
  #filter(sp %in% c("Boli_infu", "Tjal_pink")) %>% 
  # melt down the values and errors using .value sentinel
  pivot_longer(-c(sp, press, temp), names_to = c("var", ".value"), names_sep = '_') %>% 
  ggplot(
    aes(
      x = press,
      y = val,
      ymin = val - err,
      ymax = val + err,
      color = sp,
      fill = sp,
      group = sp
    )
  ) +
  facet_wrap(~var, scale = "free_y", strip.position = "left", nrow = 1) +
  geom_errorbar(
    color = "grey75",
    size = 0.5/.pt,
    width = 50,
    position = position_dodge(width = 50),
  ) +
  geom_smooth(
    data = . %>% filter(var == "dapl_var"),
    method = "lm",
    formula = y~x+0, # force the intercept
    se = FALSE,
    size = 1/.pt,
    position = position_dodge(width = 50)
  ) +
  geom_line(
    data = . %>% filter(var == "dfdr_val"),
    size = 1/.pt,
    position = position_dodge(width = 50)
  ) +
  geom_point(
    size = 3,
    shape = 21,
    color = "white",
    position = position_dodge(width = 50)
  ) +
  scale_fill_manual(values = chroma_sp) +
  scale_color_manual(values = chroma_sp) +
  theme_pubk() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    # replacing the y label with the strip text
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none"
  ) +
  guides(
    linetype = "none"
  ) +
  labs(
    x = "Pressure (bar)",
    y = "Response value",
    size = "N"
  )
slide_mdpars
ggsave(here("04-slidefigs", "mdpars_allspp_20230928a.pdf"), width = 7, height = 3)

## slide figures
panel_3ij = results_md %>% 
  #filter(sp %in% c("Boli_infu", "Tjal_pink")) %>% 
  # melt down the values and errors using .value sentinel
  pivot_longer(-c(sp, press, temp), names_to = c("var", ".value"), names_sep = '_') %>% 
  # order the vars
  mutate(var = factor(var, levels = c("dapl", "difn", "dfdr"))) %>% 
  ggplot(
    aes(
      x = press,
      y = val,
      ymin = val - err,
      ymax = val + err,
      color = sp,
      fill = sp,
      group = sp
    )
  ) +
  #facet_wrap(~var, scale = "free_y", strip.position = "left", nrow = 1) +
  facet_grid(rows = vars(var), scales = "free_y", switch = "y") +
  geom_errorbar(
    color = "grey25",
    size = 0.5/.pt,
    width = 50,
    position = position_dodge(width = 50),
  ) +
  geom_smooth(
    data = . %>% filter(var == "dapl_var"),
    method = "lm",
    formula = y~x+0, # force the intercept
    se = FALSE,
    size = 1/.pt,
    position = position_dodge(width = 50)
  ) +
  # these linregs look like hell
  #geom_smooth(
  #  data = . %>% filter(var == "dfdr"),
  #  method = "lm",
  #  formula = y~x, # don't force the intercept
  #  se = FALSE,
  #  size = 1/.pt
  #) +
  geom_line(
    data = . %>% filter(var == "dfdr_val"),
    size = 1/.pt,
    position = position_dodge(width = 50)
  ) +
  geom_point(
    size = 1.5,
    shape = 21,
    color = "white",
    position = position_dodge(width = 50)
  ) +
  scale_fill_manual(values = chroma_sp) +
  scale_color_manual(values = chroma_sp) +
  theme_tiny() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    # replacing the y label with the strip text
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none"
  ) +
  guides(
    linetype = "none"
  ) +
  labs(
    x = "Pressure (bar)",
    y = "Response value",
    size = "N"
  )
panel_3ij

# Barplots for Panel S7B
results_md %>% 
  #filter(sp %in% c("Boli_infu", "Tjal_pink")) %>% 
  # melt down the values and errors using .value sentinel
  pivot_longer(-c(sp, press, temp), names_to = c("var", ".value"), names_sep = '_') %>% 
  # order the vars
  mutate(var = factor(var, levels = c("dapl", "difn", "dfdr"))) %>%
  # order the spp
  mutate(sp = factor(sp, levels = c("Boli_vitr", "Boli_infu", "Lamp_crue", "Tjal_pink"))) %>% 
  ggplot( 
    aes(
      x = sp,
      y = val,
      fill = sp
    )
  ) +
  geom_hline(yintercept = 0, size = 0.5/.pt) +
  geom_col() +
  geom_errorbar(
    color = "grey25",
    aes(
      ymin = val-err,
      ymax = val+err
    ),
    width = 0.2
  ) +
  facet_grid(
    rows = vars(var),
    cols = vars(press),
    scales = "free_y",
    switch = 'x'
  ) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background.x = element_blank(),
    panel.background = element_rect(fill = "#D9D9D9")
  ) +
  scale_fill_manual(values = chroma_sp) +
  labs(
    x = "Pressure (bar)"
  )
ggsave(here("04-suppfigs", "FigS7_ComplexSims", "panel_s7b_20231007a.pdf"), width = 75, height = 120, units = "mm")

# do a little interpolation magic to find the pressure
# where the two dF/dRs are equal. Implementation is overkill but extensible.
interp_md = results_md %>% 
  group_by(sp) %>% 
  summarise(fun = approxfun(dfdr_val, press) %>% list())#

# _inverse_ interp function (y to x)
tjal_mod = interp_md %>% 
  filter(sp == "Tjal_pink") %>% 
  .$fun %>% .[[1]]

panel_4d_old = results_md %>% 
  ggplot(
    aes(
      x = press,
      y = dfdr_val,
      ymin = dfdr_val - dfdr_err,
      ymax = dfdr_val + dfdr_err,
      fill = sp,
      color = sp,
      group = sp
    )
  ) +
  # draw lil alignment rect
  geom_path(
    data = results_md %>% filter((sp =="Boli_infu") & (press == min(press))) %>% 
      {bind_rows(
        mutate(., press = -200),
        mutate(., press = tjal_mod(dfdr_val)),
        mutate(., press = tjal_mod(dfdr_val), dfdr_val = 0)
      )},
    linetype = "dashed",
    size = 0.5/.pt,
    color = "black"
  ) +
  geom_line(
    size = 1/.pt
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey25",
    width = 50
  ) +
  geom_point(
    size = 1,
    shape = 21,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_sp) +
  scale_color_manual(values = chroma_sp) +
  theme_tiny() +
  theme(
    legend.position = c(0.6, 0.8),
    legend.title = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2, keyheight = 0.4)) +
  # "true zoom"
  # https://stackoverflow.com/questions/25685185/limit-ggplot2-axes-without-removing-data-outside-limits-zoom
  coord_cartesian(
    xlim = c(0, NA),
    ylim = c(0.075, NA)
  ) +
  labs(
    x = "Pressure (bar)",
    y = "dF/dR"
  )
panel_4d_old
#ggsave(here("04-mainfigs", "Fig4_curvature", "panel_4d_20230728a.pdf"), width = 40, height = 40, unit = "mm")

## lateral pressure profiles!
file_pp = here("01-rawdata", "md_latpp.csv")
# this is a Google Sheets key for the curvature info file, to facilitate a nice Excel-like interface
gsht_pp = "1-lVDrejO8zvpj-nHzBSQwyRoXvYavVIkdNZBFJCJsbo"

# refresh c0 spreadsheet from online interface
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_pp}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_pp)

results_pp = file_pp %>% read_csv()

# panel S7C
results_pp %>% 
  arrange(sp, press, z) %>% 
  ggplot(
    aes(
      y = press_lat,
      x = z,
      color = press,
      group = press
    )
  ) +
  geom_line() +
  facet_wrap(~sp, nrow = 1) +
  coord_flip()

# with ambient pressure offset
results_pp %>% 
  # crop
  filter(between(z, -35, 35)) %>% 
  arrange(sp, press, z) %>% 
  # order the spp
  mutate(sp = factor(sp, levels = c("Boli_vitr", "Boli_infu", "Lamp_crue", "Tjal_pink"))) %>% 
  ggplot(
    aes(
      y = press_lat-press,
      x = z,
      color = press,
      group = -press
    )
  ) +
  geom_line() +
  #geom_smooth(span = 0.05, se = FALSE, size = 0.5) +
  facet_wrap(~sp, nrow = 1) +
  coord_flip() +
  ylim(c(-600, 300)) +
  scale_color_viridis_c(direction = -1) +
  theme_tiny() +
  theme(legend.position = "none")
ggsave(here("04-suppfigs", "FigS7_ComplexSims", "panel_s7c_raw_20231007b.pdf"), width = 120, height = 80, units = "mm")

# plot autocorrelation functions for S6C
data_autocor = list.files(here("01-rawdata"), "ct_dopc*", full.names = TRUE) %>% 
  read_csv(id = "fname") %>% 
  mutate(
    fname = basename(fname),
    press = parse_number(fname)
  )
  
# basic plots
data_autocor %>% 
  group_by(fname, press) %>% 
  arrange(Lag_time_ns) %>% 
  #head(50000) %>% 
  ggplot(
    aes(
      x = Lag_time_ns,
      ymin = Correlation + Confidence_interval_low_bound,
      y = Correlation,
      ymax = Correlation + Confidence_interval_high_bound
    )
  ) +
  facet_wrap(~press) +
  geom_ribbon(fill = "grey75") +
  geom_line() +
  theme_tiny()

# plots with "zoom"
lag_time_cutoff = 50 # time in ns at which to break the axis
data_autocor %>% 
  group_by(fname, press) %>% 
  mutate(bin = (row_number()-1) %/% 5) %>% 
  group_by(fname, press, bin) %>% 
  summarize(across(.cols = c(Lag_time_ns, Correlation, Confidence_interval_low_bound, Confidence_interval_high_bound), mean)) %>% 
  # downsample to make the vector manageable
  cross_join(tibble(part = c("beg", "end"))) %>% 
  filter(
    ((part == "beg") & (Lag_time_ns <= lag_time_cutoff)) |
      ((part == "end") & (Lag_time_ns > lag_time_cutoff))
  ) %>% 
  mutate(press_part = paste(press, part)) %>% 
  ggplot(
    aes(
      x = Lag_time_ns,
      ymin = Correlation + Confidence_interval_low_bound,
      y = Correlation,
      ymax = Correlation + Confidence_interval_high_bound
    )
  ) +
  facet_wrap(~press_part, nrow = 1, scales = "free_x") +
  geom_ribbon(fill = "grey75") +
  geom_line() +
  theme_tiny() +
  labs(
    x = "Lag time (nsec)",
    y = "Autocorrelation"
  )
ggsave(here("04-suppfigs", "Panel_s6c_autocorrels_20231012b.pdf"), width = 180, height = 40, units = "mm")
