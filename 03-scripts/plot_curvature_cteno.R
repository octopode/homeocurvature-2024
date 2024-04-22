## Plot lipidomes' mean curvature under multiple estimation schemes

library(tidyverse)
library(here)
library(ggpubr)

source(here("03-scripts", "prep_curvature.R"))
source(here("03-scripts", "plot_helpers.R"))

# in desired plotting order
sp_subset = c(
  "Boli_vitr", 
  "Leuc_pulc",
  "Boli_infu", 
  "Boli_micr", 
  "Lamp_crue", 
  "Tjal_pink"
)

# plot some species means
plot_plci_means = plcidata_wildwhole %>% 
  filter(sp %in% sp_subset) %>% 
  mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  # id-wise averaging relies on explicit zeroes!!
  #group_by(scheme, sp, class, id) %>% 
  #summarise(plci = mean(plci)) %>% 
  # consolidate headgroup blocks
  group_by(scheme, sp, eid, class) %>%
  summarize(plci = sum(plci)) %>% 
  # average those blocks
  group_by(scheme, sp, class) %>%
  summarize(
    plci = mean(plci),
    n = n()
  ) %>% 
  mutate(
    sp_n = {str_glue("{sp} N = {max(n)}")},
    sp_n = sp_n %>% factor(levels = unique(sp_n))
  ) %>% 
  gg_plcurv(
    rows = scheme,
    cols = sp_n,
    y = plci
  ) +
  labs(
    title = "Species mean c0 predicted under 3 different schemes",
    x = "Species"
  )
plot_plci_means
#ggsave(here("04-mainfigs", "Fig3_curv_spmeans_orig5_3schemes_20230712.pdf"), width = 8, height = 10)

# same as above but linreg only, formatted as panel for Fig. S7
panel_s7c = plcidata_wildwhole %>% 
  filter(sp %in% sp_subset) %>% 
  filter(scheme == "linreg") %>% 
  mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  # id-wise averaging relies on explicit zeroes!!
  #group_by(scheme, sp, class, id) %>% 
  #summarise(plci = mean(plci)) %>% 
  # consolidate headgroup blocks
  group_by(scheme, sp, eid, class) %>%
  summarize(plci = sum(plci)) %>% 
  # average those blocks
  group_by(scheme, sp, class) %>%
  summarize(
    plci = mean(plci),
    n = n()
  ) %>% 
  mutate(
    sp_n = {str_glue("{sp} N = {max(n)}")},
    sp_n = sp_n %>% factor(levels = unique(sp_n))
  ) %>% 
  gg_plcurv(
    #rows = scheme,
    cols = sp_n,
    y = plci
  ) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x = unit(4, "mm"),
    axis.text.y = element_text(angle = 90, hjust = 0.5)
  ) +
  labs(
    x = "Species"
  )
panel_s7c
ggsave(here("04-suppfigs", "panel_s7c_20231228a.pdf"), width = 60, height = 60, unit = "mm")

# plot indls for all orig spp.
plot_plci_indls = plcidata_wildwhole %>% 
  filter(sp %in% sp_subset) %>% 
  mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  arrange(sp) %>% 
  mutate(
    sp_eid = paste(sp, eid),
    sp_eid = sp_eid %>% factor(levels = unique(sp_eid))
  ) %>% 
  gg_plcurv(
    rows = scheme,
    cols = sp_eid,
    y = plci
  ) +
  theme(strip.text.x = element_text(angle = 90)) +
  labs(
    title = "Individual c0 in 5 species, predicted under 3 different schemes",
    x = "Sample"
  )
#plot_plci_indls
#ggsave(here("04-omitfigs", "curv_indls_orig5_3schemes_20230712.pdf"), width = 14, height = 8)

panel_4b_old = plot_plci_means = plcidata_wildwhole %>% 
  # remove 0-curvature entries for legend's sake
  filter(plci != 0) %>% 
  # just 1 scheme
  filter(scheme == "linreg") %>% 
  filter(sp %in% sp_subset) %>% 
  mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  # id-wise averaging relies on explicit zeroes!!
  #group_by(scheme, sp, class, id) %>% 
  #summarise(plci = mean(plci)) %>% 
  # consolidate headgroup blocks
  group_by(scheme, sp, eid, class) %>%
  summarize(plci = sum(plci)) %>% 
  # average those blocks
  group_by(scheme, sp, class) %>%
  summarize(
    plci = mean(plci),
    n = n()
  ) %>% 
  mutate(
    sp_n = {str_glue("{sp}\nN = {max(n)}")},
    sp_n = sp_n %>% factor(levels = unique(sp_n))
  ) %>% 
  gg_plcurv(
    rows = scheme,
    cols = sp_n,
    y = plci
  ) +
  theme_tiny() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.placement = "outside",
    legend.position = c(0.45, 0.2), # coords from LL on whole plot area
    legend.background = element_blank()
  ) +
  guides(
    x = "none",
    fill = guide_legend(nrow = 2, keywidth = 0.5, keyheight = 0.5)
  ) +
  labs(
    x = "Species",
    y = "Mean c0 (Å^-1)",
    fill = "PL class"
  )
panel_4b_old
ggsave(here("04-mainfigs", "Fig4_homeocurvature", "panel_4b_20230712.pdf"), width = 80, height = 40, unit = "mm")

# Panel 4C legend with just HG curvatures
# c0_allschemes comes from c0_meta_analysis.R
data_4clegend = crossing(
  scheme = "linreg",
  class = c0_allschemes$class %>% unique(),
  carbsn1 = 18,
  dbonsn1 = TRUE,
  frac_molar = 1
) %>% 
  # join linreg scheme
  mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
  left_join(
    c0_allschemes %>% 
      filter(scheme == "linreg") %>% 
      select(scheme, class, carbsn1, dbonsn1, c0, tol) %>% 
      mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
      select(-dbonsn1),
    by = c("scheme", "class"), # first, join just by class
    suffix = c('', ".linreg"),
    relationship = "one-to-many"
  ) %>% 
  # thin redundant entries
  filter(
    is.na(carbsn1.linreg) |
      (
        (carbsn1 == carbsn1.linreg) &
          (unsatsn1 == unsatsn1.linreg)
      )
  ) %>% 
  select(-contains(".linreg")) %>% 
  # finally, calculate the actual curvature contributions and error
  # `plci` stands for PhosphoLipid Curvature Index
  mutate(
    plci = c0  * frac_molar,
    ctol = tol * frac_molar # should this scale linearly or by sqrt?
  ) %>% 
  # important!
  replace_na(list(plci = 0, ctol = 0)) %>% 
  arrange(c0) %>% 
  mutate(class = class %>% factor(., levels = unique(.)))

panel_4clegend = data_4clegend %>% 
  ggplot(
    aes(
      x = c0,
      y = class,
      fill = class
    )
  ) +
  geom_vline(xintercept = 0, size = 0.5/.pt) +
  geom_col() +
  scale_fill_manual(values = chroma_cl) +
  theme_tiny() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "c0, 18:1 sn-1 chain")
panel_4clegend
ggsave(here("04-mainfigs", "Fig4_homeocurvature", "panel_4c_20240405a.pdf"), width = 25, height = 40, units = "mm")

# The same legend, but salted
data_legend_salted = crossing(
  scheme = "salted",
  class = c0_allschemes$class %>% unique(),
  carbsn1 = 18,
  dbonsn1 = TRUE,
  frac_molar = 1
) %>% 
  # join salted scheme
  mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
  left_join(
    c0_allschemes %>% 
      filter(scheme == "salted") %>% 
      select(scheme, class, carbsn1, dbonsn1, c0, tol) %>% 
      mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
      select(-dbonsn1),
    by = c("scheme", "class"), # first, join just by class
    suffix = c('', ".salted"),
    relationship = "one-to-many"
  ) %>% 
  # thin redundant entries
  filter(
    is.na(carbsn1.salted) |
      (
        (carbsn1 == carbsn1.salted) &
          (unsatsn1 == unsatsn1.salted)
      )
  ) %>% 
  select(-contains(".salted")) %>% 
  # finally, calculate the actual curvature contributions and error
  # `plci` stands for PhosphoLipid Curvature Index
  mutate(
    plci = c0  * frac_molar,
    ctol = tol * frac_molar # should this scale linearly or by sqrt?
  ) %>% 
  # important!
  replace_na(list(plci = 0, ctol = 0)) %>% 
  arrange(c0) %>% 
  mutate(class = class %>% factor(., levels = unique(.)))

panel_legend_salted = data_legend_salted %>% 
  ggplot(
    aes(
      x = c0,
      y = class,
      fill = class
    )
  ) +
  geom_vline(xintercept = 0, size = 0.5/.pt) +
  geom_col() +
  scale_fill_manual(values = chroma_cl) +
  theme_tiny() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "c0, 18:1 sn-1 chain")
panel_legend_salted
ggsave(here("04-suppfigs", "FigS5_barc0", "legend_salted_20231012a.pdf"), width = 40, height = 60, units = "mm")

# simple boxplot of bar-c0 for Bath_fost and Boli_micr
data_panels8b = plcidata_indls %>% 
  filter(
    sp %in% c(
      "Boli_micr",
      "Bath_fost"
    )
  ) %>% 
  filter(scheme == "linreg") %>% 
  mutate(sp = factor(sp, levels = c("Boli_micr", "Bath_fost")))

panel_s8b = data_panels8b %>% 
  ggplot(
    aes(
      x = sp,
      y = plci,
      fill = sp
    )
  ) +
  geom_boxplot(
    color = "white",
    fill = "grey80",
    size = 0.5/.pt,
    width = 0.5
  ) +
  geom_point(
    size = 1,
    shape = 21,
    color = "white",
  ) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_fill_manual(values = chroma_sp) +
  labs(
    x = "Species",
    y = "Mean c0 (Å^-1)"
  )
panel_s8b
ggsave(here("04-suppfigs", "FigS8_ConfocalLaurdan", "panel_S8b_20231005a.pdf"), width = 35, height = 50, unit = "mm")

# summary stats for text
data_panels8b %>% 
  group_by(sp) %>% 
  summarize(
    plci_sem = sd(plci)/sqrt(n()),
    plci_tol = max(plci) - mean(plci),
    plci_avg = mean(plci)
  )

# t-test!
data_panels8b %>% 
  ungroup() %>% 
  select(sp, plci) %>% 
  pivot_wider(names_from = "sp", values_from = "plci") %>%
  unnest(Boli_micr, Bath_fost) %>% 
  summarise(test = t.test(Boli_micr, Bath_fost) %>% tidy())

## SCRATCH ##

## scatterplot indls everyone
#bind_rows(
#  plcidata_wildwhole %>% 
#    mutate(subset = "all samples"),
#  plcidata_wildwhole %>% 
#    filter(temp_col <= 10) %>% 
#    mutate(subset = "<10 deg C")
#) %>% 
#  # QC: need to implement this at an earlier stage!
#  filter(eid != "JWL0169") %>% 
#  # get total c0 for each indl
#  group_by(scheme, sp, eid, depth_col, temp_col, subset) %>%
#  summarize(plci = sum(plci)) %>% 
#  ggplot(
#    aes(
#      x = depth_col,
#      y = plci
#    )
#  ) +
#  facet_grid(rows = vars(scheme), cols = vars(subset)) +
#  geom_point(
#    aes(
#      fill = sp,
#      shape = temp_col <= 10
#    ),
#    size = 2,
#    color = "white"
#  ) +
#  geom_smooth(method = "lm", color = "black") +
#  #geom_text(aes(label = eid)) +
#  scale_fill_manual(values = chroma_sp) +
#  scale_shape_manual(values = c(22, 21)) +
#  guides(shape = "none") +
#  theme_pubr()
#ggsave(here("04-omitfigs", "curv-vs-depth_indls_all_3schemes_20230712.pdf"), width = 6, height = 8)
