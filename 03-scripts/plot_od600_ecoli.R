## Plot OD600 data for engineered E. coli strains under pressure

library(tidyverse)
library(here)
library(httr)

source(here("03-scripts", "plot_helpers.R"))

file_od600 = here("01-rawdata", "od600_results.csv")
# this is a Google Sheets key for the curvature info file, to facilitate a nice Excel-like interface
gsht_od600 = "1mqVLOlpHpnxqvT25_MlLJPSyM7Wp8A2d14m9ehz0yTM"

# refresh c0 spreadsheet from online interface
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_od600}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_od600)

col_width = 150
jit_width = 50

# load data and calc changes
data_od600 = file_od600 %>% 
  read_csv() %>% 
  mutate(
    delta_od = od_end - od_beg,
    doubling = log(od_end / od_beg, 2),
    strain   = factor(strain, levels = c("PE", "PPE", "PC", "PG+", "PG-")),
    expt     = factor(expt, levels = c("PPE", "PC", "PG"))
  )

# get means + stderr
data_od600_summ = data_od600 %>% 
  group_by(expt, strain, press) %>% 
  summarize(
    avg_delta_od = mean(delta_od),
    ser_delta_od = sd(delta_od)/sqrt(n()),
    avg_doubling = mean(doubling),
    ser_doubling = sd(doubling)/sqrt(n())
  )

# plot growth
data_od600_summ %>% 
  ggplot(aes(x = press, y = avg_delta_od, fill = strain)) +
  facet_grid(cols = vars(expt)) +
  geom_col(
    position = position_dodge(),
    width = col_width
  ) +
  geom_point(
    data = data_od600,
    aes(y = delta_od),
    position = position_jitterdodge(
      jitter.width = jit_width,
      dodge.width = col_width
    ),
    color = "black"#,
    #shape = 'o'
  ) +
  geom_errorbar(
    aes(
      group = strain,
      ymin = avg_delta_od - ser_delta_od,
      ymax = avg_delta_od + ser_delta_od
    ),
    position = position_dodge(width = col_width),
    width = 25
  ) +
  scale_x_continuous(breaks = c(0, 250, 500)) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(
    title = "Delta OD",
    x = "Pressure (bar) for 24 h @ 37°C",
    y = "Change in OD600"
  )
 #ggsave(here("04-omitfigs", "delta_od600_20230802.pdf"), width = 8, height = 4)
 
 # plot doublings
 data_od600_summ %>% 
   ggplot(aes(x = press, y = avg_doubling, fill = strain)) +
   facet_grid(cols = vars(expt)) +
   geom_col(
     position = position_dodge(),
     width = col_width
   ) +
   geom_point(
     data = data_od600,
     aes(y = doubling),
     position = position_jitterdodge(
       jitter.width = jit_width,
       dodge.width = col_width
     ),
     color = "black"#,
     #shape = 'o'
   ) +
   geom_errorbar(
     aes(
       group = strain,
       ymin = avg_doubling - ser_doubling,
       ymax = avg_doubling + ser_doubling
     ),
     position = position_dodge(width = col_width),
     width = 25
   ) +
   scale_x_continuous(breaks = c(0, 250, 500)) +
   scale_fill_manual(values = chroma_cl) +
   theme_pubr() +
   theme(legend.position = "none") +
   labs(
     title = "Doublings",
     x = "Pressure (bar) for 24 h @ 37°C",
     y = "Change in OD600"
   )
#ggsave(here("04-omitfigs", "doublings_od600_20230802.pdf"), width = 8, height = 4)

# NORMALIZE stuff to mean of that strain at 0 bar.
data_od600_norm = data_od600 %>% 
  left_join(data_od600_summ, by = c("press", "strain", "expt")) %>% 
  group_by(strain, expt) %>% 
  arrange(press, rep) %>% 
  mutate(
    fold_dod = delta_od / avg_delta_od[[1]],
    fold_dou = doubling / avg_doubling[[1]]
  )

# calc means + stderr
data_od600_norm_summ = data_od600_norm %>% 
  # slope difference tests
  group_by(expt) %>% 
  mutate(
    lin_mod = lm(fold_dod ~ strain*press, cur_data()) %>% 
      list()
  ) %>% 
  ## pairwise t-tests
  #group_by(expt, press) %>% 
  #mutate(
  #  pval_dod = t.test(
  #    x = cur_data() %>% filter(strain != "PE") %>% .$fold_dod,
  #    y = cur_data() %>% filter(strain == "PE") %>% .$fold_dod,
  #    alternative = 't'
  #  ) %>% .$p.value,
  #  pval_dou = t.test(
  #    x = cur_data() %>% filter(strain != "PE") %>% .$fold_dou,
  #    y = cur_data() %>% filter(strain == "PE") %>% .$fold_dou,
  #    alternative = 't'
  #  ) %>% .$p.value
  #) %>% 
  #group_by(expt, strain, press, pval_dod, pval_dou) %>% 
group_by(expt, strain, press) %>% 
  summarize(
    avg_fold_dod = mean(fold_dod),
    ser_fold_dod = sd(fold_dod)/sqrt(n()),
    avg_fold_dou = mean(fold_dou),
    ser_fold_dou = sd(fold_dou)/sqrt(n()),
    lin_mod = lin_mod
  )

# get slope P-vals
# these may be important!
data_od600_norm_summ %>%
  group_by(expt) %>% 
  summarize(lin_mod = lin_mod[1]) %>% 
  .$lin_mod %>% 
  lapply(., summary)

## plot norm'd growth
#data_od600_norm_summ %>% 
#  ggplot(aes(x = press, y = avg_fold_dod, fill = strain)) +
#  facet_grid(cols = vars(expt)) +
#  geom_col(
#    position = position_dodge(),
#    width = col_width
#  ) +
#  geom_point(
#    data = data_od600_norm,
#    aes(y = fold_dod),
#    position = position_jitterdodge(
#      jitter.width = jit_width,
#      dodge.width = col_width
#    ),
#    color = "black"#,
#    #shape = 'o'
#  ) +
#  geom_errorbar(
#    aes(
#      group = strain,
#      ymin = avg_fold_dod - ser_fold_dod,
#      ymax = avg_fold_dod + ser_fold_dod
#    ),
#    position = position_dodge(width = col_width),
#    width = 25
#  ) +
#  geom_text(
#    aes(label = round(pval_dou, 3)),
#    y = 1.2,
#    check_overlap = TRUE
#  ) +
#  scale_x_continuous(breaks = c(0, 250, 500)) +
#  scale_fill_manual(values = chroma_cl) +
#  theme_pubr() +
#  theme(legend.position = "none") +
#  labs(
#    title = "Fold change: delta OD",
#    x = "Pressure (bar) for 24 h @ 37°C",
#    y = "Fold change in delta OD600"
#  )
#ggsave(here("04-omitfigs", "foldchg_od600_20230802.pdf"), width = 8, height = 4)

## plot norm'd growth
#data_od600_norm_summ %>% 
#  ggplot(aes(x = press, y = avg_fold_dou, fill = strain)) +
#  facet_grid(cols = vars(expt)) +
#  geom_col(
#    position = position_dodge(),
#    width = col_width
#  ) +
#  geom_point(
#    data = data_od600_norm,
#    aes(y = fold_dou),
#    position = position_jitterdodge(
#      jitter.width = jit_width,
#      dodge.width = col_width
#    ),
#    color = "black"#,
#    #shape = 'o'
#  ) +
#  geom_errorbar(
#    aes(
#      group = strain,
#      ymin = avg_fold_dou - ser_fold_dou,
#      ymax = avg_fold_dou + ser_fold_dou
#    ),
#    position = position_dodge(width = col_width),
#    width = 25
#  ) +
#  scale_x_continuous(breaks = c(0, 250, 500)) +
#  scale_fill_manual(values = chroma_cl) +
#  theme_pubr() +
#  theme(legend.position = "none") +
#  labs(
#    title = "Fold change: doublings",
#    x = "Pressure (bar) for 24 h @ 37°C",
#    y = "Fold change in doublings"
#  )
#ggsave(here("04-omitfigs", "foldchg_doublings_od600_20230802.pdf"), width = 8, height = 4)

# Plot for Fig 5B
panel_5b = data_od600_norm_summ %>% 
  # filter out the HDL11 PG expt
  filter(expt %in% c("PPE", "PC")) %>% 
  ggplot(aes(x = press, y = avg_fold_dod, fill = strain)) +
  facet_grid(cols = vars(expt)) +
  geom_col(
    position = position_dodge(),
    width = col_width,
    color = "white",
    linewidth = 1/.pt
  ) +
  geom_errorbar(
    aes(
      group = strain,
      ymin = avg_fold_dod - ser_fold_dod,
      ymax = avg_fold_dod + ser_fold_dod
    ),
    position = position_dodge(width = col_width),
    width = 25,
    size = 0.5/.pt,
    color = "grey25"
  ) +
  geom_point(
    data = data_od600_norm %>% filter(expt %in% c("PPE", "PC")),
    aes(y = fold_dod),
    position = position_jitterdodge(
      jitter.width = jit_width,
      dodge.width = col_width
    ),
    color = "black",
    size = 2/.pt
  ) +
  #geom_text(
  #  aes(label = round(pval_dou, 3)),
  #  y = 1.2,
  #  check_overlap = TRUE,
  #  size = 7/.pt
  #) +
  scale_x_continuous(breaks = c(0, 250, 500)) +
  scale_fill_manual(values = chroma_cl) +
  theme_tiny() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  labs(
    x = "Pressure (bar)",
    y = "ddOD600 (fold diff.)"
  )
panel_5b
ggsave(here("04-mainfigs", "Fig5_ecolimodel", "foldchg_od600_20230802a.pdf"), width = 60, height = 40, units = "mm")

# Plot for Fig S9E
panel_s9e = data_od600_norm_summ %>% 
  # filter out the HDL11 PG expt
  filter(expt %in% c("PG")) %>% 
  ggplot(aes(x = press, y = avg_fold_dod, fill = strain)) +
  facet_grid(cols = vars(expt)) +
  geom_col(
    color = "black",
    position = position_dodge(),
    width = col_width,
    color = "white",
    linewidth = 1/.pt
  ) +
  geom_errorbar(
    aes(
      group = strain,
      ymin = avg_fold_dod - ser_fold_dod,
      ymax = avg_fold_dod + ser_fold_dod
    ),
    position = position_dodge(width = col_width),
    width = 25,
    size = 0.5/.pt,
    color = "grey25"
  ) +
  geom_point(
    data = data_od600_norm %>% filter(expt %in% c("PG")),
    aes(y = fold_dod),
    position = position_jitterdodge(
      jitter.width = jit_width,
      dodge.width = col_width
    ),
    color = "black",
    size = 2/.pt
  ) +
  #geom_text(
  #  aes(label = round(pval_dou, 3)),
  #  y = 1.2,
  #  check_overlap = TRUE,
  #  size = 7/.pt
  #) +
  scale_x_continuous(breaks = c(0, 250, 500)) +
  scale_fill_manual(values = c("PG-" = "white", "PG+" = chroma_cl[["PG"]])) +
  theme_tiny() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  labs(
    x = "Pressure (bar)",
    y = "ddOD600 (fold diff.)"
  )
panel_s9e
ggsave(here("04-suppfigs", "FigS9_Ecoli", "panel_s9e_HDL11_20231008a.pdf"), width = 60, height = 40, units = "mm")
