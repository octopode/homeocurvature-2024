## Scatter/bubble-plot lipidomes' mean curvature under multiple estimation schemes
## Barplots and lipid contributions to curvature are in `plot_curvature.R`

library(tidyverse)
library(here)
library(ggpubr)
library(ggtree)

source(here("03-scripts", "prep_pgls_lipidmaps.R"))
source(here("03-scripts", "plot_helpers.R"))
source(here("03-scripts", "pgls_helpers.R"))

# Panel 3a: Depth-temp distribution plots
# Bubble and regular scatter with jitter
panel_3a_bubble = pldata_indls %>% 
  # condense lipid classes
  group_by(sp, eid, depth_col, temp_col) %>% 
  summarise() %>% 
  # order warm->cold, shal->deep?
  # condense species for bubble
  group_by(sp) %>% 
  summarize(
    n = n(),
    depth_col_serr = sd(depth_col)/sqrt(n),
    depth_col = mean(depth_col),
    temp_col_serr = sd(temp_col)/n,
    temp_col = mean(temp_col)
  ) %>% 
  arrange(-n) %>% 
  ggplot(
    aes(
      x = temp_col,
      y = depth_col
    )
  ) +
  # lines for subsets
  geom_vline(xintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 250, linetype = "dotted") +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      ymin = depth_col - depth_col_serr,
      ymax = depth_col + depth_col_serr
    )
  ) +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      xmin = temp_col - temp_col_serr,
      xmax = temp_col + temp_col_serr
    )
  ) +
  geom_point(
    shape = 21,
    color = "white",
    aes(
      size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = c(1, 3), breaks = seq(1,7,2)) +
  scale_fill_manual(values = chroma_sp) +
  scale_x_continuous(position = "top") +
  scale_y_reverse() +
  theme_tiny() +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = c(0.8, 0.375),
    #legend.position = "none",
    #legend.margin = margin(rep(c(0,30), 2)), # 10 mm "breathing room" for legend
    legend.margin = margin(rep(0, 4)),
    legend.background = element_blank(),
    legend.title.align = 0.5,
    legend.spacing.y = unit(0.5, "mm"),
  ) +
  guides(
    size = guide_legend(
      override.aes = list(fill = "black"), 
      keyheight = 0.5,
      keywidth = 0.5,
      nrow = 1
    ),
    fill = guide_legend(
      keyheight = 0.5,
      keywidth = 0.5,
      ncol = 2
    ),
    linetype = "none"
  ) +
  labs(
    x = "Collection temperature (deg C)",
    y = "Collection depth (m)",
    size = "N",
    fill = "Species"
  )
panel_3a_bubble
ggsave(here("04-mainfigs", "Fig3_ctenolipidomics", "depth_v_temp_bubble_20230919a.pdf"), width = 50, height = 45, units = "mm")

panel_3a_scatter = pldata_indls %>% 
  # condense lipid classes
  group_by(sp, eid, depth_col, temp_col) %>% 
  summarise() %>% 
  # order warm->cold, shal->deep?
  ggplot(
    aes(
      x = temp_col,
      y = depth_col
    )
  ) +
  # lines for subsets
  geom_vline(xintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 250, linetype = "dotted") +
  geom_point(
    shape = 21,
    color = "white",
    aes(fill = sp)
  ) +
  scale_fill_manual(values = chroma_sp) +
  scale_y_reverse() +
  theme_tiny() +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = c(0.6, 0.425),
    legend.margin = margin(rep(c(0,4), 2)),
    legend.background = element_blank(),
    legend.title.align = 0.5,
    legend.spacing.y = unit(0.5, "mm"),
  ) +
  guides(
    size = guide_legend(
      override.aes = list(fill = "black"), 
      keyheight = 0.5,
      keywidth = 0.5,
      nrow = 1
    ),
    fill = guide_legend(
      keyheight = 0.5,
      keywidth = 0.5,
      ncol = 2
    ),
    linetype = "none"
  ) +
  labs(
    x = "Collection temperature (deg C)",
    y = "Collection depth (m)",
    fill = "Species"
  )
panel_3a_scatter
ggsave(here("04-mainfigs", "Fig3_ctenolipidomics", "panel_3a_20230919a.pdf"), width = 50, height = 45, units = "mm")

# Bubbleplot with OLS and GLS for each lipid class
# modifying for panel S2A: each class vs. depth and temp
panel_s2a = pldata_indls_opt %>% 
  filter(
    ((subset %in% c("all samples", "<= 250 m")) & (predvar == "temp_col")) |
      ((subset %in% c("all samples", "<= 10 deg C")) & (predvar == "depth_col"))
  ) %>% 
  # since this is a class-totals analysis
  filter(class != "all") %>% 
  # non-jackknifed dataset
  filter(jknife == "nothing") %>% 
  # and just class fractions
  filter(respvar == "frac_molar") %>% 
  # calc means and std errors
  # need merged faceting var
  #mutate(predvar_subset = paste(predvar, subset)) %>% 
  mutate(subset_predvar = paste(subset, predvar)) %>% 
  group_by(class, sp, subset, predvar, subset_predvar) %>% 
  summarize(
    n = n(),
    resp_serr = sd(resp)/sqrt(n),
    resp = mean(resp),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  arrange(-n) %>% # point stacking order
  ## view sp means
  #filter(
  #  (sp %in% c("Bath_fost", "Boli_micr")) &
  #    (class == "PPE") &
  #    (subset == "all samples")
  #)
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  facet_grid(rows = vars(class), cols = vars(subset_predvar), scales = "free") +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    size = 7/.pt,
    #x = 2000,
    y = 0,
    data = allmods_multtest %>% 
      # since this is a class-totals analysis
      filter(class != "all") %>% 
      # non-jackknifed dataset
      filter(jknife == "nothing") %>% 
      # and just class fractions
      filter(respvar == "frac_molar") %>% 
      filter(term == "pred") %>% 
      # need merged faceting var
      #mutate(predvar_subset = paste(predvar, subset)) %>% 
      mutate(subset_predvar = paste(subset, predvar)) %>% 
      group_by(class, predvar, subset, subset_predvar) %>% 
      arrange(model) %>% 
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')) %>% 
      # special var to place the p-vals
      mutate(labx = ifelse(predvar == "depth_col", 2000, 10)),
    aes(
      x = labx,
      label = plabel
    ),
    hjust = 0.5,
    vjust = -1.5,
    check_overlap = TRUE
  ) +
  # best fit lines
  geom_line(
    color = "black",
    size = 0.5*0.75,
    data = modelines_lmap %>% 
      # since this is a class-totals analysis
      filter(class != "all") %>% 
      # non-jackknifed dataset
      filter(jknife == "nothing") %>% 
      # and just class fractions
      filter(respvar == "frac_molar") %>% 
      # need merged faceting var
      #mutate(predvar_subset = paste(predvar, subset)) %>% 
      mutate(subset_predvar = paste(subset, predvar)) %>% 
      ungroup(),
    aes(linetype = model)
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      ymin = resp - resp_serr,
      ymax = resp + resp_serr
    )
  ) +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      xmin = pred - pred_serr,
      xmax = pred + pred_serr
    )
  ) +
  geom_point(
    shape = 21,
    color = "white",
    aes(
      size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = 0.75*c(1, 3)) +
  scale_fill_manual(values = chroma_sp) +
  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    panel.spacing.x = unit(4, "mm"),
    panel.spacing.y = unit(2, "mm")
  ) +
  guides(fill = "none") +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    #title = "Class fractions vs. depth by dataset\nBH-corrected P-vals",
    x = "Predictor value (m or deg C)",
    y = "Mole fraction",
    size = "N"
  )
panel_s2a
ggsave(here("04-suppfigs", "FigS2_acylchains", "panel_s2a_20231012c.pdf"), width = 88, height = 200, units = "mm")

# Just PPE <= 10°C + PPE <= 200 m for fig panel
filter_panel3cd = function(x){
  x %>% 
    # non-jackknifed dataset
    filter(jknife == "nothing") %>% 
    ## temporary til we get more Leucos
    #filter(
    #  ((predvar == "depth_col") & (jknife == "nothing")) |
    #    ((predvar == "temp_col") & (jknife == "Leuc_pulc"))
    #) %>% 
    # and just class fractions
    filter(respvar == "frac_molar") %>% 
    filter(
      ((subset == "<= 10 deg C") & (class == "PPE") & (predvar == "depth_col")) |
        ((subset == "<= 250 m") & (class == "PPE") & (predvar == "temp_col"))
    )
}

# generate plot
panel_3cd = pldata_indls_opt %>% 
  filter_panel3cd() %>% 
  # calc means and std errors
  group_by(class, sp, subset, predvar) %>% 
  summarize(
    n = n(),
    resp_serr = sd(resp)/sqrt(n),
    resp = mean(resp),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  arrange(-n) %>% 
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  facet_wrap(~predvar, scales = "free_x", strip.position = "bottom", ncol = 2) +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    size = 7/.pt,
    x = c(1000, 5),
    y = c(0.4, 0.4),
    data = allmods_multtest %>% 
      filter_panel3cd() %>%
      filter(term == "pred") %>% 
      group_by(class, subset, predvar) %>% 
      arrange(model) %>% 
      mutate(plabel = str_glue({"{model} P = {format(p.value_BH, digits=3, scientific=TRUE)}"})) %>% 
      summarize(
        plabel = paste(plabel, collapse='\n')
      ),
    aes(
      label = plabel,
      x = pred
    ),
    hjust = 0.5,
    vjust = 0,
    check_overlap = TRUE
  ) +
  # best fit lines
  geom_line(
    size = 1/.pt,
    color = "black",
    data = modelines_lmap %>% 
      filter_panel3cd() %>%
      ungroup(),
    aes(linetype = model)
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      ymin = resp - resp_serr,
      ymax = resp + resp_serr
    )
  ) +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      xmin = pred - pred_serr,
      xmax = pred + pred_serr
    )
  ) +
  geom_point(
    shape = 21,
    color = "white",
    aes(
      size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = c(1, 3), breaks = seq(1,7,2)) +
  scale_fill_manual(values = chroma_sp) +
  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
  theme_tiny() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    #legend.position = "left",
    legend.position = "none",
    legend.margin = margin(rep(c(0,30), 2)), # 10 mm "breathing room" for legend
    legend.background = element_blank(),
    legend.title.align = 0.5,
    legend.spacing.y = unit(0.5, "mm"),
  ) +
  guides(
    fill = guide_legend(
      keyheight = 0.5,
      keywidth = 0.5,
      ncol = 2
    ),
    size = guide_legend(
      override.aes = list(fill = "black"), 
      keyheight = 0.5,
      keywidth = 0.5,
      nrow = 1
    ),
    linetype = "none"
  ) +
  lims(y = c(0, NA)) +
  labs(
    y = "Mole fraction",
    fill = "Species",
    size = "N"
  )
panel_3cd
ggsave(here("04-mainfigs", "Fig3_ctenolipidomics", "panel_3cd_20230928b.pdf"), width = 70, height = 40, units = "mm")

## Claims about species jackknife analysis with just PPE
jknife_pvals_lmap = allmods_multtest %>% 
  filter(subset == "<= 10 deg C") %>% 
  filter(class == "PPE") %>% 
  filter((term == "pred") & (predvar == "depth_col")) %>% 
  filter(model == "modgls") %>% 
  arrange(-p.value, model, class, jknife)

jknife_pvals_lmap
jknife_pvals_lmap$p.value %>% max()
# _This_ signal in particular counts on Tjalfie.

## Chain analyses!

filter_panel3efgh = function(x){
  x %>% 
    # since this is a class-totals analysis
    filter(class == "all") %>% 
    # non-jackknifed dataset
    filter(jknife == "nothing") %>% 
    # select readouts
    #filter(respvar %in% c("carbsn1", "dbonsn1")) %>% # just sn-1 20231002
    filter(respvar %in% c("chn", "dbi")) %>% 
    filter(
      ((predvar == "depth_col") & (subset == "<= 10 deg C")) |
        ((predvar == "temp_col") & (subset == "<= 250 m"))
    ) %>% 
    # add the faceting var why not?
    mutate(trend = paste(respvar, predvar, sep='~')) %>% 
    arrange(predvar, respvar) %>% 
    mutate(trend = trend %>% factor(levels = unique(.)))
}

# The [optional] acyl chains panel for Fig. 3
panel_3efgh = pldata_indls_opt %>% 
  filter_panel3efgh() %>% 
  # calc means and std errors
  group_by(respvar, sp, subset, predvar, trend) %>% 
  summarize(
    n = n(),
    resp_serr = sd(resp)/sqrt(n),
    resp = mean(resp),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  #facet_wrap(~trend, scales = "free", strip.position = "left", nrow = 1) +
  facet_grid(rows = vars(respvar), cols = vars(predvar), scales = "free", switch = "both") +
  ## introduce some invisible data to make axes uniform
  #geom_point(
  #  color = "transparent",
  #  data = tibble(
  #    trend = rep(c("chn~depth_col", "dbi~depth_col", "chn~temp_col", "dbi~temp_col"), 2) %>% 
  #      factor(levels = unique(.)),
  #    # max, min
  #    pred = c(c(4000, 4000, 28, 28), c(0, 0, -3.5, -3.5)),
  #    resp = c(c(20, 3.25, 20, 3.25), c(18, 1.4, 18, 1.4))
  #  )
  #) +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    size = 7/.pt,
    # no idea how this ordering works; trial and error
    x = c(1500, 7.5, 1500, 7.5),
    y = c(18.5, 18.5, 2, 2),
    data = allmods_multtest %>% 
      filter_panel3efgh() %>% 
      # OLS y-intercept for locating text
      group_by(subset, predvar, respvar, trend, jknife) %>% 
      mutate(yint = estimate[[1]]) %>% 
      filter(term == "pred") %>% 
      group_by(respvar, predvar, subset, trend, yint) %>% 
      arrange(model) %>% 
      # I don't think mult testing corr is appropriate here...are CL and DBI really a family?
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')),
    aes(label = plabel),
    hjust = 0.5,
    vjust = 0.5,
    check_overlap = TRUE
  ) +
  # best fit lines
  geom_line(
    color = "black",
    size = 1/.pt,
    data = modelines_lmap %>% 
      filter_panel3efgh() %>% 
      ungroup(),
    aes(
      linetype = model
    )
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      ymin = resp - resp_serr,
      ymax = resp + resp_serr
    )
  ) +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      xmin = pred - pred_serr,
      xmax = pred + pred_serr
    )
  ) +
  geom_point(
    shape = 21,
    color = "white",
    aes(
      size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = c(1, 3)) +
  #scale_y_continuous(breaks = c(seq(28, 40, 4), seq(3, 6))) +
  scale_fill_manual(values = chroma_sp) +
  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
  theme_tiny() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    # replacing the y label with the strip text
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none"
  ) +
  guides(
    linetype = "none"
  )
panel_3efgh
ggsave(here("04-mainfigs", "Fig3_ctenolipidomics", "panel_efgh_20231006a.pdf"), width = 60, height = 60, units = "mm")

## panel S2B: stereospecific acyl chains
filter_panels2b = function(x){
  x %>% 
    # only interested in aggregate chain trends
    filter(class == "all") %>% 
    # non-jackknifed dataset
    filter(jknife == "nothing") %>% 
    # and just sn-1, 2
    filter(str_detect(respvar, regex("sn[12]"))) %>% 
    filter(
      ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
        ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
    ) %>% 
    # add the faceting vars
    mutate(
      stereo = parse_number(respvar),
      respvar = str_remove(respvar, regex("[0-9]")),
      stereo_respvar = str_glue("sn-{stereo} {respvar}"),
      subset_predvar = paste(subset, predvar)
    ) %>% 
    ungroup() %>% 
    arrange(predvar, respvar, subset, stereo)
}

## Panel S2B
panel_s2b = pldata_indls_opt %>% 
  filter_panels2b() %>% #.$trend %>% levels()
  # calc means and std errors
  group_by(respvar, sp, subset, predvar, subset_predvar, stereo_respvar) %>% 
  summarize(
    n = n(),
    resp_serr = sd(resp, na.rm = TRUE)/sqrt(n),
    resp = mean(resp, na.rm = TRUE),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  facet_grid(rows = vars(stereo_respvar), cols = vars(subset_predvar), scales = "free") +
  #facet_wrap(~trend, scales = "free", strip.position = "left", ncol = 4) +
  # invisible points to stretch the axes
  geom_point(
    data = pldata_indls_opt %>% 
      filter_panels2b() %>% 
      group_by(respvar, stereo) %>% 
      mutate(resp = range(resp, na.rm = TRUE) %>% list()) %>% 
      unnest(resp),
    x = 0, # since we're really only dealing with vert ranges
    color = "transparent",
    fill  = "transparent"
  ) +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    size = 7/.pt,
    # no idea how this ordering works; trial and error
    x = rep(c(rep(4000, 4), rep(25, 4)), 2),
    y = c(rep(c(17.0, 20.4), 4), rep(c(0.4, 4.6), 4)),
    data = allmods_multtest %>% 
      filter_panels2b() %>% 
      # OLS y-intercept for locating text
      group_by(subset, predvar, respvar, subset_predvar, stereo_respvar, jknife) %>% 
      mutate(yint = estimate[[1]]) %>% 
      filter(term == "pred") %>% 
      group_by(respvar, predvar, subset, subset_predvar, stereo_respvar, yint) %>% 
      arrange(model) %>% 
      # I don't think mult testing corr is appropriate here...are CL and DBI really a family?
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')),
    aes(label = plabel),
    hjust = 1.0,
    vjust = 0.5,
    check_overlap = TRUE
  ) +
  # best fit lines
  geom_line(
    color = "black",
    size = 0.75*0.5,
    data = modelines_lmap %>% 
      filter_panels2b() %>% 
      ungroup(),
    aes(
      linetype = model
    )
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      ymin = resp - resp_serr,
      ymax = resp + resp_serr
    )
  ) +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey25",
    aes(
      xmin = pred - pred_serr,
      xmax = pred + pred_serr
    )
  ) +
  geom_point(
    shape = 21,
    color = "white",
    aes(
      size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = 0.75*c(1, 3)) +
  scale_fill_manual(values = chroma_sp) +
  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
  theme_tiny() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    # replacing the y label with the strip text
    strip.background = element_blank(),
    strip.text = element_text(vjust = -0.5),
    strip.placement = "outside",
    legend.position = "none"#,
    #panel.spacing = unit(-1, "mm")
  ) +
  guides(
    linetype = "none"
  ) +
  labs(
    x = "Predictor value",
    y = "Response value"
  )
panel_s2b
ggsave(here("04-suppfigs", "FigS2_acylchains", "panel_s2b_20231006a.pdf"), width = 88, height = 90, units = "mm")

## WHITE-ON-BLACK PLOTS FOR SLIDES
# implemented using code for panel S2A with a filter
filter_headgp_pgls = function(datin){
  datin %>% 
    filter(class == "PPE") %>% 
    #filter(
    #  ((predvar == "depth_col") & (subset == "<= 10 deg C")) |
    #    ((predvar == "temp_col") & (subset == "<= 250 m"))
    #  )
    #filter(subset == "<= 10 deg C")
    filter(subset == "all samples")
}

slide_headgp_pgls = pldata_indls_opt %>% 
  filter_headgp_pgls() %>% 
  filter(
    ((subset %in% c("all samples", "<= 250 m")) & (predvar == "temp_col")) |
      ((subset %in% c("all samples", "<= 10 deg C")) & (predvar == "depth_col"))
  ) %>% 
  # since this is a class-totals analysis
  filter(class != "all") %>% 
  # non-jackknifed dataset
  filter(jknife == "nothing") %>% 
  # and just class fractions
  filter(respvar == "frac_molar") %>% 
  # calc means and std errors
  # need merged faceting var
  mutate(predvar_subset = paste(predvar, subset)) %>% 
  group_by(class, sp, subset, predvar, predvar_subset) %>% 
  summarize(
    n = n(),
    resp_serr = sd(resp)/sqrt(n),
    resp = mean(resp),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  arrange(-n) %>% # point stacking order
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  facet_grid(rows = vars(class), cols = vars(predvar_subset), scales = "free") +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    color = "white",
    size = 7/.pt,
    #x = 2000,
    y = 0,
    data = allmods_multtest %>% 
      filter_headgp_pgls() %>% 
      # since this is a class-totals analysis
      filter(class != "all") %>% 
      # non-jackknifed dataset
      filter(jknife == "nothing") %>% 
      # and just class fractions
      filter(respvar == "frac_molar") %>% 
      filter(term == "pred") %>% 
      # need merged faceting var
      mutate(predvar_subset = paste(predvar, subset)) %>% 
      group_by(class, predvar, subset, predvar_subset) %>% 
      arrange(model) %>% 
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')) %>% 
      # special var to place the p-vals
      mutate(labx = ifelse(predvar == "depth_col", 2000, 10)),
    aes(
      x = labx,
      label = plabel
    ),
    hjust = 0.5,
    vjust = -1.5,
    check_overlap = TRUE
  ) +
  # best fit lines
  geom_line(
    color = "white",
    data = modelines_lmap %>% 
      filter_headgp_pgls() %>% 
      # since this is a class-totals analysis
      filter(class != "all") %>% 
      # non-jackknifed dataset
      filter(jknife == "nothing") %>% 
      # and just class fractions
      filter(respvar == "frac_molar") %>% 
      # need merged faceting var
      mutate(predvar_subset = paste(predvar, subset)) %>% 
      ungroup(),
    aes(linetype = model)
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey75",
    aes(
      ymin = resp - resp_serr,
      ymax = resp + resp_serr
    )
  ) +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey75",
    aes(
      xmin = pred - pred_serr,
      xmax = pred + pred_serr
    )
  ) +
  geom_point(
    shape = 21,
    color = "white",
    aes(
      size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = c(2, 6), breaks = seq(1,7,2)) +
  scale_fill_manual(values = chroma_sp) +
  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
  theme_pubk() +
  theme(
    legend.position = c(0.2, 1.15),
    #axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    panel.spacing.x = unit(6, "mm"),
    panel.spacing.y = unit(2, "mm")
  ) +
  guides(
    fill = "none",
    size = guide_legend(
      override.aes = list(fill = "black"), 
      keyheight = 0.5,
      keywidth = 0.5,
      nrow = 1
    )
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    x = "Predictor value (m or deg C)",
    y = "Mole fraction",
    size = ''
  )
slide_headgp_pgls
ggsave(here("04-slidefigs", "PPEvDepthTemp_allsamps_20230923a.pdf"), width = 6, height = 3, units = "in")

# for acyl chains
# based on code for panel acylch_pgls
filter_panelacylch_pgls = function(x){
  x %>% 
    # only interested in aggregate chain trends
    filter(class == "all") %>% 
    # non-jackknifed dataset
    filter(jknife == "nothing") %>% 
    # and just class fractions
    #filter(str_detect(respvar, regex("sn[12]"))) %>% 
    filter(str_detect(respvar, "sn1")) %>%
    #filter(
    #  ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
    #    ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
    #) %>% 
    filter(
      ((predvar == "depth_col") & (subset == "<= 10 deg C")) |
        ((predvar == "temp_col") & (subset == "<= 250 m"))
    ) %>% 
    # add the faceting var why not?
    mutate(
      stereo = parse_number(respvar),
      respvar = str_remove(respvar, regex("[0-9]")),
      trend = str_glue("{respvar}~{predvar}\n{subset} sn{stereo}")
    ) %>% 
    #mutate(trend = paste(respvar, predvar, sep='~')) %>% 
    ungroup() %>% 
    #arrange(respvar, predvar, subset, stereo) %>% 
    arrange(predvar, respvar, subset, stereo) %>% 
    mutate(trend = trend %>% factor(levels = unique(.)))
}

panel_acylch_pgls = pldata_indls_opt %>% 
  filter_panelacylch_pgls() %>% 
  # calc means and std errors
  group_by(respvar, sp, subset, predvar, trend) %>% 
  summarize(
    n = n(),
    resp_serr = sd(resp, na.rm = TRUE)/sqrt(n),
    resp = mean(resp, na.rm = TRUE),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  facet_grid(rows = vars(respvar), cols = vars(predvar), scales = "free") +
  #facet_wrap(~trend, scales = "free", strip.position = "left", ncol = 2) +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    color = "white",
    size = 7/.pt,
    # no idea how this ordering works; trial and error
    x = c(2000, 0   , 1000, 10),
    y = c(17.3, 19.3, 0.5 , 1 ),
    data = allmods_multtest %>% 
      filter_panelacylch_pgls() %>% 
      # OLS y-intercept for locating text
      group_by(subset, predvar, respvar, trend, jknife) %>% 
      mutate(yint = estimate[[1]]) %>% 
      filter(term == "pred") %>% 
      group_by(respvar, predvar, subset, trend, yint) %>% 
      arrange(model) %>% 
      # I don't think mult testing corr is appropriate here...are CL and DBI really a family?
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')),
    aes(label = plabel),
    hjust = 0.0,
    vjust = 1.0,
    check_overlap = TRUE
  ) +
  # best fit lines
  geom_line(
    color = "white",
    size = 1/.pt,
    data = modelines_lmap %>% 
      filter_panelacylch_pgls() %>% 
      ungroup(),
    aes(
      linetype = model
    )
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey75",
    aes(
      ymin = resp - resp_serr,
      ymax = resp + resp_serr
    )
  ) +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey75",
    aes(
      xmin = pred - pred_serr,
      xmax = pred + pred_serr
    )
  ) +
  geom_point(
    shape = 21,
    color = "white",
    aes(
      size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = chroma_sp) +
  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
  theme_pubk() +
  theme(
    #axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    # replacing the y label with the strip text
    #axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(vjust = -0.5),
    strip.placement = "outside",
    legend.position = "none",
    panel.spacing.x = unit(6, "mm"),
    panel.spacing.y = unit(2, "mm")
  ) +
  guides(
    linetype = "none"
  ) +
  labs(
    x = "Predictor value"
  )
panel_acylch_pgls
ggsave(here("04-slidefigs", "panel_acylch_subset_pgls_20230924b.pdf"), width = 6, height = 5, units = "in")

# summary of all significant trends
sigcorrs = allmods_multtest %>% 
  ungroup() %>% 
  filter(term == "pred") %>% 
  filter(model == "modgls") %>% 
  filter(p.value <= 0.05) %>% 
  arrange(p.value) %>% 
  select(subset, class, predvar, respvar, estimate, p.value) %>% 
  print()
