## Scatter/bubble-plot lipidomes' mean curvature under multiple estimation schemes
## Barplots and lipid contributions to curvature are in `plot_curvature.R`

library(tidyverse)
library(here)
library(ggpubr)

source(here("03-scripts", "prep_curvature.R"))
source(here("03-scripts", "plot_helpers.R"))
source(here("03-scripts", "pgls_helpers.R"))

# The phylogeny we're gonna use
#NTS 20230722: `20230722_iq.tre` is a hacked-up little tree with the necessary
# short names subbed in manually. Regenerate.
file_tree = here("01-rawdata", "20230722_iq.tre")
# load the phylogeny so we can fix anything that's wrong
phylodata = ape::read.tree(file_tree)
#plot(phylodata, main = "Input tree") # have a quick look

# summarize down total curvature for each individual
plcidata_indls = plcidata_wildwhole %>% 
  # EXCLUDE TECH REPLICATES
  filter(!str_detect(eid, regex("[A-Z]$"))) %>% 
  group_by(scheme, sp, eid, depth_col, temp_col) %>%
  summarize(plci = sum(plci))

# an expanded version that to test conditional exclusion
# or species jackknifing
plcidata_indls_opt = plcidata_indls %>% 
  # define expansion axes here
  cross_join(crossing(
    # try filtering to cold samples, or not
    subset = c(
      "all samples", 
      "<= 10 deg C", 
      "<= 250 m" # can comment this to run half the regressions
      ), 
    # jackknife one species at a time!
    # Don't forget we still want an all-inclusive dataset, hence "nothing".
    jknife = c("nothing", plcidata_indls$sp %>% unique()) %>% factor(levels = .)
  )) %>% 
  # define conditions for expansion here
  # the order (and prob separate filter calls) is important
  filter(!((subset == "<= 10 deg C") & (temp_col > 10))) %>% 
  filter(!((subset == "<= 250 m") & (depth_col > 250))) %>% 
  filter(as.character(sp) != as.character(jknife)) %>% 
  # and finally, melt the predictor columns
  # it's important to be able to run multiple predictors
  pivot_longer(cols = contains("_col"), names_to = "predvar", values_to = "pred") %>% 
  # remove trends I don't need
  filter(
    !(((predvar == "depth_col") & (subset == "<= 250 m")) |
        ((predvar == "temp_col") & (subset == "<= 10 deg C")))
  )


### OK, here we go w/PGLS...

### This block preps the tree to work with any subset of pheno data
# Reduce the phylo- and phenodata to the intersection of their species
# give warnings if species in phenodata are not found in phylogeny
phenophylo = match_tree(
  plcidata_indls, # don't need to use opt here, takes up extra space
  phylodata, # This tree has replicates of some species, suffixed w/a digit
  tip2sp = function(x){str_remove_all(x, "[0-9]")}
)
pheno = phenophylo$pheno
phylo = phenophylo$phylo # 1 tip per species in this tree

# Midpoint-root the phylogeny
# This is a conservative choice, mitigating signal due to long-branch attraction
message("Midpoint-rooting phylogeny")
phylo_mid = phylo %>% midpoint.root()

# ultrametrize
message("Making phylogeny ultrametric")
phylo_ult = phylo %>% chronos()

# Plot the species tree
plot(phylo_ult, main = "Species tree")

# fit using new robust wrapper func that matches the tree on each call
# this takes awhile with species jackknifing
allmods_plci = plcidata_indls_opt %>% 
  group_by(subset, scheme, predvar, jknife) %>% 
  # can specify all the models I want to fit here!
  # to fit multiple models, need to stuff data into static listcol
  summarize(thesedata = cur_data() %>% list()) %>% 
  rowwise() %>% 
  # go parallel!
  group_split() %>% 
  # do I need to add safely()?
  # see prep_pgls_lipidmaps.R
  pbapply::pblapply(
    cl = 8L,
    X = .,
    FUN = function(row){
      row %>% mutate(
        # OLS
        modols = lm(
          plci ~ pred, 
          thesedata[[1]]
        ) %>% list(),
        # GLS Brownian?
        # GLS Martins
        modgls = pgls(
          plci ~ pred, 
          thesedata[[1]], 
          ape::corMartins, 
          phylo_ult
        ) %>% list()
      )
    }
  ) %>% 
  do.call(rbind, .) %>% 
  select(-thesedata) %>% 
  pivot_longer(cols = contains("mod"), names_to = "model", values_to = "fit") %>% 
  rowwise() %>% 
  # get the standard summary and the info criteria
  mutate(fitout = fit %>% tidyglance() %>% list()) %>% 
  unnest(fitout) %>% 
  # regroup
  group_by(subset, scheme, predvar, jknife, model)

# generate prediction lines based on all the models
modelines_plci = plcidata_indls_opt %>% 
  group_by(subset, scheme, predvar, jknife) %>% 
  # filter to extreme ends of each predictor (2 points make a line!)
  filter(pred %in% range(pred)) %>% 
  distinct(pred) %>% 
  # join automatically on grouping vars
  left_join(
    allmods_plci %>% summarize(fit = fit[[1]] %>% list()),
    relationship = "many-to-many"
  ) %>% 
  group_by(subset, scheme, predvar, model, jknife) %>% 
  mutate(plci = predict(fit[[1]], cur_data())) %>% 
  # for ease of viewing
  arrange(subset, scheme, predvar, model, pred)

# Simplified exploratory bubbleplot showing effect of temperature constraint,
# c0 prediction scheme, and phylogenetic correction
# bubbleplot species means, everyone
plot_c0_pgls = plcidata_indls_opt %>% 
  # non-jackknifed dataset
  filter(jknife == "nothing") %>% 
  # just by depth
  filter(predvar == "depth_col") %>% # comment this to see temp trends too
  # calc means and std errors
  group_by(scheme, sp, subset, predvar) %>% 
  summarize(
    n = n(),
    plci_serr = sd(plci)/sqrt(n),
    plci = mean(plci),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  ggplot(
    aes(
      x = pred,
      y = plci
    )
  ) +
  facet_grid(rows = vars(scheme), cols = vars(subset)) +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    x = 4000,
    y = 0.0025,
    data = allmods_plci %>% 
      filter(jknife == "nothing") %>% 
      filter(predvar == "depth_col") %>% 
      filter(term == "pred") %>% 
      group_by(scheme, subset) %>% 
      arrange(model) %>% 
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')),
    aes(label = plabel),
    hjust = 1,
    check_overlap = TRUE
  ) +
  # best fit lines
  geom_line(
    color = "black",
    data = modelines_plci %>% 
      filter(predvar == "depth_col") %>% 
      filter(jknife == "nothing") %>% 
      ungroup(),
    aes(linetype = model)
  ) +
  geom_errorbar(
    color = "black",
    aes(
      ymin = plci - plci_serr,
      ymax = plci + plci_serr
    )
  ) +
  geom_errorbarh(
    color = "black",
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
  theme_pubr() +
  guides(fill = "none") +
  labs(
    title = "Curvature vs. depth\nby dataset and estimation scheme",
    x = "Collection depth (m)",
    y = "Mean c0 (Å^-1)",
    size = "N"
  )
plot_c0_pgls
#ggsave(here("04-omitfigs", "c0_ols_cold_and_scheme_20230726a.pdf"), width = 6, height = 8)

## Exploratory bubbleplot to look at c0 as function of temp
#plot_c0_depthtemp = plcidata_indls_opt %>% 
#  # non-jackknifed dataset
#  filter(jknife == "nothing") %>% 
#  # just by depth
#  #filter(predvar == "depth_col") %>% # comment this to see temp trends too
#  mutate(scheme = paste(predvar, scheme)) %>% 
#  # calc means and std errors
#  group_by(scheme, sp, subset, predvar) %>% 
#  summarize(
#    n = n(),
#    plci_serr = sd(plci)/sqrt(n),
#    plci = mean(plci),
#    pred_serr = sd(pred)/n,
#    pred = mean(pred)
#  ) %>% 
#  ggplot(
#    aes(
#      x = pred,
#      y = plci
#    )
#  ) +
#  facet_grid(rows = vars(scheme), cols = vars(subset)) +
#  # Use the OLS results above, rather than geom_smooth!
#  # (the intra-sp distributions are not symmetrical!)
#  geom_text(
#    x = 4000,
#    y = 0.0025,
#    data = allmods_plci %>% 
#      filter(jknife == "nothing") %>% 
#      filter(predvar == "depth_col") %>% 
#      filter(term == "pred") %>% 
#      group_by(scheme, subset) %>% 
#      arrange(model) %>% 
#      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#      summarize(plabel = paste(plabel, collapse='\n')),
#    aes(label = plabel),
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  # best fit lines
#  geom_line(
#    color = "black",
#    data = modelines_plci %>% 
#      filter(predvar == "depth_col") %>% 
#      filter(jknife == "nothing") %>% 
#      ungroup(),
#    aes(linetype = model)
#  ) +
#  geom_errorbar(
#    color = "black",
#    aes(
#      ymin = plci - plci_serr,
#      ymax = plci + plci_serr
#    )
#  ) +
#  geom_errorbarh(
#    color = "black",
#    aes(
#      xmin = pred - pred_serr,
#      xmax = pred + pred_serr
#    )
#  ) +
#  geom_point(
#    shape = 21,
#    color = "white",
#    aes(
#      size = n,
#      fill = sp
#    )
#  ) +
#  scale_size_continuous(range = c(2, 6)) +
#  scale_fill_manual(values = chroma_sp) +
#  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
#  theme_pubr() +
#  guides(fill = "none") +
#  labs(
#    title = "Curvature vs. depth\nby dataset and estimation scheme",
#    x = "Collection depth (m)",
#    y = "Mean c0 (Å^-1)",
#    size = "N"
#  )
#plot_c0_depthtemp
#ggsave(here("04-omitfigs", "c0_pgls_temp_and_subset_20230920a.pdf"), width = 6, height = 8)


#panel_4c = plcidata_indls_opt %>% 
#  # non-jackknifed dataset
#  filter(jknife == "nothing") %>% 
#  # just by depth
#  filter(predvar == "depth_col") %>% 
#  # just one panel!
#  filter(subset == "<= 10 deg C") %>% 
#  filter(scheme == "linreg") %>% 
#  # calc means and std errors
#  group_by(scheme, sp, subset, predvar) %>% 
#  summarize(
#    n = n(),
#    plci_serr = sd(plci)/sqrt(n),
#    plci = mean(plci),
#    pred_serr = sd(pred)/n,
#    pred = mean(pred)
#  ) %>% 
#  # alleviate overplotting
#  arrange(-n) %>% 
#  ggplot(
#    aes(
#      x = pred,
#      y = plci
#    )
#  ) +
#  facet_grid(rows = vars(scheme), cols = vars(subset)) +
#  # Use the OLS results above, rather than geom_smooth!
#  # (the intra-sp distributions are not symmetrical!)
#  geom_text(
#    x = 4000,
#    y = 0.001,
#    size = 7/.pt,
#    data = allmods_plci %>% 
#      filter(jknife == "nothing") %>% 
#      filter(predvar == "depth_col") %>% 
#      # just one panel!
#      filter(subset == "<= 10 deg C") %>% 
#      filter(scheme == "linreg") %>% 
#      filter(term == "pred") %>% 
#      group_by(scheme, subset) %>% 
#      arrange(model) %>% 
#      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#      summarize(plabel = paste(plabel, collapse='\n')),
#    aes(label = plabel),
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  # best fit lines
#  geom_line(
#    color = "black",
#    size = 1/.pt,
#    data = modelines_plci %>% 
#      filter(predvar == "depth_col") %>% 
#      filter(jknife == "nothing") %>% 
#      # just one panel!
#      filter(subset == "<= 10 deg C") %>% 
#      filter(scheme == "linreg") %>% 
#      ungroup(),
#    aes(linetype = model)
#  ) +
#  geom_errorbar(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      ymin = plci - plci_serr,
#      ymax = plci + plci_serr
#    )
#  ) +
#  geom_errorbarh(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      xmin = pred - pred_serr,
#      xmax = pred + pred_serr
#    )
#  ) +
#  geom_point(
#    shape = 21,
#    color = "white",
#    aes(
#      size = n,
#      fill = sp
#    )
#  ) +
#  scale_size_continuous(range = c(1, 3), breaks = seq(1, 7, 2)) +
#  scale_fill_manual(values = chroma_sp) +
#  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
#  theme_tiny() +
#  theme(
#    strip.background = element_blank(),
#    strip.text = element_blank(),
#    legend.position = c(0.4, 0.125),
#    #legend.position = "none",
#    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
#    axis.text.y = element_text(angle = 90, hjust = 0.5),
#  ) +
#  guides(
#    fill = "none",
#    linetype = "none",
#    size = guide_legend(
#      title.position = "left",
#      override.aes = list(fill = "black"), 
#      keyheight = 0.5,
#      keywidth = 0.5,
#      nrow = 1
#    ),
#  ) +
#  labs(
#    x = "Collection depth (m)",
#    y = "Mean c0 at 1 bar (Å^-1)",
#    size = "N"
#  )
#panel_4c
#ggsave(here("04-mainfigs", "Fig4_homeocurvature", "panel_4c_20231003a.pdf"), width = 40, height = 40, unit = "mm")

# New 4C 20231003: all samples with temperature
# 
filter_panel4c = function(datin){
  datin %>% 
    # non-jackknifed dataset
    filter(jknife %in% c("nothing")) %>% 
    ## just by depth
    #filter(predvar == "depth_col") %>% 
    filter(subset == "all samples") %>% 
    #filter(
    #  ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
    #    ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
    #) %>% 
    # just linreg estimation scheme
    filter(scheme == "linreg")
}

panel_4c = plcidata_indls_opt %>% 
  filter_panel4c() %>% 
  # don't want duplicated data pts in plot!
  filter(jknife %in% c("nothing")) %>% 
  # calc means and std errors
  group_by(scheme, sp, subset, predvar) %>% 
  summarize(
    n = n(),
    plci_serr = sd(plci)/sqrt(n),
    plci = mean(plci),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  mutate(correl = paste(subset, predvar)) %>% # need to order factor?
  ggplot(
    aes(
      x = pred,
      y = plci
    )
  ) +
  #facet_grid(rows = vars(subset), cols = vars(predvar)) +
  facet_wrap(~correl, ncol = 2, scales = "free_x") +
  ## Use the OLS results above, rather than geom_smooth!
  ## (the intra-sp distributions are not symmetrical!)
  geom_text(
    x = c(4000, 25),
    y = 0.002,
    alpha = 1,
    size = 7/.pt,
    data = allmods_plci %>% 
      filter_panel4c() %>% 
      filter(term == "pred") %>% 
      group_by(scheme, predvar, subset, jknife) %>% 
      arrange(model) %>% 
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')) %>% 
      mutate(correl = paste(subset, predvar)),
    aes(label = plabel),
    hjust = 1,
    check_overlap = TRUE
  ) +
  # best fit lines, separate alphas for jackknife
  geom_line(
    color = "black",
    data = modelines_plci %>% 
      filter_panel4c() %>%  
      mutate(correl = paste(subset, predvar)) %>% 
      ungroup(),
    aes(
      linetype = model,
      #group = jknife,
      #alpha = jknife
    )
  ) +
  geom_errorbar(
    color = "grey25",
    size = 0.5/.pt,
    aes(
      ymin = plci - plci_serr,
      ymax = plci + plci_serr
    )
  ) +
  geom_errorbarh(
    color = "grey25",
    size = 0.5/.pt,
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
  scale_size_continuous(range = c(1, 4), breaks = seq(1,7,2)) +
  scale_fill_manual(values = chroma_sp) +
  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
  #scale_alpha_manual(values = c("nothing" = 1, "Tjal_pink" = 0.5)) +
  theme_tiny() +
  theme(
    legend.position = c(0.20, 0.15),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1.25),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
  ) +
  guides(
    fill = "none",
    linetype = "none",
    size = guide_legend(
      override.aes = list(fill = "black"), 
      title.position = "left",
      keyheight = 0.5,
      keywidth = 0.5,
      nrow = 1
    )
  ) +
  labs(
    x = "Predictor value",
    y = "Mean c0 (Å^-1)",
    size = "N",
  )
panel_4c
ggsave(here("04-mainfigs", "Fig4_homeocurvature", "panel_4d_20240405a.pdf"), width = 60, height = 40, unit = "mm")


# Panel S7A-D: 2x2 grid with subsets and jackknives
filter_panels7abcd = function(datin){
  # just linreg estimation scheme
  datin %>% filter(scheme == "linreg") %>% 
  # constrain subsets
  filter(
    ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
      ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
  ) %>% 
    # constrain jackknives
    filter(
      ((jknife %in% c("nothing", "Tjal_pink")) & (subset %in% c("all samples", "<= 10 deg C"))) |
        ((jknife %in% c("nothing", "Coel_hali")) & (subset %in% c("all samples", "<= 250 m")))
    ) %>% 
    # add a boolean variable for subsetted or no
    mutate(subsetted = (subset != "all samples"))
}

panel_s7abcd = plcidata_indls_opt %>% 
  # non-jackknifed dataset; don't want to duplicate the actual data pts!
  filter(jknife %in% c("nothing")) %>% 
  # apply the standard filters
  filter_panels7abcd() %>% 
  # calc means and std errors
  group_by(scheme, sp, subset, subsetted, predvar) %>% 
  summarize(
    n = n(),
    plci_serr = sd(plci)/sqrt(n),
    plci = mean(plci),
    pred_serr = sd(pred)/n,
    pred = mean(pred)
  ) %>% 
  ggplot(
    aes(
      x = pred,
      y = plci
    )
  ) +
  facet_grid(rows = vars(subsetted), cols = vars(predvar), scales = "free_x") +
  # Use the OLS results above, rather than geom_smooth!
  # (the intra-sp distributions are not symmetrical!)
  geom_text(
    y = 0.002,
    alpha = 1,
    size = 7/.pt,
    data = allmods_plci %>% 
      filter_panels7abcd() %>% 
      filter(jknife == "nothing") %>% 
      filter(term == "pred") %>% 
      group_by(scheme, subset, subsetted, predvar, jknife) %>% 
      arrange(model) %>% 
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')),
    aes(
      group = predvar,
      x = ifelse(predvar=="depth_col", 4000, 26),
      label = plabel
    ),
    hjust = 1,
    check_overlap = TRUE
  ) +
  geom_text(
    y = -0.0045,
    alpha = 0.5,
    size = 7/.pt,
    data = allmods_plci %>% 
      filter_panels7abcd() %>% 
      filter(jknife == "Tjal_pink") %>% 
      filter(term == "pred") %>% 
      group_by(scheme, subset, subsetted, predvar, jknife) %>% 
      arrange(model) %>% 
      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
      summarize(plabel = paste(plabel, collapse='\n')),
    aes(
      group = predvar,
      x = ifelse(predvar=="depth_col", 4000, 26),
      label = plabel
    ),
    hjust = 1,
    check_overlap = TRUE
  ) +
  # best fit lines
  # handling alpha as an aesthetic here
    geom_line(
      color = "black",
      size = 1/.pt,
      data = modelines_plci %>% 
        filter_panels7abcd() %>% 
        filter(jknife %in% c("nothing", "Tjal_pink")) %>% 
        ungroup(),
      aes(
        group = paste(jknife, model),
        linetype = model,
        alpha = jknife
      )
    ) +
    geom_errorbar(
      color = "grey25",
      size = 0.5/.pt,
      aes(
        ymin = plci - plci_serr,
        ymax = plci + plci_serr
      )
    ) +
    geom_errorbarh(
      color = "grey25",
      size = 0.5/.pt,
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
    #scale_color_manual(values = c(chroma_sp, "nothing" = "transparent")) +
    scale_alpha_manual(values = c("nothing"=1, "Tjal_pink"=0.5)) +
    scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
    #scale_alpha_manual(values = c("nothing" = 1, "Tjal_pink" = 0.5)) +
    theme_tiny() +
    theme(
      legend.position = c(0.75, 0.075),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
    ) +
    guides(
      fill = "none",
      linetype = "none",
      alpha = "none",
      size = guide_legend(
        override.aes = list(fill = "black"), 
        title.position = "left",
        keyheight = 0.5,
        keywidth = 0.5,
        nrow = 1
      )
    ) +
    labs(
      x = "Collection depth (m)",
      y = "Mean c0 (Å^-1)",
      size = "N",
    )
panel_s7abcd
ggsave(here("04-suppfigs", "panel_s7abcd_20231004d.pdf"), width = 120, height = 120, unit = "mm")

# fancy colorful-lines version
#panel_s7abcd = plcidata_indls_opt %>% 
#  # non-jackknifed dataset; don't want to duplicate the actual data pts!
#  filter(jknife %in% c("nothing")) %>% 
#  # apply the standard filters
#  filter_panels7abcd() %>% 
#  # calc means and std errors
#  group_by(scheme, sp, subset, subsetted, predvar) %>% 
#  summarize(
#    n = n(),
#    plci_serr = sd(plci)/sqrt(n),
#    plci = mean(plci),
#    pred_serr = sd(pred)/n,
#    pred = mean(pred)
#  ) %>% 
#  ggplot(
#    aes(
#      x = pred,
#      y = plci
#    )
#  ) +
#  facet_grid(rows = vars(subsetted), cols = vars(predvar), scales = "free_x") +
#  ##facet_grid(~subset) +
#  ## Use the OLS results above, rather than geom_smooth!
#  ## (the intra-sp distributions are not symmetrical!)
#  #geom_text(
#  #  x = 4000,
#  #  y = 0.002,
#  #  alpha = 1,
#  #  size = 7/.pt,
#  #  data = allmods_plci %>% 
#  #    filter(jknife == "nothing") %>% 
#  #    filter(scheme == "linreg") %>% 
#  #    filter(predvar == "depth_col") %>% 
#  #    filter(term == "pred") %>% 
#  #    group_by(scheme, subset, jknife) %>% 
#  #    arrange(model) %>% 
#  #    mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#  #    summarize(plabel = paste(plabel, collapse='\n')),
#  #  aes(label = plabel),
#  #  hjust = 1,
#  #  check_overlap = TRUE
#  #) +
#  #geom_text(
#  #  x = 4000,
#  #  y = -0.0045,
#  #  alpha = 0.5,
#  #  size = 7/.pt,
#  #  data = allmods_plci %>% 
#  #    filter(jknife == "Tjal_pink") %>% 
#  #    filter(scheme == "linreg") %>% 
#  #    filter(predvar == "depth_col") %>% 
#  #    filter(term == "pred") %>% 
#  #    group_by(scheme, subset, jknife) %>% 
#  #    arrange(model) %>% 
#  #    mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#  #    summarize(plabel = paste(plabel, collapse='\n')),
#  #  aes(label = plabel),
#  #  hjust = 1,
#  #  check_overlap = TRUE
#  #) +
#  # best fit lines
#  geom_line(
#    size = 1/.pt,
#    data = modelines_plci %>% 
#      filter_panels7abcd() %>% 
#      ungroup(),
#    aes(
#      group = paste(jknife, model),
#      color = jknife
#    )
#  ) +
#  geom_line(
#    color = "black",
#    size = 0.5/.pt,
#    data = modelines_plci %>% 
#      filter_panels7abcd() %>% 
#      ungroup(),
#    aes(
#      group = paste(jknife, model),
#      linetype = model
#    )
#  ) +
#  geom_errorbar(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      ymin = plci - plci_serr,
#      ymax = plci + plci_serr
#    )
#  ) +
#  geom_errorbarh(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      xmin = pred - pred_serr,
#      xmax = pred + pred_serr
#    )
#  ) +
#  geom_point(
#    shape = 21,
#    color = "white",
#    aes(
#      size = n,
#      fill = sp
#    )
#  ) +
#  scale_size_continuous(range = c(1, 3), breaks = seq(1,7,2)) +
#  scale_fill_manual(values = chroma_sp) +
#  scale_color_manual(values = c(chroma_sp, "nothing" = "transparent")) +
#  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
#  #scale_alpha_manual(values = c("nothing" = 1, "Tjal_pink" = 0.5)) +
#  theme_tiny() +
#  theme(legend.position = c(0.20, 0.15)) +
#  guides(
#    fill = "none",
#    linetype = "none",
#    size = guide_legend(
#      override.aes = list(fill = "black"), 
#      title.position = "left",
#      keyheight = 0.5,
#      keywidth = 0.5,
#      nrow = 1
#    )
#  ) +
#  labs(
#    x = "Collection depth (m)",
#    y = "Mean c0 (Å^-1)",
#    size = "N",
#  )
#panel_s7abcd
#ggsave(here("04-suppfigs", "panel_s7abcd_20231004a.pdf"), width = 80, height = 80, unit = "mm")

## Panel S7: <=10°C and all samples, with and without Tjalfie
#panel_s7ab = plcidata_indls_opt %>% 
#  # non-jackknifed dataset; don't want to duplicate the actual data pts!
#  filter(jknife %in% c("nothing")) %>% 
#  # just by depth
#  filter(predvar == "depth_col") %>% 
#  # just linreg estimation scheme
#  filter(scheme == "linreg") %>% 
#  # calc means and std errors
#  group_by(scheme, sp, subset, predvar) %>% 
#  summarize(
#    n = n(),
#    plci_serr = sd(plci)/sqrt(n),
#    plci = mean(plci),
#    pred_serr = sd(pred)/n,
#    pred = mean(pred)
#  ) %>% 
#  ggplot(
#    aes(
#      x = pred,
#      y = plci
#    )
#  ) +
#  #facet_grid(rows = vars(scheme), cols = vars(jknife)) +
#  facet_grid(~subset) +
#  # Use the OLS results above, rather than geom_smooth!
#  # (the intra-sp distributions are not symmetrical!)
#  geom_text(
#    x = 4000,
#    y = 0.002,
#    alpha = 1,
#    size = 7/.pt,
#    data = allmods_plci %>% 
#      filter(jknife == "nothing") %>% 
#      filter(scheme == "linreg") %>% 
#      filter(predvar == "depth_col") %>% 
#      filter(term == "pred") %>% 
#      group_by(scheme, subset, jknife) %>% 
#      arrange(model) %>% 
#      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#      summarize(plabel = paste(plabel, collapse='\n')),
#    aes(label = plabel),
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  geom_text(
#    x = 4000,
#    y = -0.0045,
#    alpha = 0.5,
#    size = 7/.pt,
#    data = allmods_plci %>% 
#      filter(jknife == "Tjal_pink") %>% 
#      filter(scheme == "linreg") %>% 
#      filter(predvar == "depth_col") %>% 
#      filter(term == "pred") %>% 
#      group_by(scheme, subset, jknife) %>% 
#      arrange(model) %>% 
#      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#      summarize(plabel = paste(plabel, collapse='\n')),
#    aes(label = plabel),
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  # best fit lines, separate alphas for jackknife
#  geom_line(
#    color = "black",
#    alpha = 0.5,
#    data = modelines_plci %>% 
#      filter(jknife == "Tjal_pink") %>% 
#      filter(scheme == "linreg") %>% 
#      filter(predvar == "depth_col") %>% 
#      ungroup(),
#    aes(
#      linetype = model,
#      #group = jknife,
#      #alpha = jknife
#    )
#  ) +
#  geom_line(
#    color = "black",
#    data = modelines_plci %>% 
#      filter(jknife == "nothing") %>% 
#      filter(scheme == "linreg") %>% 
#      filter(predvar == "depth_col") %>% 
#      ungroup(),
#    aes(
#      linetype = model,
#      #group = jknife,
#      #alpha = jknife
#    )
#  ) +
#  geom_errorbar(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      ymin = plci - plci_serr,
#      ymax = plci + plci_serr
#    )
#  ) +
#  geom_errorbarh(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      xmin = pred - pred_serr,
#      xmax = pred + pred_serr
#    )
#  ) +
#  geom_point(
#    shape = 21,
#    color = "white",
#    aes(
#      size = n,
#      fill = sp
#    )
#  ) +
#  scale_size_continuous(range = c(1, 3), breaks = seq(1,7,2)) +
#  scale_fill_manual(values = chroma_sp) +
#  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
#  #scale_alpha_manual(values = c("nothing" = 1, "Tjal_pink" = 0.5)) +
#  theme_tiny() +
#  theme(legend.position = c(0.20, 0.15)) +
#  guides(
#    fill = "none",
#    linetype = "none",
#    size = guide_legend(
#      override.aes = list(fill = "black"), 
#      title.position = "left",
#      keyheight = 0.5,
#      keywidth = 0.5,
#      nrow = 1
#    )
#  ) +
#  labs(
#    x = "Collection depth (m)",
#    y = "Mean c0 (Å^-1)",
#    size = "N",
#  )
#panel_s7ab
#ggsave(here("04-suppfigs", "panel_s7a_20231003a.pdf"), width = 120, height = 60, unit = "mm")
#
## Panel S7 with temp as well
#curvslide = plcidata_indls_opt %>% 
#  # non-jackknifed dataset
#  filter(jknife %in% c("nothing")) %>% 
#  ## just by depth
#  #filter(predvar == "depth_col") %>% 
#  filter(
#    ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
#      ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
#  ) %>% 
#  # just linreg estimation scheme
#  filter(scheme == "linreg") %>% 
#  # calc means and std errors
#  group_by(scheme, sp, subset, predvar) %>% 
#  summarize(
#    n = n(),
#    plci_serr = sd(plci)/sqrt(n),
#    plci = mean(plci),
#    pred_serr = sd(pred)/n,
#    pred = mean(pred)
#  ) %>% 
#  mutate(correl = paste(subset, predvar)) %>% # need to order factor?
#  ggplot(
#    aes(
#      x = pred,
#      y = plci
#    )
#  ) +
#  #facet_grid(rows = vars(subset), cols = vars(predvar)) +
#  facet_wrap(~correl, ncol = 2, scales = "free_x") +
#  ## Use the OLS results above, rather than geom_smooth!
#  ## (the intra-sp distributions are not symmetrical!)
#  geom_text(
#    x = c(4000, 4000, 25, 25),
#    y = 0.002,
#    alpha = 1,
#    size = 7/.pt,
#    data = allmods_plci %>% 
#      filter(jknife == "nothing") %>% 
#      filter(scheme == "linreg") %>% 
#      #filter(predvar == "depth_col") %>% 
#      filter(
#        ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
#          ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
#      ) %>% 
#      filter(term == "pred") %>% 
#      group_by(scheme, predvar, subset, jknife) %>% 
#      arrange(model) %>% 
#      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#      summarize(plabel = paste(plabel, collapse='\n')) %>% 
#      mutate(correl = paste(subset, predvar)),
#    aes(label = plabel),
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  geom_text(
#    x = c(4000, 4000, 25, 25),
#    y = -0.0045,
#    alpha = 0.5,
#    size = 7/.pt,
#    data = allmods_plci %>% 
#      filter(jknife == "Tjal_pink") %>% 
#      filter(scheme == "linreg") %>% 
#      #filter(predvar == "depth_col") %>% 
#      filter(
#        ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
#          ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
#      ) %>% 
#      filter(term == "pred") %>% 
#      group_by(scheme, predvar, subset, jknife) %>% 
#      arrange(model) %>% 
#      mutate(plabel = str_glue({"{model} P = {format(p.value, digits=3, scientific=TRUE)}"})) %>% 
#      summarize(plabel = paste(plabel, collapse='\n')) %>% 
#      mutate(correl = paste(subset, predvar)),
#    aes(label = plabel),
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  # best fit lines, separate alphas for jackknife
#  geom_line(
#    color = "black",
#    alpha = 0.5,
#    data = modelines_plci %>% 
#      filter(jknife == "Tjal_pink") %>% 
#      filter(scheme == "linreg") %>% 
#      #filter(predvar == "depth_col") %>% 
#      filter(
#        ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
#          ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
#      ) %>% 
#      mutate(correl = paste(subset, predvar)) %>% 
#      ungroup(),
#    aes(
#      linetype = model,
#      #group = jknife,
#      #alpha = jknife
#    )
#  ) +
#  geom_line(
#    color = "black",
#    data = modelines_plci %>% 
#      filter(jknife == "nothing") %>% 
#      filter(scheme == "linreg") %>% 
#      #filter(predvar == "depth_col") %>% 
#      filter(
#        ((predvar == "depth_col") & (subset %in% c("all samples", "<= 10 deg C"))) |
#          ((predvar == "temp_col") & (subset %in% c("all samples", "<= 250 m")))
#      ) %>% 
#      mutate(correl = paste(subset, predvar)) %>% 
#      ungroup(),
#    aes(
#      linetype = model,
#      #group = jknife,
#      #alpha = jknife
#    )
#  ) +
#  geom_errorbar(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      ymin = plci - plci_serr,
#      ymax = plci + plci_serr
#    )
#  ) +
#  geom_errorbarh(
#    color = "grey25",
#    size = 0.5/.pt,
#    aes(
#      xmin = pred - pred_serr,
#      xmax = pred + pred_serr
#    )
#  ) +
#  geom_point(
#    shape = 21,
#    color = "white",
#    aes(
#      size = n,
#      fill = sp
#    )
#  ) +
#  scale_size_continuous(range = c(1, 3), breaks = seq(1,7,2)) +
#  scale_fill_manual(values = chroma_sp) +
#  scale_linetype_manual(values = c("modols" = "solid", "modgls" = "dashed")) +
#  #scale_alpha_manual(values = c("nothing" = 1, "Tjal_pink" = 0.5)) +
#  theme_tiny() +
#  theme(legend.position = c(0.20, 0.15)) +
#  guides(
#    fill = "none",
#    linetype = "none",
#    size = guide_legend(
#      override.aes = list(fill = "black"), 
#      title.position = "left",
#      keyheight = 0.5,
#      keywidth = 0.5,
#      nrow = 1
#    )
#  ) +
#  labs(
#    x = "Predictor value",
#    y = "Mean c0 (Å^-1)",
#    size = "N",
#  )
#curvslide
#ggsave(here("04-suppfigs", "curvslide_20231003b.pdf"), width = 120, height = 120, unit = "mm")

## Claims about the species jackknife analysis
jknife_pvals_plci = allmods_plci %>% 
  #filter(subset == "<= 10 deg C") %>% 
  filter(scheme == "linreg") %>% 
  filter((term == "pred") & (predvar == "depth_col")) %>% 
  #filter(model == "modgls") %>% 
  arrange(-p.value, model, scheme, jknife)

# save just the P-values to TSV for supp table
# NOTE could omit jackknifes that improve the signal to save space...
jknife_pvals_plci %>% 
  #group_by(jknife) %>% 
  arrange(-p.value) %>% 
  # refactor jknife to go in desc order of P-modgls
  mutate(jknife = jknife %>% factor(., levels = unique(.))) %>% 
  arrange(jknife, model) %>% 
  #filter(subset == "<= 10 deg C") %>% 
  filter(subset == "all samples") %>%
  #pivot_wider(
  #  id_cols = jknife, 
  #  names_from = "subset", 
  #  values_from = "p.value",
  #  names_prefix = "p.value "
  #) %>% 
  # clean up
  ungroup() %>% 
  select(jknife, model, p.value, adj.r.squared) %>% 
  write_tsv(here("04-suppfigs", "plci_allsamp_jackknife_20231003a.tsv")) %>% 
  # see most impactful omissions first
  arrange(-p.value) %>% 
  print()

# idk why this changed a little when I did away with Bero_arct
# But the outcome is the same
jknife_pvals_plci %>% filter(subset == "<= 10 deg C") %>% .$p.value %>% max()
# yup, remove anyone! :D