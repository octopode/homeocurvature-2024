library(tidyverse)
library(here)
library(ggtree)
library(ggpubr)

#source(here("03-scripts", "prep_lipidmaps.R"))
# this is for the tree in Fig. 3B
source(here("03-scripts", "prep_pgls_lipidmaps.R")) 

# "Golden 5" (or 6 or 7) species to show in high-detail plots
# in desired plotting order
sp_subset = c(
  "Boli_vitr", 
  "Boli_infu", 
  "Boli_arct", 
  "Lamp_crue", 
  "Tjal_pink"
)

# PL individuals: a spot check
pldata_wildwhole %>% 
  #filter(eid != "JWL0321") %>% # one with a LOT of PPE and little PS
  #filter(sp %in% sp_subset) %>% 
  #mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  arrange(sp) %>% 
  mutate(
    sp_eid = paste(sp, eid),
    sp_eid = sp_eid %>% factor(levels = unique(sp_eid))
  ) %>% 
  # comment to show individual compounds
  #group_by(eid, sp, class) %>% 
  #summarise(frac_molar = sum(frac_molar)) %>% 
  gg_headgp(
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  #facet_grid(cols = vars(sp), scales="free_x") +
  labs(
    title = "PL class compositions",
    x = "Sample",
    y = "Mole fraction of phospholipids"
  )
ggsave(here("04-omitfigs", "pl_indls_all_20230814a.pdf"), width = 14, height = 5)

# PL individuals: a spot check
pldata_wildwhole %>% 
  filter(sp == "Leuc_pulc") %>% 
  #filter(eid != "JWL0321") %>% # one with a LOT of PPE and little PS
  #filter(sp %in% sp_subset) %>% 
  #mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  arrange(sp) %>% 
  mutate(
    sp_eid = paste(sp, eid),
    sp_eid = sp_eid %>% factor(levels = unique(sp_eid))
  ) %>% 
  # comment to show individual compounds
  #group_by(eid, sp, class) %>% 
  #summarise(frac_molar = sum(frac_molar)) %>% 
  gg_headgp(
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  #facet_grid(cols = vars(sp), scales="free_x") +
  labs(
    title = "PL class compositions",
    x = "Sample",
    y = "Mole fraction of phospholipids"
  )
ggsave(here("04-omitfigs", "pl_indl_leucos_20230814a.pdf"), width = 5, height = 6)

# sterols
plstdata %>% 
  select(-annot) %>% 
  drop_na(sp) %>% 
  mutate(sp_eid = paste(sp, eid)) %>% 
  gg_headgp(
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "PLs + sterols, unblanked",
    x = "Sample",
    y = "Mole fraction"
  )
ggsave(here("04-suppfigs", "sterols_pls_20230809b.pdf"), width = 8, height = 5)

# barplots of sterol fractions for supp
plstdata %>% 
  select(-annot) %>% 
  drop_na(sp) %>% 
  filter(class == "ST") %>% 
  mutate(sp_eid = paste(sp, eid)) %>% 
  ggplot(
    aes(
      x = sp_eid,
      y = frac_molar,
      fill = class
    )
  ) +
  geom_col() +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_fill_manual(values = chroma_cl) +
  labs(
    title = "Cholesterol mole fractions, individual",
    x = "Sample",
    y = "Mole fraction"
  )
ggsave(here("04-omitfigs", "sterols_indls_20231012a.pdf"), width = 180, height = 80, units = "mm")

# barplots of SM fractions
plsldata %>% 
  #select(-annot) %>% 
  drop_na(sp) %>% 
  filter(class == "SM") %>% 
  mutate(sp_eid = paste(sp, eid)) %>% 
  #ggplot(
  #  aes(
  #    x = sp_eid,
  #    y = frac_molar,
  #    fill = class
  #  )
  #) +
  gg_headgp(
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  #scale_fill_manual(values = chroma_cl) +
  labs(
    title = "SM composition, individual",
    x = "Sample",
    y = "Mole fraction of SMs"
  )
ggsave(here("04-omitfigs", "SM_indls_20231012b.pdf"), width = 180, height = 80, units = "mm")

# mean/max/SEM of SM fractions
plsldata %>% 
  #select(-annot) %>% 
  drop_na(sp) %>% 
  filter(class == "PG") %>% 
  group_by(eid) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  .$frac_molar %>% mean(.)
  #{sd(.)/sqrt(length(.))}

# chol barplot by sp
plstdata %>% 
  select(-annot) %>% 
  drop_na(sp) %>% 
  group_by(sp, class) %>% 
  summarize(
    serr_molar = sd(frac_molar)/sqrt(n()),
    frac_molar = mean(frac_molar),
    n = n()
  ) %>% 
  filter(class == "ST") %>% 
  ggplot(
    aes(
      x = sp,
      y = frac_molar,
      fill = class
    )
  ) +
  geom_col() +
  geom_errorbar(
    width = 0.2,
    aes(
      ymin = frac_molar - serr_molar,
      ymax = frac_molar + serr_molar
    )
  ) +
  geom_text(
    y = 0.0925,
    color = "black",
    size = 7/.pt,
    aes(label = n)
  ) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_fill_manual(values = chroma_cl) +
  labs(
    title = "Cholesterol mole fractions by species",
    x = "Sample",
    y = "Mole fraction"
  ) +
  lims(y = c(0, 0.0925))
ggsave(here("04-suppfigs", "sterols_spp_20230809b.pdf"), width = 120, height = 80, units = "mm")

# sphingolipids
plsldata %>% 
  #select(-annot) %>% 
  mutate(sp_eid = paste(sp, eid)) %>% 
  gg_headgp(
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "PLs + sphingolipids",
    x = "Sample",
    y = "Mole fraction"
  )
ggsave(here("04-omitfigs", "sphingos_pls_20231012b.pdf"), width = 14, height = 5)

# sphingos + sterols
plstldata %>% 
  mutate(sp_eid = paste(sp, eid)) %>% 
  gg_headgp(
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  labs(
    title = "PLs + sterols + sphingolipids",
    x = "Sample",
    y = "Mole fraction"
  )
ggsave(here("04-omitfigs", "sphingos_sterols_pls_20231012.pdf"), width = 8, height = 5)

## MAIN FIG PANELS ##

# Panel 3b tree, flipping tips
# Some species sorted in a warm -> cold -> shal -> deep order
#sp_inorder = c(
#  "Boli_vitr",
#  #"Bero_ovat",
#  "Coel_hali",
#  "Mert_angu",
#  "Boli_micr",
#  #"Bero_pseu",
#  "Boli_infu", 
#  #"Bero_cucu",
#  "Lamp_crue",
#  #"Bath_fost",
#  #"Bero_abys", # might wanna snap out the outlier here
#  "Cydi_blac",
#  "Tjal_pink"
#) %>% rev() # it actually plots from bottom to top

sp_inorder = c(
  "Coel_hali",
  "Tjal_pink",
  "Mert_angu",
  "Cydi_blac",
  "Boli_vitr",
  "Boli_micr",
  "Boli_infu", 
  "Lamp_crue"
) %>% rev()

# prep lipid and envi data for the tree+barplot panel
# doing a few tricky things here
# (1) binding available cholesterol data to _all_ PL data
# (2) reordering the species factor to match the tree tips

# Prep tree here
# The write-read is a kludge to strip the "chronos" data
trefile = here("02-tidydata", "species.tre")
phylo_ult %>% write.tree(trefile)
phylo_3b = trefile %>% 
  read.tree() %>% 
  ape::drop.tip(setdiff(.$tip.label, sp_inorder)) %>% 
  # put tips in the closest-possible-to desired order
  ape::rotateConstr(sp_inorder)

# Prep data
# we could include chol data if we had it for Lampo...
# ...and all the coldwater Bolinopses >.<
# some formatting params
ratio_pheno_phylo = 1.75 # how wide is the barplot relative to the tree?
pheno_start_x     = 1.75 # where does the barplot start? (tree always ends at 1)
#data_panel3b = plstsp %>% 
data_panel3b = pldata_wildwhole %>% 
  #filter(eid != "JWL0343") %>% #TEST
  filter(sp %in% sp_inorder) %>% 
  mutate(sp = sp %>% factor(levels = tipOrder(phylo_3b))) %>% 
  group_by(eid, depth_col, temp_col, sp, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  group_by(sp, class) %>% 
  summarise(
    n = n(),
    across(contains("_col"), c("min"=min, "max"=max, "mean"=mean)),
    frac_molar = mean(frac_molar)
  ) %>% 
  # *ensure* labels all agree
  group_by(sp) %>% 
  mutate(
    n = max(n),
    across(contains("_min"), min),
    across(contains("_max"), max),
    across(contains("_col_"), round), # to integer?
    lab = str_glue("N = {n}\n{depth_col_min} - {depth_col_max} m\n{temp_col_min} - {temp_col_max}°C"),
    # labels with just the means
    mlab = str_glue("{depth_col_mean} m\n{temp_col_mean}°C")
  ) %>% 
  # stretch lipid data out to compress the tree (which has width 1)
  mutate(frac_molar = ratio_pheno_phylo*frac_molar) %>% 
  # hack to shift the whole panel to the right
  bind_rows(
    .,
    tibble(sp = .$sp %>% unique()) %>% 
      mutate(
        class = "spacer",
        frac_molar = pheno_start_x # the offset
      )
  ) %>% 
  # as last level, it plots on the left
  mutate(class = class %>% factor(levels = c(levels(pldata$class), "spacer")))

panel_3b = phylo_3b %>% 
  gg_lollitree(linewidth = 0.5/.pt) +
  geom_col(
    data = data_panel3b,
    aes(
      x = frac_molar,
      y = sp %>% as.integer(),
      fill = class
    ),
    orientation = "y",
    width = 0.6
  ) +
  # mean depth and temp labels: two separate columns
  #geom_text(
  #  data = data_panel3b,
  #  x = 1 + (pheno_start_x-1)/3, 
  #  size = 7/.pt,
  #  #fontface = "bold",
  #  hjust = 0.5,
  #  lineheight = 0.8,
  #  check_overlap = TRUE,
  #  aes(
  #    label = depth_col_mean,
  #    y = sp %>% as.integer()
  #  )
  #) +
  #geom_text(
  #  data = data_panel3b,
  #  x = 1 + 2*(pheno_start_x-1)/3, 
  #  size = 7/.pt,
  #  #fontface = "bold",
  #  hjust = 0.5,
  #  lineheight = 0.8,
  #  check_overlap = TRUE,
  #  aes(
  #    label = temp_col_mean,
  #    y = sp %>% as.integer()
  #  )
  #) +
  # trying depth/temp w/bullet separator instead
  geom_text(
    data = data_panel3b,
    x = 1 + 1.5*(pheno_start_x-1)/3, 
    size = 7/.pt,
    #fontface = "bold",
    hjust = 0.5,
    lineheight = 0.8,
    check_overlap = TRUE,
    aes(
      label = paste(depth_col_mean, temp_col_mean, sep=' • '),
      y = sp %>% as.integer()
    )
  ) +
  # N labels
  # mean depth and temp labels
  geom_text(
    data = data_panel3b,
    x = pheno_start_x + ratio_pheno_phylo + (pheno_start_x-1)/6, # just beyond end of barplot
    size = 7/.pt,
    #fontface = "bold",
    hjust = 0.5,
    lineheight = 0.8,
    check_overlap = TRUE,
    aes(
      label = n,
      y = sp %>% as.integer()
    )
  ) +
  scale_fill_manual(values = c(chroma_sp, chroma_cl, c("spacer" = "transparent"))) +
  #theme_tiny() +
  theme(
    legend.position = "none",
    plot.margin = unit(rep(0, 4), "mm"),
  ) +
  lims(x = c(0,pheno_start_x + ratio_pheno_phylo + (pheno_start_x-1)/3))
panel_3b
ggsave(ggsave(here("04-mainfigs", "Fig3_ctenolipidomics", "panel_3b_20230920a.pdf"), width = 75, height = 40, units = "mm"))

# PL species means summed by class
panel_3a = pldata_wildwhole %>% 
  filter(sp %in% sp_subset) %>% 
  mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  group_by(eid, sp, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  group_by(sp, class) %>% 
  summarise(
    frac_molar = mean(frac_molar),
    n = {str_glue("N = {n()}")}
  ) %>% 
  gg_headgp(
    x = sp,
    y = frac_molar,
    fill = class
  ) +
  geom_text(
    size = 7/.pt,
    aes(label = n),
    y = 0.5, check_overlap = TRUE
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    #title = "PL class species means",
    x = "Species",
    y = "Mole fraction of phospholipids",
    fill = "Phospholipid class"
  ) +
  coord_flip() +
  theme_tiny() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  guides(fill = guide_legend(nrow = 1, keywidth = 0.5, keyheight = 0.5))

panel_3a
#ggsave(here("04-mainfigs", "Fig3_pl_spmeans_orig5_20230727b.pdf"), width = 60, height = 60, units = "mm")
#ggsave(here("04-mainfigs", "Fig3_pl_spmeans_orig5_20230727.pdf"), width = 6, height = 6)

# Acyl chains for (P)PEs and (P)PCs, all combined
classes_combine = c("PE", "PPE", "PC", "PPC")
panel_3bc_lumped = pldata_wildwhole %>% 
  # clean things up a little
  select(sp, eid, class, id, frac_molar, carbon, dbonds, contains("_col")) %>% 
  filter(sp %in% sp_subset) %>% 
  mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  # filter to a subset of PLs
  filter(class %in% classes_combine) %>% 
  # sum all compounds in each carbon, dbonds group
  group_by(sp, eid, carbon, dbonds) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # renorm each individual. This works; each adds to 1 here.
  group_by(sp, eid) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% #summarize(frac_molar = sum(frac_molar))
  # average individuals within each sp
  group_by(sp, carbon, dbonds) %>% 
  summarize(
    frac_molar = mean(frac_molar),
    class = classes_combine %>% paste(collapse = '+')
  ) %>% 
  # melt it so we can plot carbon + dbonds together
  pivot_longer(cols = c("carbon", "dbonds"), names_to = "predvar", values_to = "pred") %>% 
  # clear the plot of dividing lines
  group_by(sp, predvar, pred, class) %>% 
  summarise(
    frac_molar = sum(frac_molar),
    # special variable to fade odd chains
    fade = ifelse((predvar == "carbon") & (pred %% 2), 0.5, 1)[[1]]
  ) %>% 
  gg_acylch(
    x = pred, 
    y = frac_molar, 
    fill = class,
    alpha = fade,
    cols = vars(predvar),
    meanline = 0.5/.pt
  ) +
    scale_fill_manual(values = c("grey50") %>% setNames(paste(classes_combine, collapse = '+'))) +
    scale_alpha_continuous(limits = c(0, 1)) +
    theme_tiny() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    scale_x_continuous(breaks = c(seq(28, 44, 8), seq(0, 12, 4))) +
    scale_y_continuous(breaks = c(0.0, 0.2, 0.4)) +
    labs(
      #title = "Acyl chain properties for main PL classes",
      x = "Total for both chains",
      y = str_glue({"Mole fraction {classes_combine %>% paste(collapse = '+')}"})
    )
panel_3bc_lumped
ggsave(here("04-mainfigs", "Fig3_acyldists_lumped_orig5_20230727b.pdf"), width = 80, height = 60, units = "mm")
#ggsave(here("04-mainfigs", "Fig3_acyldists_lumped_orig5_20230727a.pdf"), width = 6, height = 6)

# Acyl chains for (P)PEs and (P)PCs, all combined
#classes_display = c("PE", "PPE", "PC", "PPC", "PS")
#classes_display = classes_combine
classes_display = c("PE", "PPE", "PC")
panel_3bc_separate = pldata_wildwhole %>% 
  # clean things up a little
  select(sp, eid, class, id, frac_molar, carbon, dbonds, contains("_col")) %>% 
  filter(sp %in% sp_subset) %>% 
  mutate(sp = sp %>% factor(levels = sp_subset)) %>% 
  # filter to a subset of PLs
  filter(class %in% classes_display) %>% 
  # sum all compounds in each carbon, dbonds group
  group_by(sp, eid, class, carbon, dbonds) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # renorm each individual. This works; each adds to 1 here.
  group_by(sp, eid, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% #summarize(frac_molar = sum(frac_molar))
  # average individuals within each sp
  group_by(sp, class, carbon, dbonds) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # melt it so we can plot carbon + dbonds together
  pivot_longer(cols = c("carbon", "dbonds"), names_to = "predvar", values_to = "pred") %>% 
  # clear the plot of dividing lines
  group_by(sp, predvar, pred, class) %>% 
  summarise(
    frac_molar = sum(frac_molar),
    # special variable to fade odd chains
    fade = ifelse((predvar == "carbon") & (pred %% 2), 0.5, 1)[[1]]
  ) %>% 
  # combine predvar and class into a single factor for faceting
  arrange(predvar, class) %>% 
  mutate(predvar_class = paste(predvar, class) %>% factor(levels = unique(.))) %>% 
  gg_acylch(
    x = pred, 
    y = frac_molar, 
    fill = class,
    alpha = fade,
    cols = vars(predvar_class),
    meanline = 0.5/.pt
  ) +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha_continuous(limits = c(0, 1)) +
  theme_tiny() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_x_continuous(breaks = c(seq(28, 44, 8), seq(0, 12, 4))) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4)) +
  #theme_pubr() +
  labs(
    #title = "Acyl chain properties for main PL classes",
    x = "Total for both chains",
    y = str_glue("Mole fraction of PL class")
  )
panel_3bc_separate
ggsave(here("04-mainfigs", "Fig3_acyldists_separate_orig5_20230727b.pdf"), width = 80, height = 60, units = "mm")
#ggsave(here("04-mainfigs", "Fig3_acyldists_separate_orig5_20230727a.pdf"), width = 6, height = 6)

# Comparison of Arctic Boli and Beroe
# PL species means summed by class
panel_s2e = pldata_wildwhole %>% 
  filter(sp %in% c("Boli_infu", "Bero_cucu")) %>% 
  mutate(sp = sp %>% factor(levels = c("Boli_infu", "Bero_cucu"))) %>% 
  #group_by(eid, sp, class) %>% 
  #summarise(frac_molar = sum(frac_molar)) %>% 
  group_by(sp, class, annot) %>% 
  summarise(
    frac_molar = mean(frac_molar),
    n = {str_glue("N = {n()}")}
  ) %>% 
  gg_headgp(
    x = sp,
    y = frac_molar,
    fill = class
  ) +
  #geom_text(
  #  size = 7/.pt,
  #  aes(label = n),
  #  y = 0.5, check_overlap = TRUE
  #) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    #title = "PL class species means",
    x = "Species",
    y = "Mole fraction",
    fill = "PL class"
  ) +
  #coord_flip() +
  theme_tiny() +
  theme(
    #legend.position = "right",
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.5, keyheight = 0.5))

panel_s2e
ggsave(here("04-suppfigs", "FigS2_acylchains", "panel_s2e_20231006a.pdf"), width = 32.5, height = 60, units = "mm")

# Comparison of Boli and Batho
# PL species means summed by class
panel_s8c = pldata_wildwhole %>% 
  filter(sp %in% c("Boli_micr", "Bath_fost")) %>% 
  mutate(sp = sp %>% factor(levels = c("Boli_micr", "Bath_fost"))) %>% 
  #group_by(eid, sp, class) %>% 
  #summarise(frac_molar = sum(frac_molar)) %>% 
  group_by(sp, class, annot) %>% 
  summarise(
    frac_molar = mean(frac_molar),
    n = {str_glue("N = {n()}")}
  ) %>% 
  gg_headgp(
    x = sp,
    y = frac_molar,
    fill = class
  ) +
  #geom_text(
  #  size = 7/.pt,
  #  aes(label = n),
  #  y = 0.5, check_overlap = TRUE
  #) +
  theme_tiny() +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    #title = "PL class species means",
    x = "Species",
    y = "Mole fraction",
    fill = "PL class"
  ) +
  #coord_flip() +
  theme(
    #legend.position = "right",
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.5, keyheight = 0.5))

panel_s8c
ggsave(here("04-suppfigs", "FigS8_ConfocalLaurdan", "panel_s8c_20231010a.pdf"), width = 32.5, height = 60, units = "mm")

# sn-1 acyl unsat for the above
panel_s8b = pldata_wildwhole %>% 
  filter(sp %in% c("Boli_micr", "Bath_fost")) %>% 
  mutate(sp = sp %>% factor(levels = c("Boli_micr", "Bath_fost"))) %>% 
  #filter(!is.na(dbonsn1)) %>% 
  # sum all classes
  group_by(sp, eid, dbonsn1) %>% 
  summarize(frac_molar = sum(frac_molar)) %>% 
  # renorm
  group_by(sp, eid) %>% 
  # just sat and mono?
  filter(dbonsn1 <= 1) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # average indls
  group_by(sp, dbonsn1) %>% 
  mutate(mean_frac_molar = ifelse(row_number() == 1, mean(frac_molar), NA)) %>% 
  # average unsat
  group_by(sp) %>% 
  mutate(wt_mean = ifelse(row_number() == 1, weighted.mean(dbonsn1, frac_molar), NA)) %>% 
  ggplot(
    aes(x = dbonsn1)
  ) +
  geom_col(aes(
    y = mean_frac_molar,
    fill = sp
  )) +
  geom_point(aes(y = frac_molar)) +
  #geom_vline(aes(xintercept = wt_mean), color = "black") +
  facet_wrap(~sp) +
  scale_fill_manual(values = chroma_sp) +
  theme_tiny() +
  labs(
    x = "sn-1 double bonds",
    y = "Mole fraction of sn-1 chains"
  ) +
  scale_x_continuous(breaks = c(0,1)) +
  guides(fill = "none")
panel_s8b
ggsave(here("04-suppfigs", "FigS8_ConfocalLaurdan", "panel_s8b_20231011e.pdf"), width = 30, height = 60, units = "mm")


# Plot of [PUFA LPE] vs. [PPE] for everyone!
# Plot individual-wise
# OLS only!
# This assesses PPE oxidation
pldata_wildwhole %>% 
  mutate(class = as.character(class)) %>% 
  filter(class %in% c("PPE", "LPE", "PE")) %>% 
  # only PUFA LPEs
  #filter(!((class == "LPE") & (dbonds <= 1))) %>% 
  # annotate PUFA LPEs
  mutate(class = ifelse(
    class != "LPE",
    class,
    ifelse(
      dbonds <= 1, 
      "LPE S/MUFA", 
      "LPE PUFA"
    )
  )) %>% 
  group_by(eid, sp, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  group_by(sp, class) %>% 
  pivot_wider(names_from = "class", values_from = "frac_molar") %>% 
  #summarize(
  #  n = n(),
  #  PE_serr = sd(PE, na.rm = TRUE)/sqrt(n),
  #  PE = mean(PE, na.rm = TRUE),
  #  PPE_serr = sd(PPE)/n,
  #  PPE = mean(PPE)
  #) %>% 
  ggplot(
    #aes(
    #  x = PPE,
    #  y = `LPE PUFA`
    #)
    aes(
      y = `LPE PUFA`/(`LPE S/MUFA`+`LPE PUFA`),
      #y = `LPE PUFA`,
      x = PPE
    )
  ) +
  #geom_errorbar(
  #  size = 0.5/.pt,
  #  color = "grey25",
  #  aes(
  #    ymin = PE - PE_serr,
  #    ymax = PE + PE_serr
  #  )
  #) +
  #geom_abline(color = "white") +
  #geom_errorbarh(
  #  size = 0.5/.pt,
  #  color = "grey25",
  #  aes(
  #    xmin = PPE - PPE_serr,
  #    xmax = PPE + PPE_serr
  #  )
  #) +
  geom_smooth(color = "black", method = "lm", size = 0.5) +
  #geom_text(aes(label = eid)) +
  geom_point(
    size = 1, #TEST
    shape = 21,
    color = "white",
    aes(
      #size = n,
      fill = sp
    )
  ) +
  scale_size_continuous(range = c(2, 6)) +
  scale_fill_manual(values = chroma_sp) +
  theme_tiny() +
  theme(
    legend.position = c(0.7, 0.9),
    axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)
  ) +
  #coord_flip() +
  guides(
    fill = "none",
    size = guide_legend(
      override.aes = list(fill = "black"), 
      keyheight = 0.5,
      keywidth = 0.5,
      nrow = 1
    )
  ) +
  xlim(0, NA) +
  ylim(0, NA)
ggsave(here("04-suppfigs", "FigS2_acylchains", "panel_s2b_SMUFA-PUFA_LPE_20231006d.pdf"), width = 30, height = 30, units = "mm")

# Rsqs to report in the caption for those panel S2B regressions
pldata_wildwhole %>% 
  mutate(class = as.character(class)) %>% 
  filter(class %in% c("PPE", "LPE", "PE")) %>% 
  # only PUFA LPEs
  #filter(!((class == "LPE") & (dbonds <= 1))) %>% 
  # annotate PUFA LPEs
  mutate(class = ifelse(
    class != "LPE",
    class,
    ifelse(
      dbonds <= 1, 
      "LPE S/MUFA", 
      "LPE PUFA"
    )
  )) %>% 
  group_by(eid, sp, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  group_by(sp, class) %>% 
  pivot_wider(names_from = "class", values_from = "frac_molar") %>% 
  ungroup() %>% 
  summarize(
    linreg_upper = lm(`LPE PUFA`~`PPE`) %>% summary() %>% list(),
    linreg_lower = lm((`LPE PUFA`/(`LPE S/MUFA`+`LPE PUFA`))~`PPE`) %>% summary() %>% list()
  ) %>% 
  pivot_longer(cols = contains("linreg_")) %>% 
  .$value

# Plot of [PPE] vs. [PE] for everyone!
pldata_wildwhole %>% 
  group_by(eid, sp, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  filter(class %in% c("PE", "PPE")) %>% 
  group_by(sp, class) %>% 
  pivot_wider(names_from = "class", values_from = "frac_molar") %>% 
  summarize(
    n = n(),
    PE_serr = sd(PE, na.rm = TRUE)/sqrt(n),
    PE = mean(PE, na.rm = TRUE),
    PPE_serr = sd(PPE)/n,
    PPE = mean(PPE)
  ) %>% 
  ggplot(
    aes(
      x = PPE,
      y = PE
    )
  ) +
  geom_errorbar(
    size = 0.5/.pt,
    color = "grey75",
    aes(
      ymin = PE - PE_serr,
      ymax = PE + PE_serr
    )
  ) +
  geom_abline(color = "white") +
  geom_errorbarh(
    size = 0.5/.pt,
    color = "grey75",
    aes(
      xmin = PPE - PPE_serr,
      xmax = PPE + PPE_serr
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
  theme_pubk() +
  theme(legend.position = c(0.7, 0.9)) +
  coord_flip() +
  guides(
    fill = "none",
    size = guide_legend(
      override.aes = list(fill = "black"), 
      keyheight = 0.5,
      keywidth = 0.5,
      nrow = 1
    )
  ) +
  xlim(0, NA) +
  ylim(0, NA) +
  labs(
    x = "[PPE]",
    y = "[PE]"
  )
ggsave(here("04-slidefigs", "PEvsPPE_20230924b.pdf"), width = 3.25, height = 3, units = "in")

# Plot species lipidomes vs. MD systems!
#read_tsv(here("01-rawdata", "complex_md_systems.tsv")) %>% 
read_tsv(here("01-rawdata", "mdsystems_complex_revised_20231223.tsv")) %>% 
  mutate(
    #frac_molar = frac_molar_rnd,
    type = "sim"
  ) %>% 
  factorize_lipids() %>% 
  # bind to averaged and normalized PL + ST data
  bind_rows(plstsp %>% mutate(type = "act")) %>% 
  filter(
    sp %in% c(
      "Boli_vitr",
      "Boli_infu",
      "Lamp_crue",
      "Tjal_pink",
      "Boli_micr",
      "Bath_fost"
    )
  ) %>% 
  mutate(
    sp = factor(sp, levels = c(
      "Boli_vitr",
      "Boli_infu",
      "Lamp_crue",
      "Tjal_pink",
      "Boli_micr",
      "Bath_fost"
    )),
    type_sp = paste(type, sp)
  ) %>% 
  gg_headgp(
    # hack to make normalization work right
    x = type_sp,
    y = frac_molar,
    fill = class,
    darkmode = FALSE
  ) +
  facet_wrap(~sp, nrow=1, scales = "free_x") +
  #scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    #title = "PL class species means",
    x = "",
    y = "Mole fraction",
    fill = "Lipid\nclass"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.5, keyheight = 0.5))
#ggsave(here("04-slidefigs", "lipidomeVsMD_20230928a.pdf"), width = 6, height = 3, units = "in")
ggsave(here("04-suppfigs", "FigS7_ComplexSims", "PanelS7a_revisedsys_20231223a.pdf"), width = 180, height = 60, units = "mm")

# Plot ALL E. coli lipidomes
ecoli_testdata = pldata %>% 
  filter(sp == "Esch_coli") %>%
  # FILTER AND NORM OUT PPC
  group_by(eid) %>% 
  filter(!(class %in% c("PPC", "LPC"))) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar))# %>% 
  gg_headgp(
    x = eid,
    y = frac_molar,
    fill = class
  ) +
  facet_wrap(~sid, scales = "free_x", nrow=1) +
  theme_tiny() +
  theme(
    legend.position = "right"
  ) +
  # change for number of facets
  scale_x_discrete(labels = rep(c(1, 2, 3), 8)) +
  labs(
    x = "Culture replicate",
    y = "Mole fraction"
  )
ggsave(here("04-omitfigs", "all_Ecoli_lipids_20230930a.pdf"), width = 180, height = 60, unit = "mm")

# Plot just PE/PPE, PE/PPC, averaged
pldata %>% 
  filter(sp == "Esch_coli") %>%
  # FILTER AND NORM OUT PPC
  group_by(eid) %>% 
  filter(!(class %in% c("PPC", "LPC"))) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  filter(depth_col == 0) %>% 
  filter(
    str_detect(sid,  "AL") |
      str_detect(sid, "pET28") |
      str_detect(sid, "pPLs")
  ) %>% 
  # average across replicates
  group_by(sid, class, annot) %>% 
  summarise(frac_molar = mean(frac_molar)) %>% 
  group_by(sid) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # sum classes
  group_by(sid, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  # order the strains
  mutate(sid = factor(sid,
      levels = c(
        "BL21 + pET28a 1 bar",
        "AAL9256 (PE)", 
        "BL21 + pPLsCP 1 bar",
        "AAL95 + pAC (PC)"
      )
  )) %>% 
  gg_headgp(
    x = sid,
    y = frac_molar,
    fill = class
  ) +
  facet_wrap(~sid, scales = "free_y", nrow=2) +
  theme_tiny() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_flip() +
  labs(
    x = '',
    y = "Mole fraction"
  )
ggsave(here("04-suppfigs", "PanelS9d_micro_means.pdf"), width = 60, height = 30, unit = "mm")

# Panel S9A: E. coli headgroup composition by strain
pldata_ecoli %>% 
  # WHY OH WHY are the classes out of order?!?
  gg_headgp(
    x = sid, 
    y = frac_molar,
    fill = class
  ) +
  facet_wrap(~sid, scales = "free_x", nrow=1) +
  theme_tiny() +
  theme(
    #legend.position = "right",
    panel.spacing.x = unit(6, "mm"),
    legend.position = "none"
  ) +
  labs(
    y = "Mole fraction"
  )
ggsave(here("04-suppfigs", "FigS9_Ecoli", "Panel_s9a_20231008c.pdf"), width = 120, height = 60, unit = "mm")

# Panel S9B: E. coli acyl chain CLI, DBI by strain
pldata_ecoli_indls %>% 
  # only species with acyl chains resolved
  filter(str_detect(annot, regex("_|\\/"))) %>% 
  # put them all in 1 col
  pivot_longer(cols = matches(regex("sn[12]")), names_to = "varstereo", values_to = "val") %>% 
  # deidentify the chains
  mutate(varmono = str_remove(varstereo, regex("sn[12]"))) %>% 
  # strip out "trace" fatty acids that E. coli definitely cannot make
  filter(
    ((varmono == "carb") & between(val, 12, 19)) |
      ((varmono == "dbon") & between(val, 0, 1))
  ) %>% 
  # sum the classes
  group_by(sp, eid, sid, varmono, val) %>% 
  summarise(
    frac_molar = sum(frac_molar),
    class = "all"
  ) %>% 
  # normalize to total acyl chains
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  group_by(sp, sid, varmono, val) %>%
  mutate(mean_frac_molar = ifelse(row_number() == 1, mean(frac_molar), NA)) %>% 
  gg_acylch(
    x = val, 
    y = mean_frac_molar,
    meanline = FALSE
  ) +
  geom_point(
    aes(y = frac_molar),
    size = 0.5,
    position = position_jitter(width = 0.2, height = 0)
  ) +
  facet_grid(
    rows = vars(varmono),
    cols = vars(sid),
    scales = "free_y"
  ) +
  coord_flip() +
  theme_tiny() +
  theme(
    #legend.position = "right",
    #panel.spacing.x = unit(6, "mm"),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = seq(0, 19)) +
  labs(
    y = "Mole fraction of acyl chains"
  )
ggsave(here("04-suppfigs", "FigS9_Ecoli", "Panel_s9b_acylch_20231009a.pdf"), width = 120, height = 60, unit = "mm")

# Body 
pldata %>% 
  filter(sp == "Coel_hali") %>% 
  mutate(
    sp_eid = paste(sp, eid),
    sp_eid = sp_eid %>% factor(levels = unique(sp_eid))
  ) %>% 
  # comment to show individual compounds
  #group_by(eid, sp, class) %>% 
  #summarise(frac_molar = sum(frac_molar)) %>% 
  gg_headgp(
    x = sp_eid,
    y = frac_molar,
    fill = class
  ) +
  facet_wrap(~tissue, scales="free_x") +
  theme_tiny() +
  labs(
    title = "PL class compositions",
    x = "Sample",
    y = "Mole fraction of phospholipids"
  )
ggsave(here("04-suppfigs", "Coel_hawaii_bodyvtent_20231010a.pdf"), width = 90, height = 120, unit = "mm")

# Do a little reproducibility analysis on the Leuco tech replicates
data_techrep = pldata %>% 
  filter(str_detect(eid, "JWL0366")) %>% 
  mutate(rep = eid %>% as.factor() %>% as.integer())

# plot side-by-side
data_techrep %>% 
  gg_headgp(
    x = rep,
    y = frac_molar,
    fill = class
  ) +
  theme_tiny() +
  theme(legend.position = "right") +
  labs(
    x = "Technical replicate",
    y = "Mole fraction phospholipids"
  )

# get standard errors
techrep_errors = data_techrep %>% 
  group_by(class, annot) %>% 
  summarize(
    sem = sd(frac_molar)/sqrt(n()),
    mean = mean(frac_molar),
    sem_frac_of_mean = sem/mean
  ) %>% 
  filter(mean > 0) %>% 
  arrange(-sem)

techrep_errors %>% 
  ggplot(
    aes(
      x = mean,
      y = sem_frac_of_mean
    ) 
  ) +
  geom_point()

# what's the largest SEM?
techrep_errors$sem %>% summary()

pldata %>% 
  filter(str_detect(eid, regex("JWL0366[AB]"))) %>%
  group_by(eid) %>% 
  pivot_wider(names_from = "eid", values_from = "frac_molar") %>% 
  ggplot(
    aes(x = JWL0366A, y = JWL0366B)
  ) +
  geom_point()



