# plotting helper for the "cumulative stacked barplot"
# used to visualize curvature contributions

source(here("03-scripts", "seelipids_helpers.R")) # Some color palettes riff on these

# cteno species color palette!
# ooh aah
chroma_sp = c(
  "Leuc_pulc" = "#FF7F00", # orange
  "Bero_abys" = "#08529C",
  "Bero_cucu" = "#3A7DB4",
  "Bero_pseu" = "#3A7DB4",
  "Bero_ovat" = "#9ECAE1",
  "Boli_vitr" = "#984EA3", # the original purple
  "Boli_infu" = "#109292",
  "Boli_micr" = "#109292", # muted teal used in 2021 paper
  "Cest_vene" = "#68C3A4",
  "Coel_hali" = "#0D4420", # Halimeda green?
  "Lamp_crue" = "#E41A1C",
  "Bath_fost" = "#A65628", # brick red/brn
  "Llyr_deep" = "#FCC06F", # peach gold
  "Llyr_bent" = "#FCC06F", 
  "Cydi_blac" = "#000000", # straight up black
  "Tjal_pink" = "#F781BF",
  #"other"     = "#999999"
  "Mert_angu" = "#8A8B8A" # 50% grey bc related to black cyd
)

# color mapping for E. coli samples
chroma_pheno = c(
  "PPE-" = chroma_cl[["PE"]],
  "PPE+" = chroma_cl[["PPE"]],
  "PE"   = chroma_cl[["PE"]],
  "PC"   = chroma_cl[["PC"]]
)

# muted palette for lipid phases
# these are from calecopal: cal_palette("bigsur")[order(c(3,1,2))]
chroma_ph = c(
  hii = "#ECBD95", 
  la  = "#E4DECE",
  lb  = "#9BB1BB"
)

# custom "cumulative stacked barplot" wrapper
# for visualizing curvature contributions
gg_plcurv = function(
    data,
    darkmode = FALSE,
    thres_draw = 0, # can pass a mole fraction threshold below which color blocks are removed.
    label_frac = 0.015, # min fraction at which bar get labeled; set to 1 to disable labels.
    label_size = 1.5,
    baseline = 0.5/.pt, # *weight* of the baseline
    # aesthetics get passed in as naked args
    # requires `cols = [sp|eid]`
    ...
){
  # unpack the ellipsis args as strings
  mapstrs = lapply(rlang::enexprs(...), as.character)
  # set up row faceting if desired
  if("rows" %in% names(mapstrs)){
    rowvar = vars(eval(sym(mapstrs$rows)))
    rownam = sym(mapstrs$rows)
  }else{
    rowvar = NULL
    rownam = NULL
  }
  # ensure ordering
  data = data %>% arrange(class) # and maybe even id?
  this_gg = data %>%
    # apply threshold
    filter(abs(eval(sym(mapstrs$y))) >= thres_draw) %>% 
    # group by the passed x aesthetic
    #group_by(eval(sym(mapstrs$cols))) %>%
    group_by(eval(sym(mapstrs$cols)), eval(rownam)) %>%
    # sum or mean?
    #summarize(plci = sum(plci)) %>% # try leaving this for the user to do upstream
    mutate(neg = (plci<0)*0.3) %>% # 0.3 is the offset width
    arrange(neg, class) %>% 
    mutate(
      # running total
      end = cumsum(plci),
      y = end - plci/2
    ) %>% 
    ggplot(aes(x = neg)) +
    # due to plotting mechanics, columns are mapped as cols, not x
    facet_grid(
      rows = rowvar,
      cols = vars(eval(sym(mapstrs$cols))), 
      switch = 'x'
    ) +
    geom_hline(
      size = baseline,
      yintercept = 0, color = ifelse(darkmode, "white", "black")
    ) +
    geom_tile(
      aes(
        y = y,
        height = plci,
        fill = class
      ),
      width = 0.25,
      size = 0.05,
      color = ifelse(darkmode, "black", "white")
    ) +
    scale_fill_manual(values = chroma_cl) +
    ## ind'lwise error bars
    #geom_errorbar(
    #  data = data %>% 
    #    group_by(sp, sp_eid, class) %>%
    #    summarize(frac_molar = sum(frac_molar)) %>%
    #    full_join(icurv, by = c("class")) %>% 
    #    mutate(plci = frac_molar * c0) %>% 
    #    drop_na() %>% 
    #    group_by(sp, sp_eid, scenario) %>% 
    #    summarize(plci = sum(plci)) %>% 
    #    # average by species
    #    group_by(scenario, sp) %>% 
    #    summarize(
    #      plci_mean = mean(plci),
    #      plci_serr = sd(plci)/sqrt(n()),
    #      neg = TRUE * 0.3
    #    ),
    #  aes(
    #    ymin = plci_mean - plci_serr,
    #    ymax = plci_mean + plci_serr
    #  ),
    #  width = 0.05
    #) +
    theme_pubr() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    guides(color = "none") +
    labs(y = "Mean c0 (1/Å)")
  this_gg
}

# compact theme following AAAS figure prep guidelines
theme_tiny = function(...){
  theme_pubr(base_size = 7, ...) +
    theme(
      # axis options
      axis.line  = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      strip.background = element_rect(linewidth = 0),
      # bring the legend in tight
      legend.margin = margin(rep(0,4)),
      plot.margin = margin(c(0,2,0,2))
    )
}

# Use ggtree to plot a phylogeny with color-coded tips
# pretty basic atm, leans on ggtree
gg_lollitree = function(
    phylo,
    darkmode = FALSE,
    # aesthetics get passed in as naked args
    # requires `cols = [sp|eid]`
    linewidth = 0.5,
    ...
){
  phylo %>% 
    ggtree(
      mapping = aes(...),
      size = linewidth,
      ladderize = FALSE
    ) +
      geom_tippoint(
        shape = 21,
        size = 2,
        color = "white",
        aes(
          # per https://guangchuangyu.github.io/2015/09/subsetting-data-in-ggtree/
          fill = label, # label and isTip are internally recognized!
          subset=isTip 
        )
      ) +
      scale_fill_manual(values = chroma_sp)
}

# utility function to get the tip names from a phylo object, in plotting order!
tipOrder = function(phylo){
  tiplabs = phylo$tip.label
  tipnums = phylo$edge %>% .[,2] 
  tiplabs[tipnums %>% .[which(.<=length(tiplabs))]]
}