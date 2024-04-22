## Scatter/bubble-plot lipidomes' mean curvature under multiple estimation schemes
## Barplots and lipid contributions to curvature are in `plot_curvature.R`

library(tidyverse)
library(here)
library(ggpubr)

source(here("03-scripts", "prep_lipidmaps.R"))
source(here("03-scripts", "plot_helpers.R"))
source(here("03-scripts", "pgls_helpers.R"))

# The phylogeny we're gonna use
#NTS 20230722: `20230722_iq.tre` is a hacked-up little tree with the necessary
# short names subbed in manually. Regenerate.
file_tree = here("01-rawdata", "20230722_iq.tre")
# load the phylogeny so we can fix anything that's wrong
phylodata = ape::read.tree(file_tree)
plot(phylodata, main = "Input tree") # have a quick look

# summarize some stuff ind'lwise
# we want total frac_molar, sn1, sn2 chain and dbi within each class
pldata_indls = pldata_wildwhole %>% 
  # bind sterol data, but leave PL normalization for everything else
  bind_rows(
    plstdata %>% 
      filter(class == "ST") %>% 
      # JWL0185 is a tentacle sample with sterol data!
      filter(wild & (tissue %in% c("whole", "body")))
  ) %>% 
  bind_rows(
    plsmdata %>% 
      filter(class == "SM") %>% 
      # JWL0185 is a tentacle sample with sterol data!
      filter(wild & (tissue %in% c("whole", "body")))
  ) %>% 
  # bind SM data, but leave PL normalization for everything else
  # EXCLUDE TECH REPLICATES
  filter(!str_detect(eid, regex("[A-Z]$"))) %>% 
  group_by(sp, eid, depth_col, temp_col, class) %>%
  summarize(
    # weighted means
    across(contains(c("carb", "chn", "db")), function(x){stats::weighted.mean(x, frac_molar, na.rm = TRUE)}),
    # class totals
    frac_molar = sum(frac_molar)
  )

# add "all" as a class (for acyl chain properties)
# I tried to do this in the pipeline but magrittr was being weird
pldata_indls_all = bind_rows(
    pldata_indls,
    pldata_indls %>% 
      # same summarize call as above!
      summarize(
        # weighted means
        across(contains(c("carb", "chn", "db")), function(x){stats::weighted.mean(x, frac_molar, na.rm = TRUE)}),
        # class totals
        frac_molar = sum(frac_molar),
        class = "all"
      )
)  %>% 
  mutate(class = class %>% factor(levels = c("all", levels(pldata_indls$class))))

# an expanded version to test conditional exclusion
# or species jackknifing
pldata_indls_opt = pldata_indls_all %>% 
  # define expansion axes here
  cross_join(crossing(
    # try filtering to cold samples, or not
    subset = c("all samples", "<= 10 deg C", "<= 250 m"),
    # jackknife one species at a time!
    # Don't forget we still want an all-inclusive dataset, hence "nothing".
    jknife = c("nothing", pldata_indls$sp %>% unique()) %>% factor(levels = .)
  )) %>% 
  # define conditions for expansion here
  # the order (and prob separate filter calls) is important
  filter(!((subset == "<= 10 deg C") & (temp_col > 10))) %>% 
  filter(!((subset == "<= 250 m")  & (depth_col > 250))) %>% 
  filter(as.character(sp) != as.character(jknife)) %>% 
  # melt the predictor columns
  # it's important to be able to run multiple predictors
  pivot_longer(cols = contains("_col"), names_to = "predvar", values_to = "pred") %>% 
  # and finally, melt the response variables
  pivot_longer(cols = contains(c("carb", "chn", "db", "frac_molar")), names_to = "respvar", values_to = "resp")

### OK, here we go w/PGLS...

### This block preps the tree to work with any subset of pheno data
# Reduce the phylo- and phenodata to the intersection of their species
# give warnings if species in phenodata are not found in phylogeny
phenophylo = match_tree(
  pldata_indls, # don't need to use opt here, takes up extra space
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
# this takes a LONG TIME with all the classes and species jackknifing
# 720 regressions (groups*2) even if I leave out acyl chains and temperature
# 864 models if I include chain info!
# candidate for parallelization? Or just filter to what I care about? Or both?
allmods_lmap = pldata_indls_opt %>% 
  #filter((class=="SM")&(respvar=="frac_molar")&(jknife=="nothing")&(subset!="all samples")) %>%  #TEST
  # lighten the load
  filter(
    # omit jackknives!
    (jknife == "nothing") &
    (
      # for headgroup trends
      ((respvar == "frac_molar") & (class != "all")) |
        # for overall chain trends
        ((respvar %in% c("chn", "dbi", "carbsn1", "carbsn2", "dbonsn1", "dbonsn2")) & (class == "all"))
    ) &
      # 1224 groups w/o, 816 with => 33.33% savings
      !((predvar == "depth_col") & (subset == "<= 250 m")) &
      !((predvar == "temp_col") & (subset == "<= 10 deg C"))
  ) %>% 
  group_by(subset, class, predvar, respvar, jknife) %>% 
  # can specify all the models I want to fit here!
  # to fit multiple models, need to stuff data into static listcol
  summarize(thesedata = cur_data() %>% list()) %>% 
  rowwise() %>% 
  # go parallel!
  group_split() %>% 
  #.[1:5] %>% #TEST
  pbapply::pblapply(
    cl = 8L, #TEST
    X = .,
    FUN = function(row){
      row %>% 
        mutate(
          # OLS
          modols = safely(lm)(
            resp ~ pred, 
            # need to strip missing vals for accurate regressions
            thesedata[[1]] %>% filter(!is.na(pred) & !is.na(resp))
          )$result %>% list(),
          # GLS Brownian?
          # GLS Martins
          modgls = safely(pgls)(
            resp ~ pred, 
            # need to strip missing vals for accurate regressions
            thesedata[[1]] %>% filter(!is.na(pred) & !is.na(resp)), 
            ape::corMartins, 
            phylo_ult,
            #gammas = seq(-5, 15) %>% 10**. # default
            # fixes the SM problem and theoretically saves time; does it mess up others?
            gammas = seq(1, 15) %>% 10**., # spot-checked and it looks good! (same for all others) - 20231012
            debug = TRUE #TEST, needs to be put in single-thread to print
        )$result %>% list()
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
  # regroup?
  group_by(subset, class, predvar, respvar, jknife, model)

saveRDS(allmods_lmap, here("02-tidydata", "allmods_lmap.rds"))

# make a copy with adjusted p-vals
allmods_multtest = allmods_lmap %>% 
  #NTS: define a dedicated grouping variable to define overlapping tests
  # note that "class" is not in here! (that's the multiple-testing axis)
  group_by(subset, predvar, respvar, jknife, model, term) %>% 
  mutate(testgroup_headgp = ifelse(respvar == "frac_molar", cur_group_id(), 0)) %>% 
  # and "respvar" is not in here!
  group_by(subset, class, predvar, jknife, model, term) %>% 
  mutate(testgroup_acylch = ifelse(respvar %in% c("chn", "dbi"), cur_group_id(), 0)) %>% 
  # combine the two group indices (i.e. Cartesian product) %>% 
  mutate(testgroup = paste(testgroup_headgp, testgroup_acylch)) %>% 
  group_by(testgroup) %>% 
  # tidyselect-fu to make a column for every correction method
  # much more memory-efficient, don't think I need it melted
  mutate(
    across(
      .cols = matches("p.value"),
      .fns  = c(
        "holm", 
        "hochberg", 
        "hommel", 
        "bonferroni", 
        "BH", 
        "BY", 
        "fdr"
      ) %>% setNames(., .) %>% 
        lapply(., function(meth){function(ps){p.adjust(ps, meth)}})
    )
  )
  
# This is really for plotting, but it's here cuz it's kinda slow.
# generate prediction lines based on all the models
modelines_lmap = pldata_indls_opt %>% 
  group_by(subset, class, predvar, respvar, jknife) %>% 
  # filter to extreme ends of each predictor (2 points make a line!)
  filter(pred %in% range(pred)) %>% 
  distinct(pred) %>% 
  # join automatically on grouping vars
  right_join(
    allmods_multtest %>% 
      group_by(subset, class, predvar, respvar, jknife, model) %>% 
      summarize(fit = fit[[1]] %>% list()),
    relationship = "many-to-many"
  ) %>% 
  group_by(subset, class, predvar, respvar, model, jknife) %>% 
  mutate(resp = predict(fit[[1]], cur_data()))
