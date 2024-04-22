## Wrapper functions and other objects to ease the fitting of
## Phylogenetically Normalized Least-Squares (PGLS) models
##
## Assumptions chosen for this paper:
## (1) Ultrametrization by NNLS makes the tree closer to the true time tree.
## (2) A Martins-style "burst" model of trait evolution is more realistic than 
##     Brownian motion (though both are tried).
## (3) Conspecifics identified morphologically and/or with 18S or COI barcodes
##     share 100% evolutionary background (are represented as polytomies).

library(tidyverse)
library(here)
library(phytools)
library(tidymodels)
library(AICcmodavg) # provides predictSE.gls
#library(broom.mixed)

# match phenotypic data to a phylogeny
# return namedlist containing tibble and phylo
# tip2sp is the function that converts tip names to species
match_tree = function(pheno, phylo, tip2sp=identity){
  
  # species in the pheno dataset
  sp = unique(pheno$sp)
  # warn about any species not in the tree
  sp_missing = sp[which(!(sp %in% tip2sp(phylo$tip.label)))]
  if(length(sp_missing) > 0){message(paste(paste(sp_missing, collapse=", "), "not found in tree!", sep=' '))}
  
  # remove those from the table
  pheno_ok = pheno %>% filter(sp %in% phylo$tip.label)
  
  # prune taxa not in pheno
  phylo_pruned = phylo %>% keep.tip(phylo$tip.label[which(phylo$tip.label %in% pheno_ok$sp)])
  
  # add conspecifics as polytomies, use for loop because node matching recurses
  samps = pheno_ok %>% rowwise() %>% group_split()
  phylo_indl = phylo_pruned
  for(row in samps){
    phylo_indl = phylo_indl %>%
      bind.tip(tip.label = row$sp, where = match(row$sp, phylo_indl$tip.label), edge.length = 0)
  }
  
  # drop leaves that do not match samps verbatim
  phylo_ok = phylo_indl %>% 
    # keep.tip automatically drops redundantly named leaves
    keep.tip(.$tip.label[which(.$tip.label %in% pheno_ok$sp)])
  
  return(list("pheno" = pheno_ok, "phylo" = phylo_ok))
}

# this function takes a phylo object and a table of phenotypic data with column "sp"
# It consolidates duplicate tip labels, then makes polytomies to match the number of indl's/sp in pheno.
expand_tree = function(pheno, phylo, tip2sp=identity){
  # construct branch length table
  branches = cbind(phylo$edge, phylo$edge.length) %>%
    as_tibble() %>%
    set_names(c("mrca", "node", "dist")) %>%
    mutate(
      node = as.integer(node),
      mrca = as.integer(mrca)
    )
  
  # get duplicated taxa and their intraspecific MRCAs
  leaves = tibble(sp = phylo$tip.label) %>%
    mutate(node = row_number())
  
  leaves_dup = leaves %>%
    group_by(sp) %>%
    mutate(count = n()) %>%
    filter(count > 1) %>%
    mutate(mrca = phylo %>% getMRCA(node) %>% as.integer())
  
  # define duplicated taxa and their average distance from MRCA
  taxa_dup = leaves_dup %>%
    rowwise() %>%
    mutate(dist = tibble(node, mrca) %>% left_join(branches, by = c("node", "mrca")) %>% .$dist) %>%
    group_by(sp, mrca) %>%
    quietly(summarize)(
      node = min(node),
      dist = mean(dist, na.rm = TRUE)
    ) %>% .$result %>% 
    ungroup()
  
  # define leaves to drop (all but the first, per sp)
  leaves_drop = leaves_dup %>%
    anti_join(leaves_dup %>% summarize(node = first(node)), by = c("sp", "node"))
  
  # set first conspecific branch to average length
  branches_adj = bind_rows(anti_join(branches, taxa_dup, by = c("node", "mrca")), taxa_dup)
  
  # replace conspecific clades with a single tip
  phylo_pruned = phylo %>%
    compute.brlen(
      phylo$edge %>%
        as_tibble() %>%
        set_names(c("mrca", "node")) %>%
        # put in order
        left_join(branches_adj, by = c("mrca", "node")) %>%
        pull(dist)
    ) %>%
    drop.tip(leaves_drop %>% pull(node)) %>%
    # finally, drop taxa not in phenodata
    drop.tip(phylo$tip.label %>% .[which(!(. %in% pheno$sp))])
  
  # Supplement pheno table for output
  # and identify the var we'll use to name ind'l tips
  tip_id = group_vars(pheno) %>% paste(collapse='_')
  pheno_indl = pheno %>% 
    group_by(sp) %>% 
    #NTS 20230722 would be really nice to soft-code this for production
    mutate(sp_eid = paste(sp, eid, sep=':')) %>% 
    #mutate(across(matches(tip_id), )) %>% 
    group_by(sp) %>% 
    arrange(sp, eid)
    #mutate(indl = paste(sp, row_number(), sep=''))
  
  # Start generating a list of the to-be tip names (indls)
  #leaves_add = pheno %>% 
  #  summarize() %>% 
  #  group_by(sp) %>% 
  #  summarise(indl = paste(sp, seq(n()), sep='') %>% list()) %>% 
  #  unnest(indl) %>% 
  #  group_by(sp) %>% 
  #  arrange(indl) %>% 
  #  filter(row_number() > 1)
  
  # this patch is a remnant, but allows some flexibility
  leaves_add = pheno_indl
  phylo_indl = phylo_pruned
  ## concat "1" to all the existing tip labels
  #phylo_indl$tip.label = paste(phylo_indl$tip.label, '1', sep='')
  
  # iteratively add leaves to tree
  # for loop is used for same reason as before
  edge_polytomy = 0 # length of a tooth on the polytomy comb
  for(row in leaves_add %>% rowwise() %>% group_split()){
    phylo_indl = phylo_indl %>%
      bind.tip(tip.label = row$sp_eid, where = match(row$sp, tip2sp(phylo_indl$tip.label)), edge.length = edge_polytomy)
  }
  
  ## number the individuals in the tree_rooted and the table
  #phylo_indl$tip.label = tibble(sp = phylo_indl$tip.label) %>%
  #  group_by(sp) %>%
  #  mutate(indl = paste(sp, row_number(), sep="")) %>%
  #  pull(indl)
  
  # drop the original species-only (no eid) -labeled tips
  # hack: ones with no whitespace in 'em
  phylo_indl = phylo_indl %>% 
    drop.tip(phylo_indl$tip.label %>% .[which(!str_detect(., ':'))])
  
  return(list("pheno" = pheno_indl, "phylo" = phylo_indl))
}

glscorr = function(mod){
  # get that correlation matrix out of the gls
  # per https://stackoverflow.com/questions/23039638/how-to-retrieve-correlation-matrix-from-glm-models-in-r
  coef(mod$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)
}

# multiply all branches in passed tree by passed value
# and return stretched tree
phystretch = function(phylo, x){
  phylo$edge.length = phylo$edge.length * x
  return(phylo)
}

# tidy the fit and glance the model
# bind them together with nonredundant colnames
tidyglance = function(fit){
  cross_join(
    # the coefficient data take priority
    broom.mixed::tidy(fit),
    broom.mixed::glance(fit),
    suffix = c('', '.glance')
  ) %>% 
    # glance does not give rsq data for a gls, so that gets removed
    select(!contains(".glance"))
}

# Master PGLS wrapper that works just like gls()
# Give it a formula, data, a corfunc, and a phylogeny
# Matching of the phylogeny to the data is done every call,
# to make exploratory analyses easier
# gammas stretch tree branches (important for non-Brownian corr strucs)
# ... is passed to corfunc
pgls = function(formula, data, corfunc, phylo, startval=1, gammas=seq(-5, 15) %>% 10**., debug=FALSE, ...){
  # Match tree to pheno data
  # it is hardcoded here for my variable names,
  # but it would be nice to do it by grouping of pheno
  phenophylo_indl = expand_tree(
    data %>% group_by(sp, eid), 
    phylo_ult#,
    #tip2sp = function(x){str_remove_all(x, "[0-9]")}
  )
  pheno_indl = phenophylo_indl$pheno %>% ungroup()
  phylo_indl = phenophylo_indl$phylo
  
  #print(phylo_indl) #TEST
  #print(pheno_indl) #TEST
  
  # There is an issue with float handling, at least in corMartins,
  # such that unreasonably small or large branch lengths
  # quietly give bogus results.
  # To handle this, I sweep the branch multiplier (which should correspond
  # to gamma in Martins' model), then pick the highest LogLik.
  # branch mults: 1E-5 to 1E+15 by factors of 10. Maybe pass as arg?
  fits = tibble(gamma = gammas) %>%
    rowwise() %>% 
    mutate(
      phylo = phylo_indl %>% phystretch(gamma) %>% list(),
      fit   = safely(nlme::gls)(
        data = pheno_indl,
        model = as.formula(formula),
        correlation = corfunc(startval, phy = phylo, form = ~sp_eid)
      )$result %>% list()
    ) %>% 
    # strip out failed fits
    filter(!is.null(fit)) %>% 
    # spread out the summary stats
    #mutate(fitout = fit %>% glance() %>% list()) %>% 
    # can also tidyglance() if we want better debugging
    mutate(fitout = fit %>% tidyglance() %>% list()) %>% 
    unnest(fitout)
  
  # this option will print a summary tibble for ALL the fits completed
  if(debug){fits %>% print(n=nrow(.))}
  
  # return the fit with max LogLik and gamma nearest 1
  fits %>% 
    ungroup() %>% 
    # change prioritization of fits here
    arrange(-logLik, abs(log10(gamma))) %>% 
    # this has to be done in two steps to unpack properly, idk why
    .$fit %>% .[[1]]
}


### SCRATCH

## Master PGLS wrapper that works just like gls()
## Give it a formula, data, a corfunc, and a phylogeny
## Matching of the phylogeny to the data is done every call,
## to make exploratory analyses easier
## ... is passed to corfunc
#pgls_old = function(formula, data, corfunc, phylo, startval=1, debug=FALSE, ...){
#  # Match tree to pheno data
#  # it is hardcoded here for my variable names,
#  # but it would be nice to do it by grouping of pheno
#  phenophylo_indl = expand_tree(
#    data %>% group_by(sp, eid), 
#    phylo_ult#,
#    #tip2sp = function(x){str_remove_all(x, "[0-9]")}
#  )
#  pheno_indl = phenophylo_indl$pheno %>% ungroup()
#  phylo_indl = phenophylo_indl$phylo
#  
#  #print(phylo_indl) #TEST
#  #print(pheno_indl) #TEST
#  
#  # There is an issue with float handling, at least in corMartins,
#  # such that unreasonably small or large branch lengths
#  # quietly give bogus results.
#  # To handle this, I sweep the branch multiplier (which should correspond
#  # to gamma in Martins' model) until alpha starts changing.
#  #gammas = lapply(seq(0, 20, 1), function(x){10**x}) %>% unlist()
#  init_gamma = 1E-5
#  mult_gamma = 1E+1
#  max_gammas = 20
#  # compress the branch lengths initially
#  phylo_indl$edge.length = phylo_indl$edge.length * init_gamma
#  fits = c()
#  for(i in seq(max_gammas)){ # sets max # iterations
#    # run fit and PREpend to vector
#    # (overwriting object causes malloc/segfault)
#    fits = c(
#      safely(nlme::gls)(
#        data = pheno_indl,
#        model = as.formula(formula),
#        correlation = corfunc(startval, phy = phylo_indl, form = ~sp_eid)
#      )$result %>% list(),
#      fits
#    )
#    # above can append a list(NULL), so strip these
#    fits = fits %>% .[lengths(.) > 0]
#    # conditional to break loop
#    if(
#      !debug &&
#      (length(fits) > 1) && # get past the first fit
#      #(fits[[2]]$apVar[2] != 0) && # last fit not degenerate, there is correl
#      !base::all.equal(fits[[1]]$logLik, fits[[2]]$logLik) && # avoid breaking due to float error
#      (fits[[1]]$logLik < fits[[2]]$logLik) # and lower LL than current
#    ){
#      print("SNAP!")
#      break
#    }
#    print(str_glue({"{pheno_indl$jknife[[1]]} {i}"})) #TEST
#    #print(fits[[1]]) #TEST
#    # stretch the branch lengths
#    phylo_indl$edge.length = phylo_indl$edge.length * mult_gamma
#  }
#  # print ALL the fit info
#  if(debug){
#    rlang::local_options(pillar.sigfig = 7)
#    pheno_indl %>% 
#      group_by(scheme, jknife) %>% 
#      summarize() %>% 
#      cross_join(tibble(fit = fits)) %>% 
#      rowwise() %>% 
#      mutate(fitout = fit %>% tidyglance() %>% list()) %>% 
#      unnest(fitout) %>% 
#      filter(term == "pred") %>% 
#      print(n=100)
#  }
#  # return previous fit
#  # this will exist by the time the loop breaks
#  return(fits[[2]])
#}


## can check which corfunc is passed like so
# if(corfunc %>% substitute() %>% deparse() %>% str_detect("corMartins")){}

# with the old do.call data-preserving implementation
#pgls = function(formula, data, corfunc, phylo, startval=1, ...){
#  ## can check which corfunc is passed like so
#  # if(corfunc %>% substitute() %>% deparse() %>% str_detect("corMartins")){}
#  
#  # Match tree to pheno data
#  # it is hardcoded here for my variable names,
#  # but it would be nice to do it by grouping of pheno
#  phenophylo_indl = expand_tree(
#    data %>% group_by(sp, eid), 
#    phylo_ult#,
#    #tip2sp = function(x){str_remove_all(x, "[0-9]")}
#  )
#  pheno_indl = phenophylo_indl$pheno %>% ungroup()
#  phylo_indl = phenophylo_indl$phylo
#  
#  print(phylo_indl) #TEST
#  print(pheno_indl) #TEST
#  
#  # There is an issue with float handling, at least in corMartins,
#  # such that unreasonably small or large branch lengths
#  # quietly give bogus results.
#  # To handle this, I sweep the branch multiplier (which should correspond
#  # to gamma in Martins' model) until alpha starts changing.
#  #gammas = lapply(seq(0, 20, 1), function(x){10**x}) %>% unlist()
#  init_gamma = 1E-5
#  mult_gamma = 1E+1
#  max_gammas = 20
#  # compress the branch lengths initially
#  phylo_indl$edge.length = phylo_indl$edge.length * init_gamma
#  fits = c()
#  for(i in seq(max_gammas)){ # sets max # iterations
#    # using do.call stores the formula (and data?) in the fit object
#    # how mem-wasteful is this?
#    phylo_indl$edge.length %>% print() #TEST
#    # run a fit and append to vector
#    # (overwriting object causes malloc/segfault)
#    fits = c(fits, 
#             #do.call(
#             safely(nlme::gls)(#,
#               #list(
#               data = pheno_indl,
#               model = as.formula(formula),
#               correlation = corfunc(startval, phy = phylo_indl, form = ~sp_eid)
#               #)
#             )$result %>% list())
#    # stretch the branch lengths
#    phylo_indl$edge.length = phylo_indl$edge.length * mult_gamma
#    # unpack the fit
#    #fit %>% broom::tidy() %>% unnest_wider() %>% print()
#    # and have a conditional here to break the loop
#  }
#  print(fits)
#  #return(fit)
#}
