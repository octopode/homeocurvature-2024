## Estimate mean curvature from PL composition data
## Runs a few different estimation schemes,
## one of which (linreg) takes sn-1 chain structure into account

library(tidyverse)
library(here)
library(ggpubr)

source(here("03-scripts", "prep_lipidmaps.R"))
source(here("03-scripts", "c0_meta_analysis.R"))

# remove sn1 chain data for anionic PLs bc they are not predictors in the linreg model
# added 20231228 bc PS and PG were getting dropped!
pldata_anionna = pldata %>% mutate(
  carbsn1 = ifelse(class %in% c("PG", "PS"), NA, carbsn1),
  dbonsn1 = ifelse(class %in% c("PG", "PS"), NA, dbonsn1)
)

# join curvature data to lipidomic data
plcidata = bind_rows(
  pldata %>% mutate(scheme = "oleoyl"),
  pldata %>% mutate(scheme = "satsn1"),
  pldata_anionna %>% mutate(scheme = "linreg")
) %>% 
  # join headgroup-only schemes
  left_join(
    c0_allschemes %>% 
      filter(scheme %in% c("oleoyl", "satsn1")) %>% 
      select(scheme, class, c0, tol),
    by = c("scheme", "class")
  ) %>% 
  # join sn1-aware scheme
  mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
  left_join(
    c0_allschemes %>% 
      filter(scheme == "linreg") %>% 
      select(scheme, class, carbsn1, dbonsn1, c0, tol) %>% 
      mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
      select(-dbonsn1),
    by = c("scheme", "class", "carbsn1", "unsatsn1"),
    suffix = c('', ".linreg")
  ) %>% 
  mutate(
    c0  = ifelse(is.na(c0 ), c0.linreg, c0 ),
    tol = ifelse(is.na(tol), tol.linreg,  tol)
  ) %>% 
  select(-contains(".linreg")) %>% 
  # finally, calculate the actual curvature contributions and error
  # `plci` stands for PhosphoLipid Curvature Index
  mutate(
    plci = c0  * frac_molar,
    ctol = tol * frac_molar # should this scale linearly or by sqrt?
  ) %>% 
  # important!
  replace_na(list(plci = 0, ctol = 0))

# just whole and body samples
plcidata_wildwhole = plcidata %>% 
  filter(wild & (tissue %in% c("whole", "body")))

# for E. coli with QC'd data
plcidata_ecoli = bind_rows(
  pldata_ecoli %>% mutate(scheme = "oleoyl"),
  pldata_ecoli %>% mutate(scheme = "satsn1"),
  pldata_ecoli %>% mutate(scheme = "linreg")
) %>% 
  # join headgroup-only schemes
  left_join(
    c0_allschemes %>% 
      filter(scheme %in% c("oleoyl", "satsn1")) %>% 
      select(scheme, class, c0, tol),
    by = c("scheme", "class")
  ) %>% 
  # join sn1-aware scheme
  mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
  left_join(
    c0_allschemes %>% 
      filter(scheme == "linreg") %>% 
      select(scheme, class, carbsn1, dbonsn1, c0, tol) %>% 
      mutate(unsatsn1 = as.logical(dbonsn1)) %>% 
      select(-dbonsn1),
    by = c("scheme", "class", "carbsn1", "unsatsn1"),
    suffix = c('', ".linreg")
  ) %>% 
  mutate(
    c0  = ifelse(is.na(c0 ), c0.linreg, c0 ),
    tol = ifelse(is.na(tol), tol.linreg,  tol)
  ) %>% 
  select(-contains(".linreg")) %>% 
  # finally, calculate the actual curvature contributions and error
  # `plci` stands for PhosphoLipid Curvature Index
  mutate(
    plci = c0  * frac_molar,
    ctol = tol * frac_molar # should this scale linearly or by sqrt?
  ) %>% 
  # important!
  replace_na(list(plci = 0, ctol = 0))
