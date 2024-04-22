## Ingest c0 values from literature, generate estimates
## to apply to lipidomic data

library(tidyverse)
library(here)
library(httr)
library(ggpubr)

file_curv = here("01-rawdata", "curvature.csv")
# this is a Google Sheets key for the curvature info file, to facilitate a nice Excel-like interface
gsht_curv = "1ec7xtebfnTYsC9I6BAfYiGs-C-wzmpZGDI1OtyDBHcA"

# refresh c0 spreadsheet from online interface
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_curv}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_curv)

# and load the data in
data_curv = file_curv %>% read_csv()

# let's just see all of them!
plot_allc0 = data_curv %>% 
  factorize_lipids() %>% 
  pivot_longer(
    cols = c("carbsn1", "dbonsn1"),
    names_to = "var",
    values_to = "val"
  ) %>% 
  arrange(var, class) %>% 
  mutate(varclass = paste(var, class) %>% factor(levels=unique(.))) %>%
  ggplot(
    aes(
      x = val,
      y = c0
    )
  ) +
  facet_wrap(~varclass, scale = "free_x", nrow = 2) +
  geom_point(
    aes(
      color = temp,
      shape = !is.na(relax)
    )
  ) +
  geom_smooth(method="lm", se=FALSE) +
  scale_colour_viridis_c() +
  scale_shape_manual(values = c(1, 16))

#plot_allc0

# filter to 0 bar values
curv_rated = data_curv %>% 
  filter(press == 0) %>% 
  # generate a quality score of sorts; higher is better
  mutate(
    temp_std = 20,
    goodness = 
      as.integer(!est) + # a point for not being estimated
      as.integer(!is.na(relax)) + # a point for having filler
      as.integer(str_detect(buff, "H2O")) + # a point for being in water
      as.integer(npln) + # one for being at neutral plane
      as.integer((temp == temp_std) | !is.na(tslope)) # a point for "standard temp" or at least a correction factor
  )

## SCHEME 1: 

# prioritize sn-1 18:1, try to adjust to 20°C
c0_oleoyl_best = curv_rated %>% 
  # simple acyl chain criteria
  filter(
    # sn-1 needs to be oleoyl
    (carbsn1 == 18) &
      (dbonsn1 == 1 ) |#&
      ## needs to be lyso-, or
      #(str_detect(class, 'L') |
      ## sn-2 also needs to be oleoyl
      #(carbsn2 == 18) &
      #(dbonsn2 == 1 )) |
      # dioleoyl plasmalogen not available!
      str_detect(class, "PPE")
  ) %>% 
  arrange(class, -goodness) %>% 
  # grab just the top-quality data for each class
  group_by(class) %>% 
  filter(goodness == max(goodness)) %>% 
  ungroup() %>% 
  # and adjust for temperature
  mutate(
    temp_std = 20,
    c0  = c0 + ifelse(!is.na(tslope), tslope * (temp_std - temp), 0),
    tol = tol + abs(ifelse(!is.na(slopetol), ifelse(!is.na(tslope), tslope, 0) * (temp_std - temp), 0)),
    temp = ifelse(is.na(tslope), temp, temp_std)
  )
# then average
c0_oleoyl_mean = c0_oleoyl_best %>% 
  # average by class
  group_by(class) %>% 
  summarize(
    c0 = mean(c0),
    tol = max(tol, na.rm = TRUE),
    tol = ifelse(tol == -Inf, NA, tol),
    temp = mean(temp),
    across(contains("sn1"), mean),
    goodness = mean(goodness),
    ref = list(ref)
  )
# might want to add an adjustment here so PPE looks "saturated"?
# or might not, for conservatism

## SCHEME 2

# sn-1 saturated, sn-2 oleoyl, try to adjust to 20°C.
c0_satsn1_best = curv_rated %>% 
  # simple acyl chain criteria
  filter(
    (dbonsn1 == 0) &
      ((carbsn2 == 18) & (dbonsn2 == 1) |
       str_detect(class, 'L')) |
      (class == "PS") # (no saturated PSs)
  ) %>% 
  arrange(class, -goodness) %>% 
  # grab just the top-quality data for each class
  group_by(class) %>% 
  filter(goodness == max(goodness)) %>% 
  ungroup() %>% 
  # and adjust for temperature
  mutate(
    temp_std = 20,
    c0  = c0 + ifelse(!is.na(tslope), tslope * (temp_std - temp), 0),
    tol = tol + abs(ifelse(!is.na(slopetol), ifelse(!is.na(tslope), tslope, 0) * (temp_std - temp), 0)),
    temp = ifelse(is.na(tslope), temp, temp_std)
  )
# then average
c0_satsn1_mean = c0_satsn1_best %>% 
  # average by class
  group_by(class) %>% 
  summarize(
    c0 = mean(c0),
    tol = max(tol, na.rm = TRUE),
    tol = ifelse(tol == -Inf, NA, tol),
    temp = mean(temp),
    across(contains("sn1"), mean),
    goodness = mean(goodness),
    ref = list(ref)
  )

## SCHEME 3

# build carbsn1 x dbonsn1 regressions for PE, PC only
# It's evident in the plot below that chain length and unsat effects
# on the lysolipids are not statistically significant,
# and just the means should be used.
# use only really good data
linreg_data = curv_rated %>% 
  filter(
    (class %in% c("PC", "PE", "PPE", "LPC", "LPE")) &
      !is.na(relax) &
      str_detect(buff, "H2O") &
      (npln | str_detect(class, 'L')) & # none of the lysos are neutral-plane
      ((temp %in% c(20, 22)) | !is.na(tslope))
  ) %>% 
  # adjust for temperature
  mutate(
    temp_std = 20,
    c0  = c0 + ifelse(!is.na(tslope), tslope * (temp_std - temp), 0),
    temp = ifelse(is.na(tslope), temp, temp_std)
  )

# let's see those groomed data!
plot_linreg = linreg_data %>% 
  pivot_longer(
    cols = c("carbsn1", "dbonsn1"),
    names_to = "var",
    values_to = "val"
  ) %>% 
  mutate(varclass = paste(var, class) %>% factor(levels=unique(.))) %>%
  ggplot(
    aes(
      x = val,
      y = c0,
      color = class,
      group = class
    )
  ) +
  facet_wrap(~var, scale = "free_x", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE) +
  scale_color_manual(values = chroma_cl) +
  labs(
    title = "c0 by sn-1 chain length and unsat\nadjusted as necessary to 20°C"
  )
#plot_linreg

# a very simple model
# it would be nice if errors could be incorporated...oh well
linreg = linreg_data %>% 
  filter(class %in% c("PC", "PE")) %>% # leaving lysolipids out for reason given above
  # Do the carbon and dbond effects differ between classes? 
  #lm(c0 ~ (carbsn1 + dbonsn1) * class, .) # nope!
  #lm(c0 ~ (carbsn1 + dbonsn1) %in% class, .) # within-class dbond effects not significant
  # pooling chain effects across PC and PE
  lm(c0 ~ ((dbonsn1 + carbsn1) %in% carbsn1) + class, .) # this seems a pretty reasonable model! 

summary(linreg)

# what's the c0 offset associated with PPE vs PE?
offs_ppe = curv_rated %>% 
  filter(
    (class %in% c("PE", "PPE")) & 
      (carbsn1 == 18) & (dbonsn1 == 0) &
      (relax == "TS")
    ) %>% 
  arrange(-c0) %>% .$c0 %>% {(min(.) - max(.))}

# what does that model predict?
c0_chain_fx = crossing(
  class   = c("PC", "PE"),
  carbsn1 = seq(10, 22, 1),
  dbonsn1 = c(0, 1),
  temp    = 20
) %>% 
  rowwise() %>% 
  mutate(
    # run predictions w/95% CI
    c0  = predict(linreg, cur_data(), interval = "confidence", level = 0.90) %>% 
      as.list() %>% setNames(c("c0", "lo", "hi")) %>% list(),
    # insert the lm formula as ref
    ref = linreg$call %>% as.character() %>% .[[2]] %>% list()
  ) %>% 
  unnest_wider(c0) %>% 
  mutate(tol = c0-lo) %>% 
  select(-lo, -hi) %>% 
  # add PPE values with a c0 offset given by SOPE - SOPPE
  bind_rows(tibble(class="PPE", carbsn1=18, dbonsn1=0, temp=20)) %>% 
  complete(class, carbsn1, dbonsn1, temp) %>% 
  group_by(carbsn1, dbonsn1) %>% 
  arrange(class) %>% 
  mutate(
    ref = ifelse(!is.na(c0), ref, list("~SOPE-SOPPE offset")),
    tol = ifelse(!is.na(c0), tol, tol[[2]]),
    c0  = ifelse(!is.na(c0), c0, min(c0, na.rm = TRUE)+offs_ppe)
  )

# plot those results
plot_curvmodel = c0_chain_fx %>% 
  ggplot(
    aes(
      x = carbsn1,
      y = c0,
      ymin = c0-tol,
      ymax = c0+tol,
      #fill = class,
      color = class
    )
  ) +
  # perhaps should be point?
  #geom_col(position = position_dodge()) +
  geom_line() +
  geom_errorbar(position = position_dodge()) +
  facet_wrap(~dbonsn1) +
  #scale_fill_manual(values = chroma_cl) +
  scale_color_manual(values = chroma_cl) +
  theme_pubr() +
  labs(
    title = "Predicted curvature as function of headgroup, sn-1 length and unsat."
  )
#plot_curvmodel

## apply a "salt correction" on top of chain_fx
## addressing reviewer 3 comment
#saltfac_anionic = data_curv %>% 
#  filter(class == "PS") %>% 
#  filter(str_detect(buff, "10 mM MgCl2")) %>% 
#  # -0.041 is the value Dymond chooses to cite from his own paper, with 3% guest
#  # interesting, though: apparent c0 only varies with mole frac at low salt concs
#  # choosing the largest absolute value here is definitely most conservative
#  # keep in mind, though, that the prep was not relaxed!
#  arrange(c0) %>% head(1) %>% 
#  bind_rows(c0_oleoyl_mean %>% filter(class == "PS") %>% unnest(ref)) %>% 
#  # calculate the fraction of curvature to subtract
#  arrange(c0) %>% 
#  .$c0 %>% {(max(.)-min(.))/max(.)} %>% abs()
#  # i.e. for anionic lipids, subtract 6.94 * c0!

## better to compare the salted PS value to its unsalted counterpart from same paper
#saltfac_anionic = data_curv %>% 
#  filter(class == "PS") %>% 
#  filter(ref == "https://doi.org/10.1021/acs.langmuir.6b03098") %>% 
#  #filter(str_detect(buff, "10 mM MgCl2")) %>% 
#  View
#  
#  # -0.041 is the value Dymond chooses to cite from his own paper, with 3% guest
#  # interesting, though: apparent c0 only varies with mole frac at low salt concs
#  # choosing the largest absolute value here is definitely most conservative
#  # keep in mind, though, that the prep was not relaxed!
#  arrange(c0) %>% head(1) %>% 
#  bind_rows(c0_oleoyl_mean %>% filter(class == "PS") %>% unnest(ref)) %>% 
#  # calculate the fraction of curvature to subtract
#  arrange(c0) %>% 
#  .$c0 %>% {(max(.)-min(.))/max(.)} %>% abs()
## i.e. for anionic lipids, subtract 6.94 * c0!
#  
#saltfac_zwitter = data_curv %>% 
#  filter(class == "PE") %>% 
#  filter(str_detect(buff, "150 mM NaCl")) %>% 
#  # take the one without Ca2+; it's too much calcium! (25 mM)
#  arrange(-c0) %>% head(1) %>% 
#  bind_rows(c0_oleoyl_mean %>% filter(class == "PE") %>% unnest(ref)) %>% 
#  # calculate the fraction of curvature to subtract
#  arrange(c0) %>% 
#  .$c0 %>% {(max(.)-min(.))/max(.)} %>% abs()
## i.e. for zwitterionic lipids, subtract 0.15 * c0!

# bind all the schemes together into one table
c0_allschemes = bind_rows(
  c0_oleoyl_mean %>% mutate(scheme = "oleoyl"),
  c0_satsn1_mean %>% mutate(scheme = "satsn1"),
  c0_chain_fx    %>% mutate(scheme = "linreg"),
) %>% 
  ## these are all curvatures at 0 bar! redundant here?
  #filter(press == 0) %>% 
  # fill in values for the linreg scheme with an average from other two schemes
  complete(class, scheme) %>% 
  group_by(class) %>% 
  mutate(
    ref = ifelse(!is.na(c0), ref, cur_data() %>% .$ref %>% unlist() %>% unique() %>% unname()), 
    tol = ifelse(!is.na(c0), tol, cur_data() %>% .$tol %>% mean(na.rm=TRUE)), # will return NaN if all NA
    c0  = ifelse(!is.na(c0), c0,  cur_data() %>% .$c0  %>% mean(na.rm=TRUE)),
    #NTS 20230710 get model-based errors in here too!
  ) %>% 
  # finally, duplicate the linreg model and apply salt corrections
  (function(x){
    bind_rows(
      x,
      x %>% 
        filter(scheme == "linreg") %>% 
        mutate(scheme = "salted")
    )
  }) %>% 
  # join the correction factors by lipid type
  left_join(
    bind_rows(
      crossing(
        class = c("PS", "PG"),
        saltfac = saltfac_anionic
      ),
      crossing(
        class = c("PPE", "PE", "PC", "LPE", "LPC"),
        saltfac = saltfac_zwitter
      )
    ) %>% 
      mutate(scheme = "salted"),
    by = c("class", "scheme")
  ) %>% 
    # then apply the factors
    mutate(c0 = ifelse(scheme == "salted", c0 - (saltfac*abs(c0)), c0))
