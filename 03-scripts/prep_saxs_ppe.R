# Analysis and plotting for SAXS analyses of low-complexity synthetic lipid systems

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pbapply)

#source(here("02-scripts", "20220502_hpsaxs_helpers.R"))
source(here("03-scripts", "saxs_helpers.R"))

# locate data files
dir_head = here("01-rawdata", "saxs_ppe")
pat_dat = "\\.dat$"
files_data = list.files(dir_head, pat_dat, recursive = TRUE, full.names = TRUE)
files_data = files_data[which(str_detect(files_data, "01-raw"))]

# identify the samples
# can identify series by intersections of samp and temp
composition = c(
  "S1"      = "18:0/20:4 PE",
  "S2"      = "18:0/20:4 PE-PPE 50-50",
  "S3"      = "18:0/20:4 PPE",
  "JWL219A" = "18:0/20:4 PE",
  "JWL219B" = "18:0/20:4 PE",
  "JWL220A" = "18:0/20:4 PE-PPE 50-50",
  "JWL222A" = "18:0/18:1 PE",
  "JWL223A" = "18:0/18:1 PPE"
)

# load all the data
data_saxs_ppe = files_data %>% 
  #NTS 20220503: the timestamps don't parse
  # but they do for the May 2021 samples using the same load function
  # chalk it up to a change in BIOXTAS; fine for now
  read_saxsall() %>% 
  # group by run
  group_by(Time, date, pdir, press, temp, rep) %>% 
  mutate(
    press = 10*press, #convert to bar
    rep  = as.numeric(rep),
    comp = composition[samp]
  ) %>% 
  group_by(samp, temp, comp) %>% 
  # create a (P-sweep) series variable
  mutate(ser = cur_group_id())

# normalize each shot to area
data_saxs_ppe_norm = data_saxs_ppe %>% 
  # trim
  filter(q >= 0.00 & q <= 0.5) %>% 
  group_by(fname) %>% 
  mutate(
    iq = iq - min(iq),
    iq = iq / sum(iq)
  )

# inspect a series to manually classify frames
data_saxs_ppe_norm %>%
  filter(pdir == "up") %>% 
  #filter(comp == "18:0/18:1 PPE") %>% 
  filter(comp == "18:0/20:4 PPE") %>% 
  filter(ser == 10) %>% 
  filter(press %in% c(0, 1000, 1500)) %>% 
  mutate(ser = paste(ser, comp, temp)) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      color = pdir
    )
  ) +
  facet_grid(
    rows = vars(press), 
    cols = vars(ser),
    scales = "free_y"
  ) +
  geom_line() +
  #geom_vline(xintercept = 0.1) +
  theme_pubr() +
  lims(x = c(0,0.3))

# OK, now my eyeballs have done some of the work
# load in annotations
phranges = read_tsv(here("01-rawdata", "saxs_ppe", "phaseranges_20220503.tsv"))

# extract pure phase exemplars from these data
phprofs  = data_saxs_ppe_norm %>% 
  left_join(
    phranges %>% 
      # "pure" phases only
      filter(hii + la + lb == 1), 
    by = c("ser", "temp", "comp")
  ) %>% 
  drop_na(press_min, press_max) %>% 
  pivot_longer(c(hii, la, lb), names_to = "ph", values_to = "ya") %>% 
  filter(ya) %>% select(-ya) %>% 
  filter((press >= press_min) & (press <= press_max)) %>% 
  group_by(ser, temp, comp, ph, q) %>% 
  summarize(iq = mean(iq)) %>% 
  group_by(ser, temp, comp, ph) %>%
  # renormalize
  mutate(iq = iq/sum(iq))

# look at them
phprofs %>%
  ggplot(
    aes(
      x = q,
      y = iq
    )
  ) +
  facet_wrap(~paste(ser, comp, temp, ph)) +
  geom_line() +
  theme_pubr()
# That'll do for now. Would be nice to figure out logistic extrapolation!

data_saxs_ppe_frac = data_saxs_ppe_norm %>% 
  # subset data
  filter(
    (temp == 20 & str_detect(comp, "18:0/20:4")) |
      (temp == 60 & str_detect(comp, "18:0/18:1"))
  ) %>% 
  left_join(phranges, by = c("temp", "comp", "ser")) %>% 
  filter((press >= press_min) & (press <= press_max)) %>% #filter(comp == "18:0/18:1 PPE") %>% filter(press == 900) %>% .$lb %>% unique()
  full_join(
    phprofs %>% 
      pivot_wider(names_from = "ph", names_prefix = "iq_", values_from = "iq") %>% 
      replace_na(list("iq_hii" = 0, "iq_la" = 0, "iq_lb" = 0)), 
    by = c("q", "temp", "comp", "ser")
  ) %>% 
  drop_na(press) %>% # introduced by full_join
  group_by(ser, temp, comp, press, pdir, fname) %>% # include pdir or no?
  mutate(
    # assert presence/absence of terms
    iq_hii = hii*iq_hii,
    iq_la = la*iq_la,
    iq_lb = lb*iq_lb,
  ) %>% #filter(comp == "18:0/18:1 PPE") %>% filter(lb) %>% arrange(press)
  # fit model
  summarize(mod = lm(iq ~ iq_hii + iq_la + iq_lb + 0, data = cur_data()) %>% list()) %>% 
  # extract coefficients
  mutate(coeffs = mod[[1]]$coefficients %>% list()) %>% 
  unnest_wider(coeffs) %>% 
  # make them add to 1
  pivot_longer(cols = contains("iq"), names_to = "ph", values_to = "frac") %>%
  replace_na(list("frac" = 0)) %>% 
  # I don't know why fname grouping is lost
  group_by(ser, temp, comp, press, fname) %>% 
  mutate(
    ph = str_remove(ph, "iq_"),
    frac = frac-min(frac),
    frac = frac/sum(frac)
  ) %>% 
  # prep for downstream
  # these are the good data
  filter(ser %in% c(1, 4, 10, 5, 6)) %>% 
  separate(comp, into = c("chains", NA), sep= ' ', remove = FALSE, extra = "drop") %>% 
  # strip duplicates (not sure why they're present?) 30 extra rows for no good reason...
  group_by(ser, comp, chains, press, pdir, ph, frac) %>% 
  summarize() %>% 
  # need to split la for PPE into two subsets for separate logregs. foosh!
  mutate(logreg = c(1, 2) %>% list()) %>% 
  unnest(logreg) %>% 
  # push data into their respective regressions
  filter(
    (chains == "18:0/20:4") & 
      ((comp != "18:0/20:4 PPE") & (logreg == 1)) |
      ((press <= 1000) & (logreg == 2) & (comp == "18:0/20:4 PPE")) |
      ((press >= 1000) & (logreg == 1)) |
      (chains == "18:0/18:1") & 
      ((comp != "18:0/18:1 PPE") & (logreg == 1)) |
      ((press <= 1000) & (logreg == 2) & (comp == "18:0/18:1 PPE")) |
      ((press >= 1000) & (logreg == 1))
  )

#NTS 20220511 some serious refactoring of the below

# do logistic smoothing
data_saxs_ppe_mods = data_saxs_ppe_frac %>% 
  # downweight unpaired shots
  group_by(comp, press, ph, logreg) %>% 
  mutate(wts = n()) %>% # weighted, nice
  #mutate(wts = 1) %>% # unweighted
  ## try to fix slope artifact
  #mutate(wts = ifelse((comp == "18:0/20:4 PPE") & (press == 1500), -1, wts)) %>%
  # now for the regressions!
  group_by(ser, comp, ph, logreg) %>% 
  summarise(
    mod_lgt = glm(
      data = cur_data(),
      formula = frac ~ press,
      family = "binomial",
      weights = wts
    ) %>% list()
  )

# now we are gonna try to fix that slope artifact by fixing the slope
# to the mean of the others and fitting the intercept
slope_fix = data_saxs_ppe_mods %>% 
  filter(!((logreg == 1) & (ph == "lb") & (comp == "18:0/20:4 PPE"))) %>% 
  #filter(((logreg == 1) & (ph == "lb") & (comp == "18:0/20:4 PE-PPE 50-50"))) %>%
  rowwise() %>% 
  mutate(
    slope = mod_lgt %>% summary() %>% .$coefficients %>% .[[2]],
    slope = ifelse((logreg == 1) & (ph == "lb"), slope, NA)
  ) %>% .$slope %>% mean(., na.rm = TRUE)

constr_mod = data_saxs_ppe_frac %>% 
  # downweight unpaired shots
  group_by(comp, press, ph, logreg) %>% 
  mutate(wts = n()) %>% # weighted, nice
  #mutate(wts = 1) %>% # unweighted
  ## try to fix slope artifact
  #mutate(wts = ifelse((comp == "18:0/20:4 PPE") & (press == 1500), -1, wts)) %>%
  # now for the regressions!
  group_by(ser, comp, ph, logreg) %>% 
  filter((logreg == 1) & (ph == "lb") & (comp == "18:0/20:4 PPE")) %>% 
  summarise(
    mod_lgt = glm(
      data = cur_data(),
      formula = frac ~ offset(slope_fix*press),
      family = "binomial",
      weights = wts
    ) %>% list()
  )

# slip in the constrained model
data_saxs_ppe_mods = bind_rows(
  data_saxs_ppe_mods %>% filter(!((logreg == 1) & (ph == "lb") & (comp == "18:0/20:4 PPE"))),
  constr_mod
)

# run out the fits
domain_fit = seq(0, 2000, 10)
data_saxs_ppe_fits = data_saxs_ppe_mods %>% 
  group_by(ser, comp, ph, logreg) %>% 
  filter(
    !((logreg == 2) & (ph == "lb")) &
      !((logreg == 1) & (ph == "hii"))
  ) %>% 
  mutate(press = list(domain_fit)) %>% 
  unnest(press) %>% 
  mutate(frac = predict(mod_lgt[[1]], newdata = cur_data(), type = "response")) %>% 
  select(-contains("mod")) %>% 
  # "meld" back-to-back La regressions
  group_by(ser, comp, ph, press) %>% 
  filter((frac == min(frac)) | (ph != "la")) %>% 
  # normalize them so the plot winds up rectangular
  group_by(ser, comp, press) %>% 
  mutate(frac = frac/sum(frac))

##TEST sanity check!
#data_saxs_ppe_fits %>% 
#  ggplot(
#    aes(x = press, y = frac, color = ph, group = paste(ph, logreg))
#  ) +
#  facet_wrap(~comp) +
#  geom_line(aes(size = -logreg), alpha = 0.5)

# master plot
data_saxs_ppe_fits %>% 
  ggplot(
    aes(
      x = press, 
      y = frac, 
      fill = ph
    )
  ) +
  facet_grid(cols = vars(comp)) +
  geom_area() + # the models, stacked by default
  geom_point(
    data = data_saxs_ppe_frac %>% 
      # avoid redundant plotting
      filter(
        ((logreg == 1) & (ph == "lb")) |
          ((logreg == 2) & (ph == "la"))
      ) %>% 
      # recode "bottom" points
      mutate(pdir = ifelse(press == 2000, "bt", pdir)),
    aes(shape = pdir),
    #position = position_dodge(100),
    fill = "white",
    color="black"
  ) +
  theme_pubr() +
  # orientation control
  coord_flip() +
  scale_x_reverse() +
  # show bidirectional data
  scale_fill_manual(values = chroma_ph) +
  scale_shape_manual(values = c("dn" = 24, "up" = 25, "bt" = 22)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(
    fill  = "none",
    shape = "none"
  )


# plot them up the way I want
data_saxs_ppe_frac %>% 
  # just PUFA
  filter(ser %in% c(1,4,10)) %>% 
  # need to split la for PPE into two subsets for separate logregs. foosh!
  mutate(logreg = c(1, 2) %>% list()) %>% 
  unnest(logreg) %>% 
  filter(
    ((comp != "18:0/20:4 PPE") & (logreg == 1)) |
      ((press <= 1000) & (logreg == 2) & (comp == "18:0/20:4 PPE")) |
      ((press >= 1000) & (logreg == 1))
  ) %>% 
  # downweight unpaired shots
  group_by(comp, press) %>% 
  mutate(wts = n()) %>% # weighted, nice
  #mutate(wts = 1) %>% # unweighted
  # workaround overlap by ensuring points don't stack
  mutate(press = ifelse((logreg == 2) & (press <= 1000), press+1, press)) %>% 
  ggplot(
    aes(
      x = press,
      y = frac,
      fill = ph,
      group = paste(ph, logreg)
    )
  ) +
  facet_grid(rows = vars(comp)) +
  # logistic fits: prettier!
  stat_smooth(
    geom = "area",
    position = position_stack(),
    method = "glm", 
    se = FALSE, 
    mapping = aes(weight = wts),
    method.args = list(family = "binomial")
  ) +
  scale_fill_manual(values = chroma_ph) +
  geom_point(
    data = data_saxs_ppe_frac %>% 
      # just PUFA
      filter(ser %in% c(1,4,10)) %>% 
      mutate(logreg = c(1, 2) %>% list()) %>% 
      unnest(logreg) %>% 
      filter(
        ((comp != "18:0/20:4 PPE") & (logreg == 1)) |
          ((press <= 1000) & (logreg == 2) & (comp == "18:0/20:4 PPE")) |
          ((press >= 1000) & (logreg == 1))
      ) %>% 
      # avoid redundant plotting
      group_by(comp, press, logreg) %>% 
      arrange(ph) %>% 
      filter(
        ((ph[[1]] == "hii") & (frac[[1]] > 0) & (ph == "la")) |
          ((frac[[1]] == 0) & (ph == "lb"))
      ),
    aes(shape = pdir),
    #position = position_dodge(100),
    color="black"
  )  +
  theme_pubr() +
  # show bidirectional data
  scale_shape_manual(values = c("up" = 24, "dn" = 25)) +
  guides(
    fill  = "none",
    shape = "none"
  )
ggsave(here("03-mainfigs", "synsim", "20220503_panelC.pdf"), width=8, height=4)

# same thing but for MUFA.
data_saxs_ppe_frac %>% 
  # just MUFA
  filter(ser %in% c(5, 6)) %>% 
  # strip duplicates (not sure why they're present?) 30 extra rows for no good reason...
  group_by(ser, comp, press, ph, frac) %>% 
  summarize() %>% 
  # some QC
  # assert no lb below 1000 bar and renorm
  mutate(frac = ifelse(str_detect(comp, "PPE") & (ph == "lb") & (press < 1000), 0, frac)) %>% 
  group_by(ser, comp, press) %>% 
  mutate(frac = frac/sum(frac)) %>% 
  # need to split la for PPE into two subsets for separate logregs. foosh!
  mutate(logreg = c(1, 2) %>% list()) %>% 
  unnest(logreg) %>% 
  filter(
    ((comp != "18:0/18:1 PPE") & (logreg == 1)) |
      ((press <= 1000) & (logreg == 2) & (comp == "18:0/18:1 PPE")) |
      ((press >= 1000) & (logreg == 1))
  ) %>% 
  # workaround overlap by ensuring points don't stack
  mutate(press = ifelse((logreg == 2) & (press <= 1000), press+1, press)) %>% 
  arrange(press) %>% 
  ggplot(
    aes(
      x = press,
      y = frac,
      fill = ph,
      group = paste(ph, logreg)
    )
  ) +
  facet_grid(rows = vars(comp)) +
  # logistic fits: prettier!
  stat_smooth(
    geom = "area",
    position = position_stack(),
    method = "glm", 
    se = FALSE, 
    method.args = list(family = "binomial")
  ) +
  scale_fill_manual(values = chroma_ph) +
  #geom_area(position = position_stack()) + #TEST, raw area plot
  geom_point(color="black") +
  theme_pubr() +
  guides(fill = "none")
ggsave(here("03-mainfigs", "synsim", "20220503_MUFA.pdf"), width=8, height=4)

## Phase exemplar profiles for Fig. S7/8 (spectroscopy of inversion)
data_saxs_ppe_norm %>% 
  filter(q <= 0.3) %>% 
  filter(str_detect(comp, "20:4")) %>% 
  filter(!str_detect(comp, '-')) %>% 
  unite(ser, c(ser, comp)) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      color = temp
    )
  ) +
  facet_grid(col = vars(ser), rows = vars(press), scales = "free_y") +
  geom_line() +
  theme_pubr()
ggsave(here("04-suppfigs", "S7-LaurdanInversion", "20220724_synsaxs_explore.pdf"), width=8, height=14)

# OK, I picked some out
data_saxs_ppe_norm %>% 
  filter(q <= 0.3) %>% 
  filter(ser %in% c(1,2,9,10)) %>% 
  filter(
    (
      (temp == 4) &
        (press %in% c(0, 500))
    ) |
      (
        (temp == 20) &
          (press %in% c(0, 1000))
      )
  ) %>% 
  group_by(comp, temp, press) %>% 
  arrange(desc(pdir)) %>% 
  filter(fname == first(fname)) %>% 
  group_by(fname) %>% 
  mutate(iq = iq/max(iq)) %>% 
  mutate(iq = 100* iq + 1) %>% 
  unite(ser, c(comp, temp), remove = FALSE) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      color = temp,
      group = fname
    )
  ) +
  facet_grid(rows = vars(press), cols = vars(ser)) +
  geom_line() +
  scale_y_log10() +
  theme_pubr()
ggsave(here("04-suppfigs", "S7-LaurdanInversion", "20220724_synsaxs_PUFA_exemplars.pdf"), width=8, height=4)

# for minimal reformatting
# OK, I picked some out
data_saxs_ppe_norm %>% 
  filter(q <= 0.3) %>% 
  filter(ser %in% c(1,2,9,10)) %>% 
  filter(
    (
      (temp == 4) &
        (press %in% c(0, 500))
    ) |
      (
        (temp == 20) &
          (press %in% c(0, 1000))
      )
  ) %>% 
  group_by(comp, temp, press) %>% 
  arrange(desc(pdir)) %>% 
  filter(fname == first(fname)) %>% 
  group_by(fname) %>% 
  mutate(iq = iq/max(iq)) %>% 
  mutate(iq = 100* iq + 1) %>% 
  mutate(press = press>0) %>% 
  mutate(temp = c("4"="cold", "20"="warm")[as.character(temp)]) %>% 
  unite(ser, c(comp, temp), remove = FALSE) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      color = temp,
      group = fname
    )
  ) +
  facet_grid(rows = vars(press), cols = vars(ser)) +
  geom_line() +
  scale_y_log10() +
  theme_pubr()
ggsave(here("04-suppfigs", "S7-LaurdanInversion", "20220724_synsaxs_PUFA_exemplars_aligned_1000.pdf"), width=8, height=4)


