library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pbapply)

source(here("03-scripts", "saxs_helpers.R"))

# locate data files
dir_head = here("01-rawdata", "saxs_lam")
pat_dat = "\\.dat$"
files_saxsdatalam = list.files(dir_head, pat_dat, recursive = FALSE, full.names = TRUE)
files_saxsdatalam = files_saxsdatalam[which(str_detect(files_saxsdatalam, "01-raw"))]

# identify the samples
# can identify series by intersections of samp and temp
composition_lam = c(
  "JWL224A" = "18:1/18:1 PE",
  "JWL225A" = "18:1/18:1 PC"
)

# load all the data
data_saxs_lam = files_saxsdatalam %>% 
  #NTS 20220503: the timestamps don't parse
  # but they do for the May 2021 samples using the same load function
  # chalk it up to a change in BIOXTAS; fine for now
  read_saxsall() %>% 
  # group by run
  group_by(Time, date, pdir, press, temp, rep) %>% 
  mutate(
    press = 10*press, #convert to bar
    rep  = as.numeric(rep),
    comp = composition_lam[samp]
  ) %>% 
  group_by(samp, temp, comp) %>% 
  # create a (P-sweep) series variable
  mutate(ser = cur_group_id())

## normalize each shot to area
#data_saxs_lam_norm = data_saxs_lam %>% 
#  # trim
#  filter(q >= 0.00 & q <= 0.5) %>% 
#  group_by(fname) %>% 
#  mutate(
#    iq = iq - min(iq),
#    iq = iq / sum(iq)
#  )

# plot to verify first peak is highest
data_saxs_lam %>% 
  mutate(row = paste(comp, pdir)) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      color = pdir
    )
  ) +
  facet_grid(rows = vars(row), cols = vars(press), scales = "free_y") +
  geom_line()

# calculate D-spacings
# just use the q of the tallest peak!
# which is > 0.25 in DOPC, > 5.0 in DOPE
peaqs_lam = data_saxs_lam %>% 
  separate(comp, into = c("chains", "class"), sep = ' ', remove = FALSE) %>% 
  group_by(fname, comp, pdir, press) %>% 
  filter(
    ((comp == "18:1/18:1 PC") & (iq > 0.25)) |
      ((comp == "18:1/18:1 PE") & (iq > 5.0))
  ) %>% 
  arrange(q) %>% 
  filter((iq > lead(iq)) & (iq > lag(iq))) %>% 
  mutate(d = 2*pi/q) %>% 
  # there is a phase break!
  mutate(ph = ifelse(press > 1800, "lb", "la")) %>% 
  # calc delta d
  group_by(comp, ph, pdir) %>% 
  arrange(press) %>% 
  mutate(dd = d-first(d)) %>% 
  mutate(comph = paste(comp, ph))

# plot d vs. P
peaqs_lam %>% 
  ggplot(
    aes(
      x = press,
      y = dd,
      group = paste(comp, ph)
    )
  ) +
  facet_wrap(~comp) +
  geom_point(aes(shape = pdir)) +
  geom_smooth(method = "lm") +
  scale_shape_manual(values = c(up = 24, dn = 25)) +
  theme_pubr()

# plot delta d vs. P
peaqs_lam %>% 
  ggplot(
    aes(
      x = press,
      y = dd,
      group = comph
    )
  ) +
  facet_wrap(~comph, scales = "free_x") +
  geom_smooth(color = "black", method = "lm") +
  geom_point(
    size = 2,
    color = "white",
    aes(
      shape = pdir,
      fill = class
    )
  ) +
  # calc linreg slopes and print them
  geom_text(
    data = peaqs_lam %>%
      group_by(comp, ph) %>% 
      mutate(xloc = mean(range(press))) %>% 
      group_by(comp, ph, xloc, comph) %>% 
      reframe(lm(dd ~ press, data = cur_data()) %>% broom::tidy()) %>% 
      filter(term == "press") %>% 
      mutate(
        estimate  = estimate  %>% signif(3),
        std.error = std.error %>% signif(3)
      ),
    aes(
      x = xloc,
      label = str_glue("{estimate}\n+/- {std.error} Å/bar")
    ),
    y = -0.5,
    check_overlap = TRUE
  ) +
  scale_shape_manual(values = c(up = 24, dn = 25)) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(
    title = "Pressure effect on lamellar repeat spacing, T = 20 deg C",
    x = "Pressure (bar)",
    y = "Delta d (Å)"
  )
ggsave(here("04-omitfigs", "deltad_20C_20230808.pdf"), width = 6, height = 4)

# data dump for Itay
peaqs_lam %>% 
  ungroup() %>% 
  filter((comp == "18:1/18:1 PE") & (ph == "la")) %>% 
  select(press, pdir, d, dd) %>% 
  write_tsv(here("02-datadump", "deltad_DOPE.tsv"))
