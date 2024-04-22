library(tidyverse)
library(here)
library(ggpubr)

# simple parser/plotter for MCG lamellar fit results to ground-truth MD

dir_pars = here("02-tidydata", "mcg_lamellar_pars")
# I deleted all the failed DOPE fits so this would work
partabs = list.files(dir_pars, pattern = "*\\.tsv", full.names = TRUE)

pars_mcg = partabs %>% 
  read_tsv() %>% 
  mutate(fname = `@name` %>% basename()) %>% 
  separate(fname, into = c("samp", "press", "pdir", "temp")) %>% 
  # calc P-P distance
  mutate(
    ppd = 2*zH,
    ppd_err = sqrt(2)*zH_err,
    # numericalize pressure
    press = parse_number(press)*10
  )

# for DOPC (JWL0225) these appear to be a little higher than MD vals!
ppd_mcg = pars_mcg %>% select(samp, press, pdir, ppd, ppd_err) %>% 
  arrange(samp, press) %>% 
  mutate(cpd = c(
    "JWL224A" = "DOPE",
    "JWL225A" = "DOPC"
  )[samp])

ppd_md  = bind_rows(
  # DOPE MD data
  tibble(
    press   = seq(0, 1000, 250) %>% c(., 1500),
    ppd     = c(41.30, 41.81, 41.81, 41.74, 41.79, 44.96),
    ppd_err = c( 0.46, 0.44,  0.44,  0.52,  0.49,  1.00)
  ) %>% 
    mutate(cpd = "DOPE", temp = 20) %>% 
    # correction per Pluhackova et al. 2016 Table 4?
    mutate(ppd = ppd/1.08),
  # DOPC MD data
  tibble(
    press   = seq(0, 1000, 250) %>% c(., 1500),
    ppd     = c(38.57, 38.48, 38.52, 38.42, 38.46, 38.69),
    ppd_err = c( 0.40, 0.45,  0.42,  0.36,  0.36,  0.36)
  ) %>% 
    mutate(cpd = "DOPC", temp = 20) %>% 
    # correction per Pluhackova et al. 2016?
    mutate(ppd = ppd/1.015)
)

ppd_all = bind_rows(ppd_mcg, ppd_md) %>% 
  filter(press <= 1700) %>% 
  replace_na(list("pdir" = "md")) %>% 
  mutate(pdir = pdir %>% factor(., levels = c("md", "up", "dn")))

ppd_all %>% 
  #filter(cpd == "DOPE") %>% 
  #filter(press <= 1000) %>% 
  ggplot(
  aes(
    x = press,
    y = ppd
  )
) +
  facet_wrap(~cpd) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(
      color = cpd,
      group = paste(cpd, (pdir == "md"))
    )
  ) +
  geom_errorbar(
    aes(
      ymin = ppd - ppd_err,
      ymax = ppd + ppd_err,
      group = pdir
    ),
    width = 25,
    position = position_dodge2(width = 25, preserve = "single")
  ) +
  geom_point(
    aes(
      fill = cpd,
      shape = pdir
    ),
    position = position_dodge(width = 25)
  ) +
  scale_shape_manual(
    values = c("dn" = 25, "up" = 24, "md" = 21), 
    labels = c("MD", "upsweep", "downsweep")
  ) +
  scale_color_manual(values = chroma_cl[c("PC", "PE")] %>% setNames(c("DOPC", "DOPE"))) +
  scale_fill_manual( values = chroma_cl[c("PC", "PE")] %>% setNames(c("DOPC", "DOPE"))) +
  theme_pubr() +
  guides(color = "none", fill = "none") +
  labs(
    title = "HP MD sims: P-P distance ground-truthing",
    x     = "Pressure (bar) @ 20 deg C",
    y     = "P-P distance (Å)",
    shape  = ''
  )
ggsave(here("04-omitfigs", "ppgroundtruth_DOPE_20230724c.pdf"), width=6, height=4)
