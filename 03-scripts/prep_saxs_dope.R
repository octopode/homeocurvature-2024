## Load in SAXS data from temperature sweeps on E. coli polar lipids
## and cross-reference them by timestamp to the temperature log.

source(here("03-scripts", "saxs_helpers.R"))

# identify the samples
samp2pheno = c(
  "JWL0311" = "PPE+",
  "JWL0312" = "PPE-",
  "JWL0315A" = "PE",
  "JWL0316A" = "PC"
)

# location of the SAXS profiles
dir_profiles = here("01-rawdata", "saxs_dope")

# load all the profiles
saxsprof_dope = list.files(
  path = dir_profiles,
  pattern = "*\\.dat",
  full.names = TRUE
) %>% 
  #.[1:3] %>% 
  read_saxsall() %>% 
  # convert pressure to bar
  mutate(press = 10*press)

# Panel S6C: DOPE profiles from 700-1500 bar
saxsprof_dope %>% 
  filter(between(press, 700, 1500)) %>% 
  # average up- and downsweeps
  group_by(temp, press, q) %>% 
  summarise(iq = mean(iq)) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      color = press,
      group = press
    )
  ) +
  ## or view sweeps separately
  #facet_wrap(~pdir) +
  geom_line() +
  theme_pubr() +
  scale_y_log10()

# Panel S6C: DOPE profiles from 700-1500 bar
# offset version
saxsprof_dope %>% 
  ungroup() %>% 
  filter(between(press, 700, 1500)) %>% 
  filter(between(q, 0.05, 0.5)) %>% 
  mutate(
    # transform
    iq = log10(iq),
    # offset!
    iq = iq + 1.7*(press %>% as.factor() %>% as.integer())
  ) %>% 
  # average up- and downsweeps
  group_by(temp, press, q) %>% 
  summarise(iq = mean(iq)) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      #color = press,
      group = press
    )
  ) +
  ## or view sweeps separately
  #facet_wrap(~pdir) +
  geom_line() +
  #geom_text(
  #  x = 0.55,
  #  nudge_y = 0.5,
  #  check_overlap = TRUE,
  #  aes(label = str_glue("{press} bar"))
  #) +
  theme_pubr() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  xlim(c(0, 0.6)) +
  labs(
    x = "q (Å^-1)",
    y = "I(q) (arb. log units)"
  )
ggsave()