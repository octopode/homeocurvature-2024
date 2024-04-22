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
dir_profiles = here("01-rawdata", "saxs_dope_dopc")

# load all the profiles
saxsprof_dope_dopc = list.files(
  path = dir_profiles,
  pattern = "*\\.dat",
  full.names = TRUE
) %>% 
  #.[1:3] %>% 
  read_saxsall() %>% 
  # convert pressure to bar
  mutate(press = 10*press)

# Panel S3D: DOPE profiles from 700-1500 bar
saxsprof_dope_dopc %>% 
  #filter(between(press, 700, 1500)) %>% 
  # average up- and downsweeps
  #group_by(temp, press, q) %>% 
  #summarise(iq = mean(iq)) %>% 
  ggplot(
    aes(
      x = q,
      y = iq,
      color = press,
      group = paste(press, pdir)
    )
  ) +
  ## or view sweeps separately
  facet_wrap(~samp) +
  geom_line() +
  theme_pubr() +
  scale_y_log10()

# S3D offset version
saxsprof_dope_dopc %>% 
  filter(between(q, 0.05, 0.5)) %>% 
  filter(!((samp == "DOPC") & (press == 100))) %>% 
  filter(pdir == "up") %>% 
  # make amplitudes uniform
  group_by(samp) %>% 
  mutate(
    iq = iq-min(iq)+0.01,
    iq = iq/max(iq)
  ) %>% #summarise(iq=max(iq))
  # specify where to put the label
  group_by(samp, press) %>% 
  arrange(q) %>% 
  mutate(
    # stack with pressure increasing downward
    ploty = log10(iq) + 0.01 * (-press),
    texty = ifelse(row_number() == 1, ploty, NA)
  ) %>% 
  cross_join(tibble(color = c("wht", "blk") %>% factor(., levels = .))) %>% 
  arrange(samp, press, color) %>% 
  mutate(grp = str_glue("{samp} {press} {color}") %>% factor(., levels = unique(.))) %>% 
  ggplot(
    aes(
      x = q,
      y = ploty,
      group = grp
    )
  ) +
  ## or view sweeps separately
  facet_wrap(~samp) +
  geom_line(
    aes(
      color = color,
      size = color
    )
  ) +
  geom_text(
    aes(
      x = 0.025,
      y = texty,
      label = str_glue("{press} bar")
    ),
    check_overlap = TRUE,
    hjust = 1,
    size = 7/.pt
  ) +
  theme_pubr() +
  scale_color_manual(values = c(wht = "white", blk = "black")) +
  scale_size_manual(values = c(wht = 0.5, blk = 0.25)) +
  scale_x_continuous(
    breaks = seq(0, 0.5, 0.2),
    limits  = c(-0.2, 0.5)
  ) +
  labs(
    x = "q (1/Å)",
    y = "I(q) (arb. log units)"
  ) +
  theme_tiny() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
  #scale_y_log10()
ggsave(here("04-suppfigs", "panel_s3d_20231016a.pdf"), width = 80, height = 80, unit = "mm")
