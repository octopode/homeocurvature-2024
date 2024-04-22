## Plot grids of E. coli SAXS profiles arranged by temperature

library(httr)

source(here("03-scripts", "plot_helpers.R"))
source(here("03-scripts", "prep_saxs_ecoli.R"))

# CHANGING THIS STUFF WILL BREAK THE PEAK ANNOTATOR BELOW. Fork the script first!
# q-window to crop to
qmin = 0.05
qmax = 0.35 # short of the Kapton peak, also easy to see compact plot

# temperature bin width for averaging profiles
binwid = 10/3

# set temp intervals to be used for the plots
# exceptions can be made below
tstart = 0
tstep  = 10

# average the scattering data into temperature bins
# and then select the ones we want to display
# it's pretty fast so simple to do all at once
profs_tbins = profs_wtemp %>%
  # crop profiles (and speed things up a lot!)
  filter(between(q, qmin, qmax)) %>%
  # make pheno a factor for plotting order
  #mutate(pheno = pheno %>% factor(levels = c("PPE-", "PPE+", "PE", "PC"))) %>% 
  mutate(pheno = pheno %>% factor(levels = c("PPE+", "PC", "PPE-", "PE"))) %>% 
  # round to integer
  mutate(temp = binwid*round(temp/binwid)) %>% 
  group_by(temp, samp, pheno, q) %>% 
  summarize(iq = mean(iq)) %>% 
  arrange(samp, temp) %>% 
  # filter to the intervals we want to work with
  filter(
    # get t-intervals
    !((temp-tstart) %% tstep) |
      # hack to get in a PC profile < 20°C!
      ((pheno == "PC") & (temp %in% c(40/3)))
  ) %>% 
  # round remaining temps to integer
  # this avoids joining on values with repeating decimals
  # and gives the correct # sig figs (the bins are >1°C wide)
  mutate(temp = round(temp)) %>% 
  # reverse temperatures so they rise from the bottom in a facet
  ungroup() %>% 
  arrange(-temp) %>% 
  mutate(tempfac = temp %>% factor(levels = unique(temp)))

# "Offset plots" instead of a facet, to save some time in Illustrator
# to do this, I need to log-transform all the I(q)s
# and specify a vertical offset...in log(I(q)) units per °C
offs = 0.1
profs_offs = profs_tbins %>% 
  # derive plotting x-values
  # use temp here (increases in right direction) rather than tempfac
  mutate(hiq = log(iq) + offs*temp) %>% 
  arrange(samp, temp, q) 

#profs_offs %>% 
#  ggplot(
#    aes(
#      x = q,
#      y = hiq,
#      group = temp
#    )
#  ) +
#  facet_grid(col = vars(pheno)) +
#  geom_line(
#    aes(color = pheno)
#  ) +
#  # automatic labels!
#  geom_text(
#    data = profs_offs %>% 
#      group_by(samp, temp) %>% 
#      arrange(q) %>% 
#      slice(1),
#    x = 0.04,
#    aes(label = str_glue("{temp}°C")),
#    color = "black",
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  scale_color_manual(values=chroma_pheno) +
#  theme_pubr() +
#  guides(color = "none", y = "none") +
#  lims(x = c(0, qmax)) +
#  labs(
#    title = "E. coli polar lipids: SAXS profiles along temperature ramp",
#    x = "q (1/Å)",
#    y = "I(q) (arb. log units)"
#  )
#ggsave(file = here("04-mainfigs", "Ecoli_all_offset_20230726a.pdf"), width = 10, height = 4)

# Now gonna try a peak assignment workflow...
# defining characteristic spacings
# a nice symbolic guide to q-spacings in review by Seddon, BBA 1990
pkspacings = bind_rows(
  # for hii
  tibble(
    ord = seq(7),
    spc = c(1, 3, 4, 7, 9, 12, 13) %>% sqrt()
  ) %>% 
    mutate(ph = "hii"),
  # and lamellar phases
  tibble(
    ord = seq(7),
    spc = ord
  ) %>% 
    mutate(ph = "la"),
  tibble(
    ord = seq(7),
    spc = ord
  ) %>% 
    mutate(ph = "lb")
) %>% 
  group_by(ph) %>% 
  arrange(ord) %>% 
  filter(row_number() <= 4) # only to 4th-order

# pick top n peaks (oversensitive top-down algorithm)
# DO NOT CHANGE ANY OF THESE CRITERIA without forking the entire script!
# It will break downstream analyses!
npk = 6
wid = 5 # number of points across which to average for simple smoothing
pksall = profs_tbins %>% 
  # crop
  filter(between(q, qmin, qmax)) %>% 
  # do the q-bin smoothing
  group_by(pheno, temp, tempfac) %>% 
  arrange(q) %>% 
  mutate(binnum = as.integer((row_number()-1)/wid)) %>%
  group_by(pheno, temp, tempfac, binnum) %>% 
  summarise(across(contains('q'), mean)) %>% 
  filter((iq > lead(iq)) & (iq > lag(iq))) %>% 
  arrange(-iq) %>% 
  filter(row_number() <= npk) %>% 
  arrange(pheno, temp, q) %>% 
  mutate(pknum = row_number())
  
## facet showing picked peaks
## good enough. NOTE: doesn't pick up shoulders
#profs_tbins %>% 
#  # crop
#  filter(between(q, qmin, qmax)) %>% 
#  # get t-intervals
#  filter(!((temp-tstart) %% tstep)) %>% 
#  ggplot(
#    aes(
#      x = q,
#      y = iq,
#      color = pheno
#    )
#  ) +
#  facet_grid(col = vars(pheno), row = vars(tempfac)) +
#  geom_line() +
#  geom_point(
#    shape = 6,
#    color = "red",
#    data = pksall,
#    aes(y = iq*1.5)
#  ) +
#  scale_y_log10() +
#  scale_color_manual(values=chroma_pheno) +
#  theme_pubr() +
#  guides(color = "none") +
#  labs(
#    title = "E. coli -/+ PPE:\nprofiles along temperature ramp",
#    x = "q (1/Å)",
#    y = "I(q) (arb. units)"
#  ) +
#  lims(x = c(0, 0.35))

# Make peak assignments here using a GSheets interface.
# As of 20230725, this allows the assignment of one peak per phase per profile,
# and the others get filled in.
# refresh c0 spreadsheet from online interface
gsht_assgts = "1GvZuPtJSCouGXZarkgNTjZVF258DWOAa2s0OgoGaIac"
file_assgts = here("01-rawdata", "pkassign_ecoli.csv")
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_assgts}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_assgts)
# load it back from the disk
assgts = file_assgts %>% 
  read_csv() %>% 
  filter(show)

# fill in the other peaks for each profile
pksfit = assgts %>% 
  left_join(pkspacings, by = "ph", suffix = c('_assg', '_pred'), relationship = "many-to-many") %>% 
  mutate(pknum = ifelse(ord_assg == ord_pred, pknum, NA)) %>% 
  # predict other peak positions
  left_join(pksall) %>% # avoid float matching trouble
  group_by(pheno, temp, ph) %>% 
  mutate(
    basis = q/spc,
    # there should only be one known val per group
    q = mean(basis, na.rm = TRUE)*spc,
    tempfac = first(tempfac, na_rm = TRUE)
  ) %>% 
  # make pheno a factor for plotting order
  mutate(pheno = pheno %>% factor(levels = c("PPE-", "PPE+", "PE", "PC")))

## now visualize those predictions
#profs_tbins %>% 
#  # crop
#  filter(between(q, qmin, qmax)) %>% 
#  ggplot(
#    aes(
#      x = q,
#      y = iq,
#      color = pheno
#    )
#  ) +
#  facet_grid(col = vars(pheno), row = vars(tempfac)) +
#  geom_line() +
#  # can show all picked peaks
#  geom_point(
#    shape = 6,
#    color = "red",
#    data = pksall,
#    aes(y = iq*1.5)
#  ) +
#  geom_vline(
#    data = pksfit,
#    aes(
#      xintercept = q,
#      linetype = ph
#    )
#  ) +
#  scale_y_log10() +
#  scale_color_manual(values=chroma_pheno) +
#  theme_pubr() +
#  guides(color = "none") +
#  labs(
#    title = "E. coli -/+ PPE: profiles along temperature ramp",
#    x = "q (1/Å)",
#    y = "I(q) (arb. units)"
#  ) +
#  lims(x = c(0, 0.35))
##ggsave(file = here("04-omitfigs", "Ecoli_phaselabels_20230725b.pdf"), width = 8, height = 4)

# OK, now bundle it all together: make a nice offset profile plot
# with colored labels done by geom_line
# which will require a funny color scale
chroma_phenoph = c(chroma_pheno, chroma_ph)

# construct a little tick running under the tip of each peak
pkoffs = pksfit %>% 
  select(-iq) %>% 
  left_join(
    profs_offs,
    by = c("pheno", "temp"), 
    relationship = "many-to-many",
    suffix = c("_pk", '')
  ) %>% 
  group_by(pheno, temp, ph, ord_pred, q_pk) %>% 
  summarize(hiqfun = approxfun(q, hiq) %>% list()) %>% 
  rowwise() %>% 
  mutate(
    hiq = hiqfun(q_pk),
    hiq_hi = hiq + 0.25, # how far to go vertically above
    hiq_lo = hiq - 0.25  # and below
  ) %>% 
  select(-hiqfun) %>% 
  pivot_longer(cols = contains("hiq"), names_to = "pttype", values_to = "hiq") %>% 
  group_by(pheno) %>% 
  # to plot vertically w geom_path
  arrange(hiq) %>%
  dplyr::rename(q = q_pk)

#profs_offs %>% 
#  mutate(cmap = pheno) %>% 
#  ggplot(
#    aes(
#      x = q,
#      y = hiq,
#      color = cmap,
#      group = temp
#    )
#  ) +
#  facet_wrap(~pheno, ncol = 2) +
#  # phase peak labels
#  geom_path(
#    linewidth = 2,
#    data = pkoffs %>% mutate(cmap = ph),
#    aes(group = paste(ph, ord_pred))
#  ) +
#  # these are actual profile data
#  geom_line(color = "white", linewidth = 1.25) + # white underlay
#  geom_line() +
#  # automatic labels!
#  geom_text(
#    data = profs_offs %>% 
#      group_by(samp, temp) %>% 
#      arrange(q) %>% 
#      slice(1),
#    x = 0.04,
#    aes(label = str_glue({"{round(temp)}°C"})),
#    color = "black",
#    hjust = 1,
#    check_overlap = TRUE
#  ) +
#  scale_color_manual(values=chroma_phenoph) +
#  theme_pubr() +
#  guides(color = "none", y = "none") +
#  lims(x = c(0, qmax)) +
#  labs(
#    title = "E. coli polar lipids: SAXS profiles along temperature ramp",
#    x = "q (1/Å)",
#    y = "I(q) (arb. log units)"
#  )
#ggsave(file = here("04-mainfigs", "Ecoli_all_offset_peaks_20230726c.pdf"), width = 6, height = 6)

# tiny version for Fig 5C
panel_5c = profs_offs %>% 
  # for color mapping
  mutate(cmap = pheno) %>% 
  ggplot(
    aes(
      x = q,
      y = hiq,
      color = cmap,
      group = temp
    )
  ) +
  facet_wrap(~pheno, ncol = 2) +
  # phase peak labels
  geom_path(
    linewidth = 2/.pt,
    data = pkoffs %>% mutate(cmap = ph),
    aes(group = paste(ph, ord_pred))
  ) +
  # these are actual profile data
  geom_line(color = "white", linewidth = 1/.pt) + # white underlay
  geom_line(size = 0.75/.pt) +
  # automatic labels!
  geom_text(
    data = profs_offs %>% 
      # show "even" temps (mults of 20) only,
      mutate(templab = ifelse(
        !(temp %% 20), 
        str_glue({"{round(temp)}°C"}), 
        NA
      )) %>% 
      group_by(samp, temp) %>% 
      arrange(q) %>% 
      slice(1),
    x = 0.04,
    aes(label = templab),
    color = "black",
    hjust = 1,
    check_overlap = TRUE,
    size = 7/.pt
  ) +
  scale_color_manual(values=chroma_phenoph) +
  theme_tiny() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing.y = unit(-2, "mm"),
    panel.background = element_blank()
  ) +
  guides(color = "none", y = "none") +
  lims(x = c(-0.035, qmax)) +
  labs(
    x = "q (Å^-1)",
    y = "I(q) (arb. log units)"
  )
panel_5c
ggsave(here("04-mainfigs", "Fig5_ecolimodel", "panel_5c_20230802c.pdf"), width = 60, height = 40, units = "mm")

# The Real Panel 5D
# In which I work up some phase diagrams

# Make peak assignments here using a GSheets interface.
# As of 20230725, this allows the assignment of one peak per phase per profile,
# and the others get filled in.
# refresh c0 spreadsheet from online interface
gsht_phcall = "1beN2M_0BXQ_6_ev3qRXwVduVpWMIrpx7jjJlAu2yXs4"
file_phcall = here("01-rawdata", "phcall_ecoli.csv")
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_phcall}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_phcall)
# load it back from the disk
phcall = file_phcall %>% 
  read_csv()

# use the Clausius-Clapeyron values from So et al. 1993 (4.4E-2 +/- 3E-3 °C/bar)
cc_la2hii = 4.4E-2
# and, for now, from Stamatoff et al. 1978, which So cites. But that's for DPPC!
# 48 atm/°C = 48.64 bar/°C = 2.05E-2°C/bar
cc_lb2la  = 2.05E-2
t_range = c(7.5, 82.5)
p_range = c(0, 500)

# just the points where profiles were inspected
pts_panel5d = phcall %>% 
  select(-ph) %>% 
  distinct()

# super super basic:
panel_5d = pts_panel5d %>% 
  bind_rows(
    # stating transition temps here
    tibble(
      pheno = c("PPE+", "PPE+", "PPE-", "PPE-", "PE", "PE", "PC", "PC"),
      txn   = rep(c("lb2la", "la2hii"), 4),
      ccc   = rep(c(cc_lb2la, cc_la2hii), 4), # CC coeffs
      temp  = c(25, 40, 30, 75, 30, 75, 20, NA)
    )
  ) %>% #View()
  mutate(pheno = factor(pheno, levels = c("PPE-", "PE", "PPE+", "PC"))) %>% 
  ggplot(
    aes(
      x = press,
      y = temp,
    )
  ) +
  facet_wrap(~pheno) +
  geom_abline(
    size = 0.5/.pt,
    aes(
      group = txn,
      slope = ccc,
      intercept = temp
    )
  ) +
  geom_point(size = 2/.pt) +
  theme_tiny() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_x_continuous(breaks = c(0, 250)) +
  # handles axis reversal and safe cropping
  coord_flip(
    ylim = c(7.5  , 82.5),
    xlim = c(500,  0)
  ) +
  labs(
    y = "Temperature (deg C)",
    x = "Est. pressure (bar)"
  )
panel_5d
