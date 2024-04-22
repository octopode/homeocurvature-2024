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

# location of the temperature log
file_templog = here("01-rawdata", "templog_20230219.tsv")

# read it in right away
trace_temp = read_tsv(file_templog)

# location of the SAXS profiles
dir_profiles = here("01-rawdata", "saxs_ecoli")

# load all the profiles
saxsprof_ecoli = list.files(
  path = dir_profiles,
  pattern = "*\\.dat",
  full.names = TRUE
) %>% 
  #.[1:3] %>% 
  read_saxsall() %>% 
  # strip out T and P etc.; they'll be added back from the logfile
  select(-temp, -press, -pdir, -rep)

# generate a time -> temp interpolation function
time2temp = approxfun(x = trace_temp$clock, y = trace_temp$temp_ext_a)

# join the temperature to the SAXS data
profs_wtemp = saxsprof_ecoli %>% 
  mutate(temp = time2temp(date)) %>% 
  # and annotate with phenotype
  filter(samp %in% c("JWL0311", "JWL0312", "JWL0315A", "JWL0316A")) %>% 
  mutate(pheno = samp2pheno[samp])
