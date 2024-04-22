## Do some basic prep and grooming on the raw data, including
## PPE standard correction
## metadata joining
## subsetting (PLs, TLs, etc.) and renormalization
## sterol GCMS data joining and renorm
## all model fitting happens in plot_lipidmaps.R

## CHANGELOG
# originated 20230706 JRW
# added sterol data cross-ref 20230709 JRW

library(tidyverse)
library(here)
library(httr) # to load the metadata table

source(here("03-scripts", "seelipids_helpers.R"))
#source(here("03-scripts", "parse_lipidmaps.R")) # load all the way from XLSX

# QC list of samples that failed (clearly contaminated or too dilute)
blacklist = c(
  "JWL0169"
)

# load the parsed LCMS data
file_tidy = here("02-tidydata", "lipidmaps_raw.tsv")
# skip if parsing was just done and lmapdata_long is already in the global environment
lmapdata_raw = read_tsv(file_tidy)

# this file contains corrections for erroneous response factors in pre-2022 data
file_corr = here("01-rawdata", "20220609_lcms_corrfacs.tsv")

# this file joins EID (JWL#) to species, SID (e.g. D1234-SS1), and environmental data for ALL samples
file_meta = here("01-rawdata", "lipidomics_metadata.tsv")
# this is a Google Sheets key for the metadata file, to facilitate a nice Excel-like interface
gsht_meta = "1afeFHB3qRPNarH26FyAh5yZXU_yXlKS6p4XrFdmGjLg"

# sterol concentrations prenormalized to nmol/µg (consistent with Lipid Maps)
# as of 20230709 upstream code is in the projects cteno-cholesterol (=tidychrom)
# and cteno-pl-transnorm (=seelipids). Could include here, but none of the other
# chrom traces are included.
# Because I did the integration, this tsv is already tidy.
file_ster = here("01-rawdata", "gcms", "sterol_molconcs_merge.tsv")

# make JWL extract ID compliant
# funky braces! per https://stackoverflow.com/questions/58575581/use-stringrstr-glue-with-pipes
cln_jwl = function(eid){
  # allow one letter at the end
  suff = str_extract(eid, regex("[A-Z]$")) %>% 
    ifelse(is.na(.), '', .)
  eid %>% 
    parse_number() %>% 
    str_pad(4, pad='0') %>% 
    {str_glue("JWL{.}{suff}")} %>% 
    as.character()
}

lmapdata_cln = lmapdata_raw %>% 
  mutate(
    # ensure eids are correctly padded
    eid = eid %>% cln_jwl(),
    # rename plasmalogens for ease of labeling
    class = class %>% str_remove('-'),
    # linking key for pre-2022 files
    orig = str_detect(fname, "2020")
  ) %>% 
  # then make compound ID a factor so it plots in a uniform order
  #NTS 20230706: figure out whether dividing the DG, TG abundances by 10 is just a pre-2022 thing.
  # apply Ossie's headgroup response factor corrections to the pre-2022 data
  left_join(
    read_tsv(file_corr) %>% mutate(orig = TRUE),
    by = c("class", "orig")
  ) %>% 
  replace_na(list("corrfac" = 1)) %>% 
  mutate(rab = rab*corrfac) %>%
  # then, multiply by standard molar masses
  # to get absolute abundances (aab) in nmol/sample
  #left_join(file_molm %>% read_tsv(), by = "class") %>% 
  # transforming to put mw in units of (ng/nmol)
  # aab [10 ng/sample] / MW [ng/nmol] = aab [10 nmol/sample] => 
  # divide by 10 to get nmol/sample
  #mutate(aab = (rab / mw)/10) %>% 
  mutate(aab = rab) %>% 
  group_by(eid) %>%
  factorize_lipids()

# refresh metadata from online interface
GET(str_glue("https://docs.google.com/spreadsheet/ccc?key={gsht_meta}&output=csv")) %>% 
  content("raw") %>% 
  write_file(file_meta)

# join metadata to lipid data
lmapdata = lmapdata_cln %>%
  ungroup() %>% 
  left_join(file_meta %>% read_csv(na = c("", "NA", "#N/A")), by = "eid") %>% 
  replace_na(list(mass = 1)) %>% 
  # normalize all aab values to *ratio per sample*
  mutate(aab = aab*mass) %>% 
  # sort the rows for easy viewing
  # throws a *harmless* error, which is an R or Rstudio bug
  # (Error in exists(cacheKey, where = .rs.WorkingDataEnv, inherits = FALSE))
  arrange(eid, id) %>% 
  filter(!(eid %in% blacklist)) %>% 
  # normalize/introduce 'frac_molar'
  # though note that at this point the denominators will be different:
  # some samples are comprehensive glycerolipid, some are PL only.
  mutate(
    frac_molar = aab,
    frac_molar = frac_molar/sum(frac_molar),
    # 20230709 think I should leave aab alone here
    # also calculate per-chain dbi
    dbi = dbonds/tails,
    chn = carbon/tails
  )

# check for missing metadata entries
eids_missing = lmapdata %>% 
  filter(is.na(sp)) %>% 
  .$eid %>% 
  unique()
if(length(eids_missing > 0)){
  message("Metadata entries are missing for:")
  print(eids_missing)
}

# store it
lmapdata %>% write_tsv(here("02-tidydata", "lipidmaps_wmeta.tsv"))

# filter for PLs and renormalize each sample
pldata = lmapdata %>% 
  filter(
    str_detect(as.character(class), 'P') &
      !str_detect(as.character(class), 'Cer')
  ) %>% 
  group_by(eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # every sample was analyzed for PLs, so I can make zeroes explicit
  # this guards against some pitfalls down the line when calculating mean c0,
  # DBI, mean chain length etc.
  ungroup() %>% 
  complete(
    # depth, temp, other predictors NEED to be here, or trouble down the line!
    nesting(sp, eid, sid, tissue, wild, depth_col, temp_col, lat, lon, remarks), 
    # do NOT include id here; it includes RT
    nesting(class, annot, carbon, dbonds, dbi, tails, carbsn1, dbonsn1, carbsn2, dbonsn2), 
    fill = list(frac_molar = 0, aab = 0)
  )

# just whole and body samples
pldata_wildwhole = pldata %>% 
  # filter out tech reps for downstream analyses
  filter(!str_detect(eid, regex("[ABC]"))) %>% 
  filter(wild & (tissue %in% c("whole", "body")))

# QC the E. coli data and average by strain
pldata_ecoli_indls = pldata %>% 
  filter(sp == "Esch_coli") %>%
  # We can do some annotation QC based on genotypes
  # use mutate(ifelse()) instead of filter(), which is for some reason messing up the stack order
  mutate(
    frac_molar = ifelse((sid != "AAL95 + pAC (PC)") & str_detect(class, "PC"), 0, frac_molar),
    frac_molar = ifelse(!str_detect(sid, "pPLsCP") & (class == "PPE"), 0, frac_molar)
  ) %>% 
  # there isn't any anyway, so this gets rid of the legend item
  filter(class != "PPC") %>% 
  # renorm
  group_by(eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # 1 bar cultures only; might wanna remove
  filter(depth_col == 0)

pldata_ecoli = pldata_ecoli_indls %>% 
  # average compounds across replicates
  replace_na(list("carbsn2" = 0, "dbonsn2" = 0)) %>% 
  group_by(sid, class, annot, carbon, carbsn1, carbsn2, dbonds, dbonsn1, dbonsn2) %>% 
  summarise(
    frac_molar = mean(frac_molar)#,
    #eid = eid[[1]] #hack
  ) %>% 
  # and renorm
  group_by(sid) %>%
  mutate(
    frac_molar = frac_molar/sum(frac_molar)#,
    #eid = eid[[1]]
  ) %>% 
  arrange(sid, class, carbon, dbonds, annot) %>% 
  group_by(sid)

## exploring
#pldata_wildwhole %>% 
#  filter(parse_number(eid) > 300) %>% 
#  #filter(aab <500) %>% 
#  ggplot(aes(x=aab)) + geom_histogram() + facet_wrap(~eid)#scale_y_log10()

## Here, join GCMS sterol data and norm to (PLs + sterols)
plstdata = pldata %>% 
  bind_rows(
    read_tsv(file_ster) %>% 
      select(eid, class, id, aab) %>% 
      mutate(annot = id)# %>% 
      #mutate(aab = aab/1E5)
  ) %>% 
  # only samples for which both sterols and PL/GLs were measured
  group_by(eid) %>% 
  filter(max(str_detect(class, "ST")) %>% as.logical()) %>% 
  filter(max(str_detect(class, "PC")) %>% as.logical()) %>%
  arrange(sp) %>% 
  #mutate(sp = sort(sp)[[1]]) %>% 
  mutate(across(c(sp, depth_col, temp_col), function(x){sort(x)[[1]]})) %>% 
  #left_join(file_meta %>% read_csv(na = c("", "NA", "#N/A")), by = c("eid", "sp")) %>% 
  # renorm
  mutate(frac_molar = aab/sum(aab, na.rm = TRUE)) %>% 
  # fill in the 'wild' and 'tissue' columns
  group_by(eid) %>% 
  arrange(wild, tissue) %>% 
  mutate(
    wild = wild %>% .[which(!is.na(.))] %>% .[[1]],
    tissue = tissue %>% .[which(!is.na(.))] %>% .[[1]]
  ) %>% 
  ungroup() %>% 
  factorize_lipids()

## having a look at dist of chol
#plstdata %>% 
#  filter(id == "chol") %>% 
#  ggplot(aes(
#    x = eid,
#    y = frac_molar
#  )) +
#  geom_col() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
## low EIDs show very large values; aab units must be different
## could be factor of 10
## could be consequence of correction factors!

# and here's a potential useful dataframe of species means
# binding the (smaller) sterol dataset to all available PL data by species
# why is this so goshdarn slow?! -20230926
plstsp = pldata_wildwhole %>% 
  ## sum within classes
  #group_by(eid, sp, class) %>% 
  #summarise(frac_molar = sum(frac_molar)) %>% 
  filter(eid != "JWL0171") %>% # 20231223 for the MD systems plot, since
  # that sample is degraded
  # average across individuals
  # do NOT group by id; it includes RT
  group_by(sp, class, annot, carbon, dbonds) %>% 
  summarise(frac_molar = mean(frac_molar)) %>% 
  # now bind the cholesterol data
  bind_rows(
    plstdata %>% 
      filter(class == "ST") %>% 
      group_by(sp, class, id, annot) %>% 
      summarize(frac_molar = mean(frac_molar)) %>% 
      # propose same for Boli_micr as for Boli_infu
      bind_rows(
        {.} %>% filter(sp == "Boli_infu") %>% 
            mutate(
              sp = "Boli_micr",
              annot = "chol?"
            )
      )
  ) %>% 
  group_by(sp) %>% 
  # drop all spp lacking sterol data
  filter("ST" %in% class) %>% 
  # shrink everything else to fit the niche left by cholesterol
  mutate(
    frac_molar = ifelse(
      class == "ST", 
      frac_molar, 
      frac_molar * (1 - (cur_data() %>% filter(class == "ST") %>% .$frac_molar %>% sum()))
    )
  )

# and here's a potential useful dataframe of species means
# binding the (smaller) sterol dataset to all available PL data by species
plstspmean = pldata_wildwhole %>% 
  # sum within classes
  group_by(eid, sp, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  # average across individuals
  group_by(sp, class) %>% 
  summarise(
    frac_molar = mean(frac_molar),
    n = {str_glue("N = {n()}")}
  ) %>% 
  # now bind the cholesterol data
  bind_rows(
    plstdata %>% 
      filter(class == "ST") %>% 
      group_by(sp, class) %>% 
      summarize(
        frac_molar = mean(frac_molar),
        n = {str_glue("N = {n()}")}
      )
  ) %>% 
  # shrink everything else to fit the niche left by cholesterol
  mutate(
    frac_molar = ifelse(
      class == "ST", 
      frac_molar, 
      frac_molar * (1 - (cur_data() %>% filter(class == "ST") %>% .$frac_molar %>% sum()))
    )
  ) %>% 
  # drop all spp lacking sterol data
  filter("ST" %in% class)

## PLs and sphingolipids for all available samples
plsldata = lmapdata %>% 
  filter(wild & (tissue %in% c("whole", "body"))) %>% 
  filter(!(class %in% c("DG", "TG", "AC", "CE"))) %>% 
  group_by(sp, eid) %>% 
  filter(max(str_detect(class, "Cer")) %>% as.logical()) %>% 
  # renorm
  mutate(frac_molar = aab/sum(aab, na.rm = TRUE)) %>% 
  factorize_lipids()

## PLs and just sphingomyelin for all available samples
plsmdata = lmapdata %>% 
  filter(wild & (tissue %in% c("whole", "body"))) %>% 
  filter(!(class %in% c("DG", "TG", "AC", "CE"))) %>% 
  group_by(sp, eid) %>% 
  # detect SL analysis; ensures the zeroes get thru
  filter(max(str_detect(class, "Cer")) %>% as.logical()) %>% 
  # just SM!
  filter(!str_detect(class, "Cer")) %>% 
  # renorm
  mutate(frac_molar = aab/sum(aab, na.rm = TRUE)) %>% 
  factorize_lipids()

# PLs, sterols and sphingos
plstldata = plstdata %>% 
  bind_rows(
    plsldata %>%
      filter(class %in% c("Cer", "Cer1P", "SM"))
  ) %>% 
  group_by(eid) %>% 
  filter(("ST" %in% class) & ("Cer" %in% class)) %>% 
  # renorm
  group_by(sp, eid) %>% 
  mutate(frac_molar = aab/sum(aab, na.rm = TRUE)) %>% 
  factorize_lipids()

# spit out Table S1
pldata_wildwhole %>% 
  select(eid, sp, depth_col, temp_col, lat, lon) %>% 
  distinct() %>% 
  drop_na() %>% 
  # arrange them in a particular way
  group_by(sp) %>% 
  mutate(
    across(contains("_col"), list("mean" = mean)),
    shal = depth_col_mean <= 250,
    cold = temp_col_mean <= 10
  ) %>% 
  mutate(temp_col = round(temp_col, 1)) %>% 
  arrange(-shal, -temp_col_mean, -cold, depth_col_mean, -temp_col, depth_col, eid) %>% 
  mutate(eid = str_remove(eid, "JWL")) %>% 
  select(-contains("_mean"), -cold, -shal) %>% 
  print(n = 66) %>% 
  write_tsv(here("04-suppfigs", "TabS1.tsv"))

# get mean z and T
ctenodata_summ = pldata_wildwhole %>% 
  group_by(eid, depth_col, temp_col, sp, class) %>% 
  summarise(frac_molar = sum(frac_molar)) %>% 
  group_by(sp, class) %>% 
  summarise(
    n = n(),
    across(contains("_col"), c("min"=min, "max"=max, "mean"=mean)),
    frac_molar = mean(frac_molar)
  ) %>% 
  # *ensure* labels all agree
  group_by(sp) %>% 
  mutate(
    n = max(n),
    across(contains("_min"), min),
    across(contains("_max"), max),
    across(contains("_col_"), round), # to integer?
    lab = str_glue("N = {n}\n{depth_col_min} - {depth_col_max} m\n{temp_col_min} - {temp_col_max}°C"),
    # labels with just the means
    mlab = str_glue("{depth_col_mean} m\n{temp_col_mean}°C")
  ) %>% 
  # stretch lipid data out to compress the tree (which has width 1)
  mutate(frac_molar = ratio_pheno_phylo*frac_molar) %>% 
  # hack to shift the whole panel to the right
  bind_rows(
    .,
    tibble(sp = .$sp %>% unique()) %>% 
      mutate(
        class = "spacer",
        frac_molar = pheno_start_x # the offset
      )
  ) %>% 
  # as last level, it plots on the left
  mutate(class = class %>% factor(levels = c(levels(pldata$class), "spacer")))

ctenodata_summ %>% 
  filter(
    sp %in% c(
      "Boli_vitr",
      "Leuc_pulc",
      "Boli_micr",
      "Boli_infu",
      "Lamp_crue",
      "Tjal_pink"
    )
  ) %>% 
  select(contains("col_mean")) %>% 
  distinct() %>% 
  mutate(readout = str_glue("{depth_col_mean} m • {temp_col_mean}°C"))

# revised MD systems for Tjal, Lampo, B. vitrea
# 20231011
# adding Boli micr and Batho 20231223
mdsys_revised = pldata_wildwhole %>% 
  filter(
    sp %in% c(
      "Boli_vitr",
      "Boli_infu",
      "Lamp_crue",
      "Tjal_pink",
      "Boli_micr",
      "Bath_fost"
    )
  ) %>% 
  #that one bad Batho
  filter(eid != "JWL0171") %>% 
  # sum lipid species
  group_by(sp, eid, class) %>% 
  arrange(-frac_molar) %>% 
  summarize(
    annot = annot[[1]],
    frac_molar = sum(frac_molar)
  ) %>% 
  # average indls
  group_by(sp, class) %>% 
  arrange(-frac_molar) %>% 
  summarize(
    frac_molar = mean(frac_molar),
    annot = annot[[1]]
  ) %>% 
  # renorm
  group_by(sp) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # at least 5%! But always include PS.
  filter((class == "PS") | (frac_molar > 0.05)) %>% 
  filter(!((sp == "Boli_infu") & (class %in% c("PI")))) %>% # hack to keep # components to 7 or less
  filter(!((sp == "Boli_micr") & (class %in% c("LPE")))) %>% # hack to keep # components to 7 or less
  # no more than 6 classes plus cholesterol
  # scale to 90% if Boli_vitr; else 95%
  mutate(
    frac_molar = ifelse(
      sp == "Boli_vitr", 
      0.90*frac_molar/sum(frac_molar),
      0.95*frac_molar/sum(frac_molar)
    )
  ) %>% 
  bind_rows(crossing(class = "ST", annot = "chol", sp = {.$sp %>% unique()})) %>% 
  mutate(
    frac_molar = ifelse(class == "ST",
                          ifelse(sp == "Boli_vitr",
                            0.10,
                            0.05
                          ),
                          frac_molar
                        )
  ) %>% 
  factorize_lipids() %>% 
  # deduce chains on unfragmented PPCs
  mutate(
    annot = ifelse(annot == "39:5", "P-19:0/20:5", annot),
    annot = ifelse(annot == "39:2", "19:0/20:2", annot), # for Bath_fost
    annot = ifelse(annot == "38:5", "P-18:0/20:5", annot),
    annot = ifelse(annot == "36:7", "P-14:0/22:6", annot),
    # idk why this is happening, but the hack does not corrupt anything
    annot = ifelse(annot == "22:6/19:0", "18:0/22:6", annot),
  )
  
# store it
mdsys_revised %>% 
  arrange(sp, class) %>% 
  print() %>% 
  write_tsv(here("01-rawdata", "mdsystems_complex_revised_20231223.tsv"))

## DATA S1: THE BIG LIPIDOMIC TABLE
# Metadata: eid, sp, depth, temp, tissue, remarks
# Data columns: each lipid species, arranged by class then chain.
# Normalized to total PLs. Chol and SM are normalized to PLs + chol or SM.
# NAs where chol or SM was not measured.

# a list to expand the codenames to spnames used in the paper
spname_aliases = c(
  "Bath_fost" = "Bathocyroe aff. fosteri", 
  "Bero_abys" = "Beroe abyssicola", 
  "Bero_cucu" = "Beroe cucumis", 
  "Bero_ovat" = "Beroe ovata", 
  "Bero_pseu" = "Beroe pseudocucumis", 
  "Boli_infu" = "Beroe infundibulum",
  "Boli_micr" = "Bolinopsis microptera",
  "Boli_vitr" = "Bolinopsis vitrea", 
  "Cest_vene" = "Cestum veneris", 
  "Coel_hali" = "Coeloplana sp.", 
  "Cydi_blac" = "Mertensiidae sp. K",
  "Lamp_crue" = "Lampocteis cruentiventer",
  "Leuc_pulc" = "Leucothea pulchra", 
  "Llyr_bent" = "Lobata sp. L1", 
  "Llyr_deep" = "Lobata sp. L2", 
  "Mert_angu" = "Mertensiidae sp. T",
  "Tjal_pink" = "Platyctenida sp. T"
  )

# this code from prep_pgls_lipidmaps.R
# varname is changed
# summarize some stuff ind'lwise
# we want total frac_molar, sn1, sn2 chain and dbi within each class
dataS1 = pldata %>% 
  # but no E. coli or captive specimens
  filter(wild) %>% 
  # bind SM data, but leave PL normalization for everything else
  bind_rows(
    plsmdata %>% 
      filter(class == "SM")
  ) %>% 
  # bind sterol data, but leave PL normalization for everything else
  bind_rows(
    plstdata %>% 
      filter(class == "ST")# %>% 
      # JWL0185 is a tentacle sample with sterol data!
      #filter(wild & (tissue %in% c("whole", "body")))
  ) %>% 
  # make tissue ordered
  mutate(subsample = factor(tissue, levels = c("whole", "body", "tentacle"))) %>% 
  # fill missing metadata for SM, chol rows
  group_by(sp, eid) %>% 
  fill(depth_col, temp_col, lat, lon, subsample, remarks) %>% 
  # give the lipids names without RTs and order them
  factorize_lipids() %>% 
  #arrange(id) %>% 
  mutate(
    lipidsp = paste(annot, class),
    # just rename cholesterol
    lipidsp = ifelse(lipidsp == "chol ST", "chol", lipidsp),
    lipidsp = lipidsp %>% factor(., levels = unique(.)),
  ) %>% 
  # strip out lipid species that are just empty factor levels
  group_by(lipidsp) %>% 
  filter(max(frac_molar) > 0) %>% 
  # metadata columns we don't want to pivot
  group_by(sp, eid, depth_col, temp_col, lat, lon, subsample, remarks) %>%
  # check that PLs are properly normalized
  #filter(str_detect(lipidsp, 'P')) %>% 
  #summarise(frac_molar = sum(frac_molar)) # looks right!
  # group columns get left in, of course
  select(lipidsp, frac_molar) %>% 
  pivot_wider(
    id_cols = group_cols(), 
    names_from = lipidsp, 
    values_from = frac_molar,
    values_fn = sum # turns out there are some isobaric lipid spp with different RTs; sum them
  ) %>% 
  # arrange in a way that matches Table S1
  group_by(sp) %>% 
  mutate(
    across(contains("_col"), list("mean" = mean)),
    shal = depth_col_mean <= 250,
    cold = temp_col_mean <= 10
  ) %>% 
  mutate(temp_col = round(temp_col, 1)) %>% 
  arrange(-shal, -temp_col_mean, -cold, depth_col_mean, subsample, -temp_col, depth_col, eid) %>% 
  mutate(eid = str_remove(eid, "JWL")) %>% 
  select(-contains("_mean"), -cold, -shal) %>% #.$sp %>% unique() %>% dput()
  # translate the spnames
  mutate(sp = spname_aliases[sp]) %>% 
  print()

dataS1 %>% 
  write_tsv(here("04-suppfigs", "DataS1_20231016.tsv"))
  
# 20231228 what is this?
#  
#
## sort the rows
#  arrange(sp_placeholder, tissue)
#
#bind_rows(
#  
#)
#  