# code to prepare `sdtm_paramcd` dataset goes here
library(tidyverse)
library(janitor)
library(assertr)

ld_sdtm_raw <- rio::import_list("https://evs.nci.nih.gov/ftp1/CDISC/SDTM/SDTM%20Terminology.xls")

idx_ld_sdtm_raw <- which(startsWith(names(ld_sdtm_raw), "SDTM"))
stopifnot(length(idx_ld_sdtm_raw) == 1)
d_sdtm_term <-
  ld_sdtm_raw[[idx_ld_sdtm_raw]] |>
  clean_names()

choose_pptest <- function(x) {
  l_str <- strsplit(x = x, split = "; *")
  # Choose the longest description when there are multiple
  sapply(X = l_str, FUN = \(x) x[nchar(x) == max(nchar(x))][1])
}

# d_sdtm_param <-
#   d_sdtm_term |>
#   filter(codelist_name %in% "PK Parameters") |>
#   filter(!is.na(codelist_code)) |>
#   rename(
#     PPTEST = cdisc_submission_value
#   )
sdtm_paramcd <-
  d_sdtm_term |>
  filter(codelist_name %in% "PK Parameters Code") |>
  filter(!is.na(codelist_code)) |>
  mutate(
    PPTEST = choose_pptest(cdisc_synonym_s)
  ) |>
  select(
    PPTESTCD = cdisc_submission_value,
    PPTEST,
    CDISC_definition = cdisc_definition
  ) |>
  arrange(PPTESTCD) |>
  verify(!duplicated(PPTESTCD)) |>
  verify(!duplicated(PPTEST))

usethis::use_data(sdtm_paramcd, overwrite = TRUE)

sdtm_pkunit <-
  d_sdtm_term |>
  filter(codelist_name %in% "PK Units of Measure") |>
  filter(!is.na(codelist_code)) |>
  select(
    submission = cdisc_submission_value,
    synonym = cdisc_synonym_s
  ) |>
  mutate(
    synonym = strsplit(x = synonym, split = "; *")
  ) |>
  verify(!duplicated(submission)) |>
  group_by(submission) |>
  reframe(synonym = unlist(synonym)) |>
  ungroup()

usethis::use_data(sdtm_pkunit, overwrite = TRUE)
