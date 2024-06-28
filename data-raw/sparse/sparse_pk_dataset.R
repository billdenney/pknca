library(tidyverse)
library(rio)
library(PKNCA)

d_raw <- import("sparse_pk_dataset.xlsx")

d_nca_source <-
  d_raw %>%
  select(
    `Source filename`, DOI, `Drug name`, `Dose Route`, `Dose unit`, Specimen,
    AUC, `AUC unit`, `AUC SE`, `AUC degrees of freedom`
  ) %>%
  filter(!is.na(AUC))

d_conc <-
  d_raw %>%
  select(
    `Source filename`, DOI, `Drug name`, `Dose Route`, `Dose unit`, Specimen,
    `Animal ID`, Time, `Time unit`, Concentration, `Concentration unit`
  ) %>%
  filter(!is.na(Concentration))

stopifnot(all(is.na(d_conc$`Animal ID`)))
d_conc$`Animal ID` <- seq_len(nrow(d_conc))

o_conc <-
  d_conc %>%
  PKNCAconc(Concentration~Time|`Source filename`+DOI+`Drug name`+Specimen+`Animal ID`, sparse=TRUE)
d_intervals <-
  d_conc %>%
  group_by(`Source filename`, DOI, `Drug name`, `Dose Route`, `Dose unit`, Specimen) %>%
  summarize(
    start=min(Time),
    end=max(Time),
    sparse_auclast=TRUE
  )
o_data <-
  PKNCAdata(o_conc, intervals=d_intervals)
o_nca <- pk.nca(o_data)
summary(o_nca, drop.group = character())
as.data.frame(as.data.frame(o_nca))

