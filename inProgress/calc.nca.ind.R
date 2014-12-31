############################################################
##  This script is for calculating Individual Non-Compartmental PK
##  parameters 
############################################################

source(file.path(pathRunType, "nca.calc.R"))
source(file.path(pathRunType, "auto.pk.days.R"))

## Make sure that all NTAFD == 0 concentrations are 0
if (!filetype.options$allow.nonzero.predose.pk) {
  ## allow missing or 0 pk
  mask.nonzero <- with(d.pk, NTAFD %in% 0 & !(DV %in% c(NA, 0)))
  if (any(mask.nonzero)) {
    nonzero.data <- d.pk[mask.nonzero,c("RANDID", "PERD.name")]
    stop(sprintf("The following subjects have nonzero pk at time 0.  To remove this error, set filetype.options$allow.nonzero.predose.pk=TRUE\n%s",
                 with(nonzero.data,
                      paste("ID", RANDID, PERD.name, collapse="\n"))))
  }
}

## Ensure that the PK data is appropriately sorted.
d.pk <- d.pk[with(d.pk, order(RANDID, PERD, NTAFD)),]

## Determine the times to find the dense PK.
if (given.filename(filenames, "dense.pk.times")) {
  auc.intervals <- normalize.pk.days(filenames$dense.pk.times, d.pk)
} else {
  if (identical(d.dosing, NA)) {
    warning("FOR INFO: Attempting automatic detection of single dose PK.")
    auc.intervals <- auto.pk.days.single(d.pk)
  } else {
    auc.intervals <- auto.pk.days(d.pk, d.dosing)
  }
  ## Add AUC0_24 when AUC0_Inf is requested
  tmp.intervals <- subset(auc.intervals, (NTAFD.START %in% 0 &
                                          NTAFD.END %in% Inf))
  if (nrow(tmp.intervals) > 0) {
    tmp.intervals$NTAFD.END <- 24
    auc.intervals <- rbind(auc.intervals, tmp.intervals)
  }
  ## Compute the AUCs in a pretty order
  auc.intervals <- auc.intervals[with(auc.intervals,
                                      order(RANDID, PERD,
                                            NTAFD.START, NTAFD.END)),]
}

## Determine the rows to have in the end of the computation
if (filetype.options$minimize.pk.calculations) {
  names.to.match <- c("RANDID", "PERD", "NTAFD.START", "TRT.CDE")
} else {
  names.to.match <- c("RANDID", "PERD", "NTAFD.START", "NTAFD.END", "TRT.CDE")
}
d.ind.nca <- my.merge(unique(auc.intervals[,names.to.match]),
                      randcodes[,c("RANDID", "PERD", "TRT.CDE",
                                   "DOSE", "TRTG")])

## Find times to exclude for this subject/period/ntafd
mask.keep <- rep(TRUE, nrow(d.pk))
if (!identical(d.pk.exclusions, NA)) {
  for (i in 1:nrow(d.pk.exclusions)) {
    mask.exclude <- rep(TRUE, nrow(d.pk))
    for (n in names(d.pk.exclusions)) {
      mask.exclude <- (mask.exclude &
                       (d.pk[,n] %in% d.pk.exclusions[i,n]))
    }
    mask.keep <- mask.keep & !mask.exclude
  }
}
if (any(!mask.keep)) {
  warning(sprintf("As listed in pk.exclusion file, excluding %d rows of PK data from NCA.",
                  sum(!mask.keep)))
}

## Do the PK calculations for each subject/period combination
for (i in 1:nrow(d.ind.nca)) {
  ## What AUC intervals are requested for this subject/start time (and
  ## possibly /end time)
  mask.auc.intervals <-
    (auc.intervals$RANDID %in% d.ind.nca$RANDID[i] &
     auc.intervals$PERD %in% d.ind.nca$PERD[i] &
     auc.intervals$NTAFD.START %in% d.ind.nca$NTAFD.START[i])
  if ("NTAFD.END" %in% names.to.match) {
    mask.auc.intervals <-
      (mask.auc.intervals &
       auc.intervals$NTAFD.END %in% d.ind.nca$NTAFD.END[i])
  }
  tmp.auc.intervals <- auc.intervals[mask.auc.intervals,]
  start.time <- min(tmp.auc.intervals$NTAFD.START)
  ## make the end time Inf if "last" is given
  mask.last <- tolower(tmp.auc.intervals$NTAFD.END) %in% "last"
  if (any(mask.last)) {
    end.time <- Inf
  } else {
    end.time <- max(tmp.auc.intervals$NTAFD.END)
  }
  ## Assumes that an ID, dose, period, and day define a treatment
  df <- subset(d.pk, (RANDID %in% d.ind.nca$RANDID[i] &
                      PERD %in% d.ind.nca$PERD[i] &
                      TRT.CDE %in% d.ind.nca$TRT.CDE[i] &
                      start.time <= NTAFD &
                      NTAFD <= end.time &
                      ## PK exclusions specified by the user
                      mask.keep))
  ## Convert times to starting at the beginning of the period if
  ## applicable.
  if (filetype.options$pk.start.zero) {
    df$NTAFD <- df$NTAFD - start.time
    tmp.auc.intervals$NTAFD.START <-
      tmp.auc.intervals$NTAFD.START - start.time
    ## The as.numeric(as.character(...)) verifies that this will work
    ## with factors and characters.
    tmp.auc.intervals$NTAFD.END[!mask.last] <-
      as.numeric(as.character(tmp.auc.intervals$NTAFD.END[!mask.last])) -
        start.time
  }
  tmp.auc.intervals <- renameCol(tmp.auc.intervals,
                                 c("NTAFD.START", "NTAFD.END"),
                                 c("START", "END"))
  tmp <- nca.calc(df,
                  DEPVAR="DV",
                  TIMEVAR="NTAFD",
                  DOSE=d.ind.nca$DOSE[i],
                  AUC.intervals=tmp.auc.intervals,
                  min.hl.points=pds.pk$min.hl.points)
  ## Add all missing columns from d.ind.nca with NAs before putting
  ## the data in.
  if (!all(names(tmp) %in% names(d.ind.nca)))
    d.ind.nca[,names(tmp)[!(names(tmp) %in% names(d.ind.nca))]] <- NA
  ## Add in all calculated parameters
  d.ind.nca[i,names(tmp)] <- tmp
}

d.ind.nca$DAY <- d.ind.nca$NTAFD.START/24 + 1

## Find the rows that are not all NA.
mask.interesting <- rep(FALSE, nrow(d.ind.nca))
for (i in setdiff(names(d.ind.nca), c("RANDID", "PERD", "DOSE")))
  mask.interesting <- mask.interesting | !is.na(d.ind.nca[,i])

d.ind.nca <- d.ind.nca[mask.interesting,]
d.ind.nca <- d.ind.nca[with(d.ind.nca, order(RANDID, PERD, DOSE)),]

output.data[["csv"]][["indiv.nca.dirty"]] <- d.ind.nca
