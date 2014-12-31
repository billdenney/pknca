## Generate a NONMEM dosing/pk dataset with some reasonable defaults.

if (!identical(d.dosing, NA) & !identical(d.pk, NA)) {
  ## Verify the dosing dataset is usable
  if (any(d.dosing$SCREENING))
    stop("Cannot generate a NONMEM dataset with dosing during screening")
  if (any(d.dosing$UNSCHEDULED))
    stop("Cannot generate a NONMEM dataset with unscheduled dosing")
  if (any(d.dosing$FOLLOWUP))
    stop("Cannot generate a NONMEM dataset with dosing during followup")
  ##############################
  ## Organize the dosing records
  tmp.dosing <- d.dosing
  if ("DONE" %in% names(tmp.dosing))
    tmp.dosing <- subset(tmp.dosing, DONE %in% "DONE")
  if ("PRIMARY.DOSE.AGENT" %in% names(tmp.dosing))
    tmp.dosing <- subset(tmp.dosing, PRIMARY.DOSE.AGENT)
  dosing.keep.cols <-
    intersect(c("STUDY", "RANDID", "PERD", "ATAFD", "NTAFD", "FORM",
                "DOSE", "TRT.CDE", "TRTG", "ROUTE", "DATETIME"),
              names(tmp.dosing))
  tmp.dosing <- tmp.dosing[,dosing.keep.cols]
  ## Rename the ROUTE and DATETIME and then merge in the dosing CMT.
  tmp.dosing <-
    merge(renameCol(tmp.dosing,
                    c("ROUTE", "DATETIME"),
                    c("CMT.text", "DATETIME.DOSE")),
          filetype.options$nonmem.dosing.cmt, all.x=TRUE)
  ## Fill in the CMT if it is missing one way or the other
  if (!("CMT" %in% names(tmp.dosing))) {
    tmp.dosing$CMT <- filetype.options$nonmem.dosing.cmt[1]
  } else {
    tmp.dosing$CMT[is.na(tmp.dosing$CMT)] <-
      filetype.options$nonmem.dosing.cmt$CMT[1]
  }
  tmp.dosing$EVID <- 1
  tmp.dosing$AMT <- tmp.dosing$DOSE

  ##############################
  ## Organize the PK records

  ## Verify the pk dataset is usable
  if (any(d.pk$SCREENING %in% TRUE))
    stop("Cannot generate a NONMEM dataset with PK during screening")
  if (any(d.pk$UNSCHEDULED %in% TRUE))
    stop("Cannot generate a NONMEM dataset with unscheduled PK")
  if (any(d.pk$FOLLOWUP %in% TRUE))
    stop("Cannot generate a NONMEM dataset with PK during followup")

  tmp.pk <- d.pk
  pk.keep.cols <-
    intersect(c("STUDY", "RANDID", "PERD", ## Subject and period vars
                "ATAFD", "NTAFD", "NTPD", "ATPD", "DATETIME", ## Time vars
                "DV.name", "DV.unit", "DV", "BLQ", "MATRIX", ## Conc vars
                "DOSE", "TRT.CDE", "TRTG"), ## General vars
              names(tmp.pk))
  tmp.pk <- tmp.pk[,pk.keep.cols]
  if ("DATETIME" %in% names(tmp.pk))
    tmp.pk <- renameCol(tmp.pk, "DATETIME", "DATETIME.PK")
  tmp.pk$EVID <- 0
  tmp.pk$CMT <- filetype.options$nonmem.pk.cmt
  ## Prep to merge them together
  pk.missing.cols <- setdiff(names(tmp.dosing), names(tmp.pk))
  tmp.pk[,pk.missing.cols] <- NA
  dosing.missing.cols <- setdiff(names(tmp.pk), names(tmp.dosing))
  tmp.dosing[,dosing.missing.cols] <- NA
  d.nonmem <- rbind(tmp.pk, tmp.dosing)
  ## Make sure that we have both ATAFD and NTAFD to simplify finding
  ## the first record.
  time.missing.cols <- setdiff(c("ATAFD", "NTAFD"), names(d.nonmem))
  d.nonmem[,time.missing.cols] <- NA
  ## Find the first record for each period and add a reset if
  ## requested.
  if (filetype.options$nonmem.reset.periods) {
    first.record <- summaryBy(ATAFD+NTAFD~STUDY+RANDID+PERD, data=d.nonmem,
                              FUN=my.min)
    first.record <- renameCol(first.record, c("ATAFD.my.min", "NTAFD.my.min"),
                              c("ATAFD", "NTAFD"))
    first.record$EVID <- 3
    first.missing.cols <- setdiff(names(d.nonmem), names(first.record))
    first.record[,first.missing.cols] <- NA
    d.nonmem <- rbind(d.nonmem, first.record)
  }
  ## Merge in the demographics.
  if (!identical(d.demog, NA) &
      !identical(filetype.options$nonmem.demog.cols, NA)) {
    tmp.demog <-
      d.demog[,unique(c("RANDID", filetype.options$nonmem.demog.cols))]
    d.nonmem <- merge(d.nonmem, tmp.demog, all.x=TRUE)
  }
  ## Merge in the vitals with interpolation as applicable.
  if (!identical(d.vitals, NA) &
      !identical(filetype.options$nonmem.vitals.interp, NA)) {
    for (n in filetype.options$nonmem.vitals.interp) {
      tmp.vitals <- subset(d.vitals, DV.name %in% n & !is.na(DV))
      if (nrow(tmp.vitals) > 0) {
        ## Is there more than one of the measurement for a subject?
        if (nrow(tmp.vitals) ==
            length(unique(tmp.vitals[,c("RANDID", "PERD")]))) {
          ## There is only a single measurement per subject or
          ## subject/period; simply merge it in.
          tmp.vitals <- renameCol(tmp.vitals[,c("RANDID", n)], "DV", n)
          d.nonmem <- my.merge(d.nonmem, tmp.vitals, all.x=TRUE)
        } else {
          ## There is more than one measurement per subject/period.
          ## Within each period, linearly interpolate the value using
          ## ATAFD (or NTAFD).
          ## FIXME: This should use ATAFD if available
          time.col <- "NTAFD"
          ## time.col <- intersect(c("ATAFD", "NTAFD"),
          ##     intersect(names(d.nonmem), names(tmp.vitals)))
          if (length(time.col) == 0) {
            warning(sprintf("Cannot interpolate %s into nonmem dataset.", n))
          } else {
            time.col <- time.col[1]
            u.id.perd <- unique(d.nonmem[,c("RANDID", "PERD")])
            d.nonmem[,n] <- NA
            for (i in 1:nrow(u.id.perd)) {
              ## Find the input rows to estimate from
              mask.id.perd.in <- 
                (tmp.vitals$RANDID %in% u.id.perd$RANDID[i] &
                 tmp.vitals$PERD %in% u.id.perd$PERD[i])
              ## Find the output rows to estimate for
              mask.id.perd.out <-
                (d.nonmem$RANDID %in% u.id.perd$RANDID[i] &
                 d.nonmem$PERD %in% u.id.perd$PERD[i] &
                 !is.na(d.nonmem[,time.col]))
              if (sum(mask.id.perd.in) == 0) {
                ## Do nothing, we can't get the value for this
                ## subject/time
              } else if (sum(mask.id.perd.in) == 1) {
                d.nonmem[mask.id.perd.out,n] <- tmp.vitals$DV[mask.id.perd.in]
              } else {
                try({
                  ## Get the estimation function
                  tmp.vitals.fun <- approxfun(tmp.vitals[mask.id.perd.in, time.col],
                                              tmp.vitals$DV[mask.id.perd.in],
                                              rule=2, ties=mean)
                  ## Estimate the values
                  d.nonmem[mask.id.perd.out,n] <-
                    tmp.vitals.fun(d.nonmem[mask.id.perd.out, time.col])
                })
              }
            }
          }
        }
      }
    }
  }
  ##############################
  ## Done creating the dataset, now make it usable
  d.nonmem$MDV <- as.numeric(d.nonmem$DV %in% c(NA, 0))
  ## Merge back in some columns so that they are carried through a
  ## period for a subject.
  u.id.perd <- unique(d.nonmem[,c("RANDID", "PERD")])
  for (n in filetype.options$nonmem.period.cols) {
    for (i in 1:nrow(u.id.perd)) {
      mask.id.perd <- (d.nonmem$RANDID %in% u.id.perd$RANDID[i] &
                       d.nonmem$PERD %in% u.id.perd$PERD[i])
      value <- na.omit(unique(d.nonmem[mask.id.perd,n]))
      if (length(value) == 1)
        d.nonmem[mask.id.perd,n] <- value
    }
  }
  ## Carry values forward in the dataset.
  for (n in filetype.options$nonmem.carry.forward.cols)
    d.nonmem[,n] <- carry.forward(d.nonmem[,n])
  ## Convert any text columns into character and then convert anything
  ## other than letters and numbers into underscores.
  for (n in names(d.nonmem)) {
    if (!(all(is.na(d.nonmem[,n]))) &
        (is.factor(d.nonmem[,n]) | is.character(d.nonmem[,n]))) {
      d.nonmem[,n] <- gsub("[^A-Za-z0-9]", "_", d.nonmem[,n])
    } else if (!(all(is.na(d.nonmem[,n]))) &
               is.logical(d.nonmem[,n])) {
      d.nonmem[,n] <- as.numeric(d.nonmem[,n])
    }
  }
  ## Convert requested covariates from text to numbers
  for (n in names(filetype.options$nonmem.character.conversion)) {
    if (n %in% names(d.nonmem)) {
      if (!(n %in% names(filetype.options$nonmem.character.conversion[[n]])))
        stop(sprintf("Invalid NONMEM character conversion for %s, missing column named %s", n, n))
      d.nonmem <- merge(d.nonmem,
                        filetype.options$nonmem.character.conversion[[n]],
                        all.x=TRUE)
    }
  }
  ## Remove placebo
  if (filetype.options$nonmem.remove.placebo) {
    d.pbo <- unique(subset(d.nonmem, DOSE %in% 0)[,c("RANDID", "TRTG")])
    d.pbo$randid.trtg <- sprintf("%s(%s)", d.pbo$RANDID, d.pbo$TRTG)
    warning(sprintf("Removing placebo subjects(treatment) from NONMEM file: %s\nTo keep placebo subjects, set filetype.options$nonmem.remove.placebo <- FALSE",
                    paste(d.pbo$randid.trtg, collapse=", ")))
    d.nonmem <- subset(d.nonmem, !(DOSE %in% 0))
  }
  ## Remove no reading subjects
  if (filetype.options$nonmem.remove.noreading) {
    d.noreading <- summaryBy(DV~STUDY+RANDID, data=d.nonmem,
                             FUN=function(x) {all(is.na(x))},
                             keep.names=TRUE)
    id.noreading <- subset(d.noreading, DV %in% TRUE)$RANDID
    if (length(id.noreading) > 0) {
      warning(sprintf("Removing subjects without any non-NA readings from NONMEM file: %s\nTo keep these subjects, set filetype.options$nonmem.remove.noreading <- FALSE",
                      paste(id.noreading, collapse=", ")))
      d.nonmem <- subset(d.nonmem, !(RANDID %in% id.noreading))
    }
  }
  ## Add subject identifier column
  d.nonmem$ID <- with(d.nonmem, as.numeric(as.factor(paste(STUDY, RANDID))))
  ## Place the columns in an order that makes sense for NONMEM
  col.order <- c("ID", "STUDY", "RANDID",
                 "PERD", "NTAFD", "ATAFD", "NTPD", "ATPD",
                 "EVID", "AMT", "DOSE", "DV", "BLQ", "MDV")
  character.cols <- sapply(d.nonmem, is.character) %in% TRUE
  factor.cols <- sapply(d.nonmem, is.factor) %in% TRUE
  numeric.cols <- !(character.cols | factor.cols)
  col.order <- c(
    col.order,
    setdiff(names(d.nonmem)[numeric.cols], col.order),
    setdiff(names(d.nonmem)[factor.cols | character.cols], col.order))
  d.nonmem <- d.nonmem[,col.order]
  ## Sort the output to a usable order.
  d.nonmem <-
    d.nonmem[with(d.nonmem,
                  order(STUDY, RANDID, PERD, NTAFD, ATAFD, -EVID, CMT)),]
  ## Convert all NAs into "."
  for (n in names(d.nonmem)) {
    d.nonmem[,n] <- as.character(d.nonmem[,n])
    d.nonmem[is.na(d.nonmem[,n]),n] <- "."
  }
  output.data[["csv"]][["nonmem.dataset"]] <- d.nonmem
} else {
  print("Missing dosing and/or pk data, cannot create NONMEM dataset.")
}
