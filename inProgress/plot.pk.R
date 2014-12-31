############################################################
##  This Script is for Plotting the Individual PK Data
############################################################

## These lines are not necessary during normal operation, but can be
## used for troubleshooting or using this script fully independent of
## other parts of the scripts.

## Make sure that treatments sort well
d.pk$TRTG <- as.factor.smarttext(d.pk$TRTG)
## Unique subject/period combinations
u.id.perd <- unique(d.pk[,c("RANDID", "PERD", "TRTG")])

## FIXME: exclude points that were excluded somewhere in the plotting
## FIXME: plot above LOQ points differently.

## How to plot the individual data
panel.individual.pk <- function(x, y, loq=NA, col, lty, lwd, pch,
                                fill, type, hl.est=NA, use.log=FALSE,
                                ...) {
  ## BLQ will be 0 if not on the log scale and -Inf if on the log
  ## scale
  if (use.log) {
    mask.blq <- y %in% -Inf & !is.na(y)
    loq <- log(loq)
  } else {
    mask.blq <- y %in% 0 & !is.na(y)
  }
  mask.alq <- !mask.blq & !is.na(y)
  ## Plot a line through all the data
  panel.xyplot(x, y, lty=1, lwd=1, col="black", type="l", ...)
  ## If we are plotting the half-life estimation
  if (!identical(hl.est, NA)) {
    for (i in 1:nrow(hl.est)) {
      cur.hl.est <- hl.est[i,]
      cur.hl.est$TFIRST.LZ <- cur.hl.est$TFIRST.LZ + cur.hl.est$NTAFD.START
      cur.hl.est$TLAST <- cur.hl.est$TLAST + cur.hl.est$NTAFD.START
      if (!is.na(cur.hl.est$TFIRST.LZ)) {
        if (use.log) {
          hl.points <- data.frame(x=c(cur.hl.est$TFIRST.LZ, cur.hl.est$TLAST))
        } else {
          ## we will want to represent the curve instead of the line
          hl.points <-
            data.frame(x=seq(cur.hl.est$TFIRST.LZ, cur.hl.est$TLAST, length.out=30))
        }
        hl.points$y <- exp(cur.hl.est$YINT-
                           cur.hl.est$LAMBDA.Z*(hl.points$x - cur.hl.est$NTAFD.START))
        if (use.log)
          hl.points$y <- log10(hl.points$y)
        ## Plot the half-life estimate line
        if (!all(is.na(hl.points$y)))
          panel.xyplot(hl.points$x, hl.points$y, type="l", lwd=2, lty=2,
                   col="orange")
        ## Highlight the points used for the half-life (plot a slightly
        ## larger point in orange to circle them)
        mask.hl.points <- (!is.na(x) &
                           cur.hl.est$TFIRST.LZ <= x &
                           x <= cur.hl.est$TLAST)
        panel.xyplot(x[mask.hl.points], y[mask.hl.points],
                     pch=19, cex=1.3, col="orange", ...)
      }
    }
  }
  ## Plot the blq points as open circles
  panel.xyplot(x[mask.blq], y[mask.blq], pch="o", col="black",
               fill=TRUE, ...)
  ## Plot the above LOQ points as closed circles
  panel.xyplot(x[mask.alq], y[mask.alq], pch=19, col="black",
               fill=FALSE, ...)
  ## Plot the LOQ
  if (!is.na(loq)) {
    if (use.log) {
      panel.abline(h=log10(loq), lty=2, col="red", ...)
    } else {
      panel.abline(h=loq, lty=2, col="red", ...)
    }
  }
}

## Set the limit of quantification
loq <- NA
if ("BLQ" %in% names(d.pk)) {
  loq <- setdiff(d.pk$BLQ, c(NA, 0))
  if (length(loq) > 1)
    loq <- NA
}

## Determine the text to use for the concentration (like
## "Concentration (ng/ml)")
conc.unit <- unique(na.omit(d.pk$DV.unit))
if (length(conc.unit) == 1) {
  conc.text <- sprintf("Concentration (%s)", tolower(conc.unit))
} else {
  conc.text <- "Concentration"
}

############################################################
## Treatment Group PK

## Plot all individuals by treatment on the linear scale
output.data[["unblinded"]][["Summary PK"]] <- list()
output.data[["unblinded"]][["Individual PK by Treatment"]] <- list()
tmp <- subset(d.pk, !is.na(DV))
tmp$TRTG <- as.factor.smarttext(tmp$TRTG)

## Plot median by treatment, linear scale
output.data[["unblinded"]][["Summary PK"]] <-
  append(output.data[["unblinded"]][["Summary PK"]],
         list(xyplot(DV~NTAFD,
                     data=tmp,
                     groups=TRTG,
                     panel=panel.superpose,
                     panel.groups=panel.errbars,
                     point.est.fun=business.median,
                     lb.fun=NA, ub.fun=NA,
                     type=NA,
                     xlab="Time (hr)",
                     ylab=paste("Median", conc.text),
                     main=compound,
                     auto.key=list(columns=2),
                     page=watermark("Summary PK",
                       length(output.data[["unblinded"]][["Summary PK"]])+1))))

## Plot median by treatment, log scale
output.data[["unblinded"]][["Summary PK"]] <-
  append(output.data[["unblinded"]][["Summary PK"]],
         list(xyplot(DV~NTAFD,
                     data=tmp,
                     groups=TRTG,
                     panel=panel.superpose,
                     panel.groups=panel.errbars,
                     point.est.fun=business.median,
                     lb.fun=NA, ub.fun=NA,
                     type=NA,
                     xlab="Time (hr)",
                     ylab=paste("Median", conc.text),
                     main=compound,
                     auto.key=list(columns=2),
                     scales=list(relation='free',
                       alternating=FALSE,
                       y=list(log=TRUE)),
                     yscale.components=yscale.components.log10ticks,
                     page=watermark("Summary PK",
                       length(output.data[["unblinded"]][["Summary PK"]])+1))))

## Plot all individuals by treatment with PK parameters
for (i in levels(tmp$TRTG)) {
  tmp2 <- subset(tmp, TRTG %in% i)
  tmp2$RANDID <- as.factor.smarttext(tmp2$RANDID)
  tmp.pk.params <- subset(output.data[["unblinded"]][["pk"]], TRTG %in% i)
  tmp.pk.params <-
    tmp.pk.params[,setdiff(names(output.data[["unblinded"]][["pk"]]), "TRTG")]
  annotation <- paste(names(tmp.pk.params), tmp.pk.params, sep=":  ")
  ## Linear
  output.data[["unblinded"]][["Individual PK by Treatment"]] <-
    append(output.data[["unblinded"]][["Individual PK by Treatment"]],
           list(xyplot(DV~NTAFD,
                       data=tmp2,
                       type="l",
                       groups=RANDID,
                       xlab="Time (hr)",
                       ylab=conc.text,
                       main=paste(compound, i, sep="\n"),
                       auto.key=list(space="bottom", columns=4),
                       page=watermark("Individual PK by Treatment",
                         length(output.data[["unblinded"]][["Individual PK by Treatment"]])+1))))
  ## key=simpleKey(space="bottom",
  ##   text=annotation,
  ##   columns=2,
  ##   points=FALSE)))
  ## Log
  if (!all(tmp2$DV %in% c(0, NA)))
    output.data[["unblinded"]][["Individual PK by Treatment"]] <-
      append(output.data[["unblinded"]][["Individual PK by Treatment"]],
             list(xyplot(DV~NTAFD,
                         data=tmp2,
                         type="l",
                         groups=RANDID,
                         xlab="Time (hr)",
                         ylab=conc.text,
                         main=paste(compound, i, sep="\n"),
                         scales=list(y=list(log=TRUE)),
                         yscale.components=yscale.components.log10ticks,
                         auto.key=list(space="bottom", columns=4),
                         page=watermark("Individual PK by Treatment",
                           length(output.data[["unblinded"]][["Individual PK by Treatment"]])+1))))
  ## key=simpleKey(space="bottom",
  ##   text=annotation,
  ##   columns=2,
  ##   points=FALSE)))
}

## Provide the summary PK in the blinded and unblinded outputs
output.data[["blinded"]][["Summary PK"]] <-
  output.data[["unblinded"]][["Summary PK"]]

############################################################
## Individual PK

## Plot all individuals' PK on the linear scale
output.data[["unblinded"]][["Individual PK"]] <- list()
for (i in 1:nrow(u.id.perd)) {
  tmp <- subset(d.pk, (RANDID %in% u.id.perd$RANDID[i] &
                       PERD %in% u.id.perd$PERD[i]))
  ## Determine which subset of the data to plot.
  if (pk.plot.times %in% "measured") {
    tmp <- subset(tmp, !is.na(DV))
  } else if (pk.plot.times %in% "all") {
    ## Do nothing
  } else if (pk.plot.times %in% "above loq") {
    tmp <- subset(tmp, !is.na(DV) & !(DV %in% 0))
  } else {
    stop(sprintf("Invalid value for pk.plot.times: %s", pk.plot.times))
  }
  if (nrow(tmp) > 0 & !all(is.na(tmp$DV))) {
    tmp.nca <- subset(d.ind.nca, (RANDID %in% u.id.perd$RANDID[i] &
                                  PERD %in% u.id.perd$PERD[i]))
    annotation <- "No PK calculations available"
    if (nrow(tmp.nca) > 0) {
      annotation <-
        c(paste("Dose =", tmp.nca$DOSE),
          sprintf("Lambda z = %0.3g", tmp.nca$LAMBDA.Z),
          sprintf("Adjusted R-Squared = %0.3g", tmp.nca$adj.r.squared),
          sprintf("Terminal t1/2 (h) = %0.3g", tmp.nca$THALF_T),
          sprintf("Span Ratio = %0.3g", tmp.nca$SPANR))
    }
    ## FIXME: multiple dosing needs a fix here to split troughs by
    ## shingles.
    ##
    ## Setup the plot arguments for the general (non-log, non-annotated)
    ## case.
    for (j in pk.plot.types) {
      ## Plot all the different variants of individual PK requested
      use.log <- grepl("semilog", j)
      use.annotation <- grepl("annotated", j)
      ## Reload tmp so that log and linear plots can both be plotted
      ## while still setting BLQ values to NA for log.
      tmp2 <- tmp
      if (use.log)
        tmp2$DV[tmp$DV %in% 0] <- NA
      plot.args <-
        list(DV~NTAFD,
             data=tmp2,
             panel=panel.individual.pk, loq=loq,
             main=paste(compound, "\nfor ID", u.id.perd$RANDID[i],
               "with", u.id.perd$TRTG[i]),
             ylab=conc.text,
             xlab="Time (hr)",
             hl.est=tmp.nca,
             page=watermark("Individual PK",
               length(output.data[["unblinded"]][["Individual PK"]])+1))
      if (use.annotation) {
        plot.args$auto.key <- list(space="bottom",
                                   text=annotation,
                                   columns=2,
                                   points=FALSE)
      }
      if (use.log) {
        plot.args$yscale.components <- yscale.components.log10ticks
        plot.args$scales <- list(y=list(log=TRUE))
        plot.args$use.log <- TRUE
      }
      ## This will catch semi-log plots with all BLQ PK
      if (!all(is.na(tmp2$DV))) {
        output.data[["unblinded"]][["Individual PK"]] <-
          append(output.data[["unblinded"]][["Individual PK"]],
                 list(do.call(xyplot, plot.args)))
      }
    }
  }
}
