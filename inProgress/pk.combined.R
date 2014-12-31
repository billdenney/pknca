## Functions to calculate combinations of PK parameters or to perform
## many calculations on a PK data set.

#' Compute simple (observed) PK parameters
#'
#' Extract the observed PK parameters from the \code{conc}entration
#' and \code{time} for a subject.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @return A data frame with Cmax, Tmax, Tlast, and Cmin
#' @export
pk.calc.simple <- function(conc, time, tmax.first=TRUE) {
  data.frame(## Maximum plasma concentration
             CMAX=pk.calc.cmax(conc),
             ## (first) Time of maximum plasma concentration
             TMAX=pk.calc.tmax(conc, time, use.first=tmax.first),
             ## Time of last measurable concentration
             TLAST=pk.calc.tlast(conc, time),
             ## Minimum concentration in the interval
             CMIN=pk.calc.cmin(conc))
}

## Compute all NCA parameters
##
## Inputs
##   data - data frame with raw data included
##   DEPVAR - dependent variable column name
##   TIMEVAR - independent variable column name
##   DOSE - dose to use for normalization of dose-normalized PK parameters
##   C.times - Concentration summarization times
##   AUC.intervals - interval (from TIMEVAR) for computing AUCs
##   method - AUC computation interval ("lin up/log down" (default) or
##     "linear"
##   adj.r.squared.factor - value to consider adjusted r-squared
##     values to be within tolerance.
##   min.hl.points - minimum number of points to use for half-life
##     calculations.
##   reset.times - should AUC intervals be computed starting at 0 or
##     starting at the current TIMEVAR column value?
## Output
##   A data frame containing all computed values is returned.
nca.calc <- function(data, DEPVAR="DV", TIMEVAR="NTAFD", DOSE=NA,
                     C.times=c(NA),
                     AUC.intervals=data.frame(START=c(0, 0),
                       END=c("last", Inf),
                       stringsAsFactors=FALSE),
                     method="lin up/log down",
                     adj.r.squared.factor=0.0001,
                     min.hl.points=3,
                     reset.times=TRUE) {
  data$DV <- as.numeric(data[,DEPVAR])
  data$TIME <- as.numeric(data[,TIMEVAR])
  ## Minimize the AUC intervals to just the columns required
  AUC.intervals <- AUC.intervals[,c("START", "END")]
  ## Start preparation of the output data
  OUT <- data.frame(DOSE=DOSE)
  ## Only keep positive dose values (negative indicate a dose of TBD)
  ## because they should not be used for dose-normalizing
  ## computations.
  OUT$DOSE.usable <- OUT$DOSE
  OUT$DOSE.usable[!is.na(OUT$DOSE.usable) & OUT$DOSE.usable < 0] <- NA
  ## Only compute the result if not all DV are missing.
  if(!all(is.na(data$DV))){
    data <- data[!is.na(data$DV)&!is.infinite(data$DV),]
    data <- data[order(data$TIME),]
    tmp.hl <- half.life.calc(data$DV, data$TIME,
                             min.hl.points=min.hl.points,
                             adj.r.squared.factor=adj.r.squared.factor)
    OUT <- cbind(OUT, tmp.hl)
    ## Set all outputs to missing in case they are not calculable.
    OUT$CL <- OUT$AUCPEXT <- OUT$AUMCPEXT <- OUT$MRT <- OUT$KEL <-
      OUT$THALF_EFF <- OUT$VZ <- OUT$VSS <- NA
    ## Should we try to compute AUCs?
    if (!is.na(OUT$TMAX) & (nrow(data) > 1)) {
      ## Add AUC0last if AUC0Inf is there
      inf.row <- subset(AUC.intervals,
                        Inf %in% END)
      if (nrow(inf.row) > 0) {
        inf.row$END <- "last"
        AUC.intervals <-
          unique(rbind(AUC.intervals, inf.row))
      }
      ## Compute all requested AUCs
      for (i in 1:nrow(AUC.intervals)) {
        interval <- c(AUC.intervals$START[i], AUC.intervals$END[i])
        ## Clean up the word last to lower case if applicable
        if (tolower(interval[2]) %in% "last") {
          interval[2] <- "last"
        }
        interval.text <- paste(interval, collapse="_")
        OUT$AUCTMP <-
          pk.calc.auc(data$DV, data$TIME, interval=interval,
                      OUT$LAMBDA.Z,
                      method="lin up/log down")
        OUT$AUMCTMP <-
          pk.calc.aumc(data$DV, data$TIME, interval=interval,
                       OUT$LAMBDA.Z,
                       method="lin up/log down")
        OUT <- renameCol(OUT, c("AUCTMP", "AUMCTMP"),
                         paste(c("AUC", "AUMC"), interval.text, sep=""))
      }
      if ("AUC0_Inf" %in% names(OUT)) {
        ## FIXME: This should also detect if the AUC is 0-tau for
        ## steady state.
        ## FIXME: This assumes mg and time*ng/mL for units
        OUT$CL <- pk.calc.cl(OUT$DOSE.usable, OUT$AUC0_Inf, unitconv=1000)
      }
      if (all(c("AUC0_Inf", "AUC0_last") %in% names(OUT))) {
        ## AUCInf percent extrapolated
        OUT$AUCPEXT <- pk.calc.aucpext(OUT$AUC0_last, OUT$AUC0_Inf)
        ## AUMCInf percent extrapolated
        OUT$AUMCPEXT <- pk.calc.aucpext(OUT$AUMC0_last, OUT$AUMC0_Inf)
        ## Mean residence time
        OUT$MRT <- pk.calc.mrt(OUT$AUMC0_Inf, OUT$AUC0_Inf)
        ## Elimination rate
        OUT$KEL <- pk.calc.kel(OUT$MRT)
        ## Effective half-life
        OUT$THALF_EFF <- pk.calc.thalf.eff(OUT$MRT)
        ## Terminal volume
        ## FIXME: Assumes dose in mg and concentration in ng/mL
        OUT$VZ <- pk.calc.vz(OUT$DOSE.usable, OUT$AUC0_Inf, OUT$KEL, 1000)
        ## Steady-state volume
        OUT$VSS <- pk.calc.vss(OUT$CL, OUT$MRT)
      }
      ##############################
      ## Summarize concentration at a time values.
      if (!all(is.na(C.times))) {
        for (i in C.times)
          OUT[,paste("Conc", i)] <-
            interpolate.conc(data$DV, data$TIME, C.times[i],
                             lambda.z=OUT$LAMBDA.Z,
                             interp.method=method)
      }
      ##############################
      ## Compute dose-normalized parameters
      if (!is.na(OUT$DOSE.usable)) {
        ## Dose normalize for CMAX, CMIN, CLAST, AUC, and AUMC
        for (i in grep("^(CMAX|CMIN|CLAST|AUM?C[0-9])", names(OUT),
                       value=TRUE)) {
          OUT[,paste(i, "D", sep="_")] <- OUT[,i]/OUT$DOSE.usable
        }
      }
    }
  }
  ## This column is only used for computation but not for reporting.
  ## It changed the negative doses from TBD dose descriptions to NA.
  OUT$DOSE.usable <- NULL
  OUT
}
