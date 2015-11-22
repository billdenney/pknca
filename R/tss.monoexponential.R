#' Compute the time to steady state using nonlinear, mixed-effects
#' modeling of trough concentrations.
#'
#' Trough concentrations are selected as concentrations at the time of
#' dosing.  An exponential curve is then fit through the data with a
#' different magnitude by treatment (as a factor) and a random
#' steady-state concentration and time to stead-state by subject (see
#' \code{random.effects} argument).
#' 
#' @param \dots See \code{\link{pk.tss.data.prep}}
#' @param tss.fraction The fraction of steady-state required for
#' calling steady-state
#' @param output Which types of outputs should be produced?
#' \code{population} is the population estimate for time to
#' steady-state (from an nlme model), \code{popind} is the individual
#' estimate (from an nlme model), \code{individual} fits each
#' individual separately with a gnls model, and \code{single} fits all
#' the data to a single gnls model.
#' @param check See \code{\link{pk.tss.data.prep}}.
#' @param verbose Describe models as they are run, show convergence of
#' the model (passed to the nlme function), and additional details
#' while running.
#' 
#' @return A scalar float for the first time when steady-state is
#' achieved or \code{NA} if it is not observed.
#' @seealso \code{\link{pk.tss}}, \code{\link{pk.tss.stepwise.linear}}
#' @references
#' Maganti L, Panebianco DL, Maes AL.  Evaluation of Methods for
#' Estimating Time to Steady State with Examples from Phase 1 Studies.
#' AAPS Journal 10(1):141-7. doi:10.1208/s12248-008-9014-y
#' @export
pk.tss.monoexponential <- function(...,
                                   tss.fraction=0.9,
                                   output=c(
                                     "population",
                                     "popind",
                                     "individual",
                                     "single"),
                                   check=TRUE,
                                   verbose=FALSE) {
  ## Check inputs
  modeldata <- pk.tss.data.prep(..., check=check)
  if (is.factor(tss.fraction) |
      !is.numeric(tss.fraction))
    stop("tss.fraction must be a number")
  if (!length(tss.fraction) == 1) {
    warning("Only first value of tss.fraction is being used")
    tss.fraction <- tss.fraction[1]
  }
  if (tss.fraction <= 0 | tss.fraction >= 1) {
    stop("tss.fraction must be between 0 and 1, exclusive")
  } else if (tss.fraction < 0.8) {
    warning("tss.fraction is usually >= 0.8")
  }
  ## Note that this will by default choose "population" if nothing is
  ## requested.
  output <- match.arg(output, several.ok=TRUE)
  if (!("subject" %in% names(modeldata))) {
    if (any(c("population", "popind", "individual") %in% output)) {
      warning("Cannot give 'population', 'popind', or 'individual' ",
              "output without multiple subjects of data")
      output <- setdiff(output, c("population", "popind", "individual"))
    }
  }
  ## Set the tss.constant so that exp(tss.constant) == tss.fraction so
  ## that the model below solves for the requested tss.
  modeldata$tss.constant <- log(1-tss.fraction)
  ret <- data.frame()
  if (any(c("population", "popind") %in% output)) {
    ret <-
      pk.tss.monoexponential.population(
        modeldata,
        output=intersect(c("population", "popind"), output),
        verbose=verbose)
  }
  if (any(c("individual", "single") %in% output)) {
    ret <-
      merge(ret,
            pk.tss.monoexponential.individual(
              modeldata,
              output=intersect(c("individual", "single"), output),
              verbose=verbose),
            all=TRUE)
  }
  ret
}

#' A helper function to generate the formula and starting values for
#' the parameters in monoexponential models.
#'
#' @param data The data used for the model
#' @return a list with elements for each of the variables
tss.monoexponential.generate.formula <- function(data) {
  ## Setup the correct ctrough.ss by treatment or not.
  if ("treatment" %in% names(data)) {
    ctrough.by <-
      list("Ctrough.ss by treatment"=list(
             formula=ctrough.ss~treatment-1,
             ## Set the starting values for the ctrough.ss as the mean
             ## concentration by treatment
             start=doBy::summaryBy(conc~treatment,
               data=data,
               FUN=mean)$conc.mean))
  } else {
    ctrough.by <-
      list("Single Ctrough.ss"=list(
             formula=ctrough.ss~1,
             ## Set the starting values for the ctrough.ss as the mean
             ## concentration.
             start=mean(data$conc)))
  }
  ## Try all combinations of Tss by treatment and single value
  tss.by <- list(
    "Single value for Tss"=list(
      formula=tss~1,
      start=stats::median(unique(data$time))))
  if ("treatment" %in% names(data))
    tss.by[["Tss by treatment"]] <- list(
      formula=tss~treatment,
      start=rep(stats::median(unique(data$time)),
        length(unique(data$treatment))))
  ## Try combinations of random effects with Ctrough.ss and
  ## Tss. (These are returned even if they aren't always used)
  ranef.by <-
    list("Ctrough.ss and Tss"=list(formula=ctrough.ss+tss~1|subject),
         "Tss"=list(formula=tss~1|subject),
         "Ctrough.ss"=list(formula=ctrough.ss~1|subject))
  ## Return the list of parameters to test
  list(ctrough.by=ctrough.by,
       tss.by=tss.by,
       ranef.by=ranef.by)
}

#' A helper function to estimate population and popind outputs for
#' monoexponential time to steady-state.
#'
#' This function is not intended to be called directly.  Please use
#' \code{pk.tss.monoexponential}.
#' 
#' If no model converges, then the
#' \code{tss.monoexponential.population} column will be set to NA. If
#' the best model does not include a random effect for subject on Tss
#' then the \code{tss.monoexponential.popind} column of the output
#' will be set to NA.
#'
#' @param data a data frame as prepared by
#' \code{\link{pk.tss.data.prep}}.  It must contain at least columns
#' for \code{subject}, \code{time}, \code{conc}, and
#' \code{tss.constant}.
#' @param output a character vector requesting the output types.
#' @param verbose Show verbose output.
#' @return A data frame with either one row (if \code{population}
#' output is provided) or one row per subject (if \code{popind} is
#' provided).  The columns will be named
#' \code{tss.monoexponential.population} and/or
#' \code{tss.monoexponential.popind}.
pk.tss.monoexponential.population <- function(data,
                                              output=c(
                                                "population",
                                                "popind"),
                                              verbose=FALSE) {
  output <- match.arg(output, several.ok=TRUE)
  test.formula <- tss.monoexponential.generate.formula(data)
  ## Loop over creating all models
  models <- list()
  for (myctrough.ss in names(test.formula$ctrough.by)) {
    for (mytss in names(test.formula$tss.by)) {
      for (myranef in names(test.formula$ranef.by)) {
        ## Describe the current model
        current.desc <-
          sprintf("Fixed effects of %s and %s; random effects of %s",
                  myctrough.ss, mytss, myranef)
        if (verbose)
          print(current.desc)
        ## Be prepared for failure
        current.model <- NA
        current.aic <- NA
        current.model.summary <- "Did not converge"
        try({
          ## Test the current model
          current.model <-
            nlme::nlme(conc~ctrough.ss*(1-exp(tss.constant*time/tss)),
                       fixed=list(
                         test.formula$ctrough.by[[myctrough.ss]]$formula,
                         test.formula$tss.by[[mytss]]$formula),
                       random=test.formula$ranef.by[[myranef]]$formula,
                       start=c(
                         test.formula$ctrough.by[[myctrough.ss]]$start,
                         test.formula$tss.by[[mytss]]$start),
                       data=data,
                       verbose=verbose)
          ## If the model converges, get the summary and AIC out.
          current.model.summary <- summary(current.model)
          current.aic <- stats::AIC(current.model)
        }, silent=!verbose)
        ## Put the current model (or model attempt) into the list of
        ## models.
        models <- append(models,
                         list(
                           list(desc=current.desc,
                                model=current.model,
                                summary=current.model.summary,
                                AIC=current.aic)))
      }
    }
  }
  ## Find the best model of the set and return the output from that one.
  all.model.summary <- AIC.list(lapply(models, function(x) x$model))
  rownames(all.model.summary) <- sapply(models, function(x) x$desc)
  if (verbose)
    print(all.model.summary)
  if (all(is.na(all.model.summary$AIC)) |
      length(all.model.summary) == 0) {
    warning("No population model for monoexponential Tss converged, no results given")
    ret <- merge(
      data.frame(tss.monoexponential.population=NA,
                 tss.monoexponential.popind=NA),
      data[,"subject",drop=FALSE],
      all=TRUE)
  } else {
    best.model <-
      models[all.model.summary$AIC %in%
             min(all.model.summary$AIC, na.rm=TRUE)][[1]]$model
    ret <- data.frame(tss.monoexponential.population=nlme::fixef(best.model)[["tss"]])
    best.ranef <- nlme::ranef(best.model)
    if ("tss" %in% names(best.ranef)) {
      ret <- merge(
        ret,
        data.frame(
          tss.monoexponential.popind=(best.ranef$tss +
                                      ret$tss.monoexponential.population),
          subject=rownames(best.ranef)),
        all=TRUE)
    } else if ("popind" %in% output) {
      warning("tss.monoexponential.popind was requested, but the best model did not include a random effect for tss.  Set to NA.")
      ret <- merge(ret,
                   data.frame(tss.monoexponential.popind=NA,
                              subject=unique(data$subject)),
                   all=TRUE)
    }
  }
  ## Return the requested columns
  ret[,intersect(c("subject", "treatment",
                   paste("tss.monoexponential", output, sep=".")),
                 names(ret))]
}

#' A helper function to estimate individual and single outputs for
#' monoexponential time to steady-state.
#'
#' This function is not intended to be called directly.  Please use
#' \code{pk.tss.monoexponential}.
#' 
#' If no model converges, then the \code{tss.monoexponential.single}
#' and/or \code{tss.monoexponential.individual} column will be set to
#' NA.
#'
#' @param data a data frame as prepared by
#' \code{\link{pk.tss.data.prep}}.  It must contain at least columns
#' for \code{subject}, \code{time}, \code{conc}, and
#' \code{tss.constant}.
#' @param output a character vector requesting the output types.
#' @param verbose Show verbose output.
#' @return A data frame with either one row (if \code{population}
#' output is provided) or one row per subject (if \code{popind} is
#' provided).  The columns will be named
#' \code{tss.monoexponential.population} and/or
#' \code{tss.monoexponential.popind}.
pk.tss.monoexponential.individual <- function(data,
                                              output=c(
                                                "individual",
                                                "single"),
                                              verbose=FALSE) {
  fit.tss <- function(d) {
    tss <- NA
    try({
      current.model <-
        nlme::gnls(conc~ctrough.ss*(1-exp(tss.constant*time/tss)),
                   params=list(
                     ctrough.ss~1,
                     tss~1),
                   start=c(
                     mean(d$conc),
                     stats::median(unique(d$time))),
                   data=d,
                   verbose=verbose)
      tss <- stats::coef(current.model)[["tss"]]
    }, silent=!verbose)
    tss
  }
  output <- match.arg(output, several.ok=TRUE)
  ## Run by treatment or a single value for the full data frame
  ret <- sapplyBy("treatment", data,
                  FUN=list(tss.monoexponential.single=fit.tss))
  if ("subject" %in% names(data) &
      "individual" %in% output) {
    ret.sub <-
      sapplyBy(c("treatment", "subject"),
               data,
               FUN=list(tss.monoexponential.individual=fit.tss))
    ret <- merge(ret, ret.sub, all=TRUE)
  }
  ## Return the requested columns
  ret[,c(intersect(names(ret),
                   c("subject", "treatment")),
         paste("tss.monoexponential", output, sep="."))]
}
