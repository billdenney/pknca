#' Compute NCA parameters for each interval for each subject.
#'
#' The computations assume that all calculation options are set in
#' \code{\link{PKNCA.options}}.
#'
#' @param data A PKNCAdata object
#' @return A data frame with a row for each interval for each grouping
#' in the concentration data.
#' @seealso \code{\link{PKNCAdata}}, \code{\link{PKNCA.options}}
#' @export
pk.nca <- function(data) {
  if (nrow(data$intervals) == 0) {
    warning("No intervals given; no calculations done.")
    results <- data.frame()
  } else {
    if (identical(NA, data$dose)) {
      # If no dose information is given, add NULL dose information.
      message("No dose information provided, assuming default dosing information.")
      tmp.dose.data <- unique(getGroups(data$conc))
      data$dose <-
        PKNCAdose(data=tmp.dose.data,
                  formula=as.formula(
                    paste0(".~.|",
                           paste(names(tmp.dose.data), collapse="+"))))
    }
    if (identical(all.vars(parseFormula(data$dose)$lhs), character())) {
      ## If dose amount is not given, give a false dose column with all 
      ## NAs to simplify dose handling in subsequent steps.
      col.dose <- paste0(max(names(data$dose$data)), "X")
      data$dose$data[,col.dose] <- NA
      data$dose$formula <-
        stats::update.formula(data$dose$formula, paste0(col.dose, "~."))
    }
    ## Merge the options into the default options.
    tmp.opt <- PKNCA.options()
    tmp.opt[names(options)] <- data$options
    data$options <- tmp.opt
    splitdata <- split.PKNCAdata(data)
    ## Calculate the results
    tmp.results <-
      parallel::mclapply(X=splitdata,
                         FUN=pk.nca.intervals,
                         options=data$options)
    ## Put the group parameters with the results
    for (i in seq_len(length(tmp.results))) {
      ## If no calculations were performed, the results are NULL.
      if (!is.null(tmp.results[[i]])) {
        ## If calculations were performed, the results are non-NULL, add
        ## the grouping information to the results, if applicable.
        keep.group.names <- setdiff(names(attr(splitdata, "groupid")),
                                    names(tmp.results[[i]]))
        if (length(keep.group.names) > 0) {
          tmp.results[[i]] <-
            cbind(attr(splitdata, "groupid")[i,keep.group.names,drop=FALSE],
                  tmp.results[[i]])
        }
      }
    }
    ## Generate the outputs
    results <- dplyr::bind_rows(tmp.results)
    rownames(results) <- NULL
  }
  PKNCAresults(result=results,
               data=data,
               exclude="exclude")
}

## Subset data down to just the times of interest and then pass it
## further to the calculation routines.
##
## This is simply a helper for pk.nca
pk.nca.intervals <- function(conc.dose, intervals, options) {
  if (is.null(conc.dose$conc)) {
    ## No data; potentially placebo data (the warning would have
    ## already been generated from making the PKNCAdata object.
    return(NULL)
  }
  ret <- data.frame()
  all.intervals <- conc.dose$intervals
  pformula.conc <- parseFormula(conc.dose$conc)
  pformula.dose <- parseFormula(conc.dose$dose)
  shared.names <- all.vars(pformula.conc$groups)
  ## Column names to use
  col.conc <- all.vars(pformula.conc$lhs)
  col.time <- all.vars(pformula.conc$rhs)
  col.dose <- all.vars(pformula.dose$lhs)
  col.time.dose <- all.vars(pformula.dose$rhs)
  # Insert NA doses and dose times if they are not given
  if (!(col.dose %in% names(conc.dose$dose$data))) {
    col.dose <- paste0(max(names(conc.dose$dose$data)), "X")
    conc.dose$dose$data[[col.dose]] <- NA
  }
  if (!(col.time.dose %in% names(conc.dose$dose$data))) {
    col.time.dose <- paste0(max(names(conc.dose$dose$data)), "X")
    conc.dose$dose$data[[col.time.dose]] <- NA
  }
  for (i in seq_len(nrow(all.intervals))) {
    ## Subset the data down to the group of current interest, and make 
    ## the first column of each the dependent variable and the second 
    ## column the independent variable.
    tmpconcdata <-
      merge(conc.dose$conc$data,
            all.intervals[i, intersect(shared.names, names(all.intervals)), drop=FALSE])[,c(col.conc, col.time, conc.dose$conc$exclude)]
    tmpdosedata <-
      merge(conc.dose$dose$data,
            all.intervals[i, intersect(shared.names, names(all.intervals)), drop=FALSE])[,c(col.dose, col.time.dose, conc.dose$dose$exclude)]
    ## Choose only times between the start and end.
    mask.keep.conc <- (all.intervals$start[i] <= tmpconcdata[[col.time]] &
                         tmpconcdata[[col.time]] <= all.intervals$end[i] &
                         is.na(tmpconcdata[[conc.dose$conc$exclude]]))
    tmpconcdata <- tmpconcdata[mask.keep.conc,]
    mask.keep.dose <- (is.na(tmpdosedata[,col.time.dose]) |
                         (all.intervals$start[i] <= tmpdosedata[[col.time.dose]] &
                            tmpdosedata[[col.time.dose]] < all.intervals$end[i]) &
                         is.na(tmpdosedata[[conc.dose$dose$exclude]]))
    tmpdosedata <- tmpdosedata[mask.keep.dose,]
    ## Sort the data in time order
    tmpconcdata <- tmpconcdata[order(tmpconcdata[[col.time]]),]
    tmpdosedata <- tmpdosedata[order(tmpdosedata[[col.time.dose]]),]
    ## Setup for detailed error reporting in case it's needed
    error.preamble <-
      paste("Error with interval",
            paste(c(shared.names, c("start", "end")),
                  c(unlist(conc.dose$conc$data[1,shared.names]),
                    unlist(all.intervals[i,c("start", "end")])),
                  sep="=", collapse=", "))
    if (nrow(tmpconcdata) == 0) {
      warning(paste(error.preamble, "No data for interval", sep=": "))
    } else {
      tryCatch(
        ## Try the calculation 
        calculated.interval <-
          pk.nca.interval(conc=tmpconcdata[,col.conc],
                          time=tmpconcdata[,col.time],
                          dose=tmpdosedata[,col.dose],
                          time.dose=tmpdosedata[,col.time.dose],
                          interval=all.intervals[i,],
                          options=options),
        error=function(e) {
          e$message <- paste(error.preamble, e$message, sep=": ")
          stop(e)
        })
      ## Add all the new data into the output
      ret <- rbind(ret,
                   cbind(all.intervals[i,c("start", "end")],
                         conc.dose$conc$data[1, shared.names, drop=FALSE],
                         calculated.interval,
                         row.names=NULL))
    }
  }
  ret
}

#' Compute all PK parameters for a single concentration-time data set
#'
#' For one subject/time range, compute all available PK parameters.
#' All the internal options should be set by
#' \code{\link{PKNCA.options}} prior to running.  The only part that
#' changes with a call to this function is the \code{conc}entration
#' and \code{time}.
#'
#' @param conc Concentration measured
#' @param time Time of concentration measurement
#' @param dose Dose amount (may be a scalar or vector)
#' @param time.dose Time of the dose (must be the same length as
#'   \code{dose})
#' @param interval One row of an interval definition (see
#'   \code{\link{check.interval.specification}} for how to define the
#'   interval.
#' @param options List of changes to the default
#'   \code{\link{PKNCA.options}} for calculations.
#' @return A data frame with the start and end time along with all PK
#'   parameters for the \code{interval}
#' 
#' @seealso \code{\link{check.interval.specification}}
#' @export
pk.nca.interval <- function(conc, time,
                            dose, time.dose,
                            interval, options=list()) {
  if (!is.data.frame(interval))
    stop("interval must be a data.frame")
  if (nrow(interval) != 1)
    stop("interval must be a one-row data.frame")
  ## Prepare the return value using SDTM names
  ret <- data.frame(PPTESTCD=NA, PPORRES=NA)[-1,]
  ## Determine exactly what needs to be calculated in what order.
  ## Start with the interval specification and find any dependencies
  ## that are not listed for calculation.  Then loop over the
  ## calculations in order confirming what needs to be passed from a
  ## previous calculation to a later calculation.
  all.intervals <- get.interval.cols()
  ## Set the dose to NA if its length is zero
  if (length(dose) == 0) {
    dose <- NA
    time.dose <- NA
  }
  ## Make sure that we calculate all of the dependencies.  Do this in
  ## reverse order for dependencies of dependencies.
  for (n in rev(names(all.intervals)))
    if (interval[1,n])
      for (deps in all.intervals[[n]]$depends)
        interval[1,deps] <- TRUE
  ## Do the calculations
  for (n in names(all.intervals))
    if (interval[1,n] & !is.na(all.intervals[[n]]$FUN)) {
      call.args <- list()
      ## Prepare to call the function by setting up its arguments.
      ## Ignore the "..." argument if it exists.
      for (arg in setdiff(names(formals(get(all.intervals[[n]]$FUN))),
                          "...")) {
        if (arg == "conc") {
          call.args[[arg]] <- conc
        } else if (arg == "time") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg]] <- time - interval$start[1]
        } else if (arg == "dose") {
          call.args[[arg]] <- dose
        } else if (arg == "time.dose") {
          ## Realign the time to be relative to the start of the
          ## interval
          call.args[[arg]] <- time.dose - interval$start[1]
        } else if (arg %in% c("start", "end")) {
          ## Provide the start and end of the interval if they are requested
          call.args[[arg]] <- interval[1,arg]
        } else if (arg == "options") {
          call.args[[arg]] <- options
        } else if (any(mask.arg <- ret$PPTESTCD %in% arg)) {
          call.args[[arg]] <- ret$PPORRES[mask.arg]
        } else {
          ## Give an error if there is not a default argument.
          ## FIXME: checking if the class is a name isn't perfect.  
          if (class(formals(get(all.intervals[[n]]$FUN))[[arg]]) == "name")
            stop(sprintf(
              "Cannot find argument '%s' for NCA function '%s'",
              arg, all.intervals[[n]]$FUN))
        }
      }
      tmp.result <- do.call(all.intervals[[n]]$FUN, call.args)
      ## If the function returns a data frame, save all the returned
      ## values, otherwise, save the value returned.
      if (is.data.frame(tmp.result)) {
        ret <- rbind(ret,
                     data.frame(PPTESTCD=names(tmp.result),
                                PPORRES=unlist(tmp.result, use.names=FALSE),
                                exclude=NA_character_,
                                stringsAsFactors=FALSE))
      } else {
        ret <- rbind(ret,
                     data.frame(PPTESTCD=n,
                                PPORRES=tmp.result,
                                exclude=NA_character_,
                                stringsAsFactors=FALSE))
      }
    }
  ret
}
