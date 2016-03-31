#' Create a PKNCAconc object
#'
#' @param data A data frame with concentration, time, and the groups
#' defined in \code{formula}.
#' @param formula The formula defining the
#' \code{concentration~time|groups}
#' @param subject The column indicating the subject number (used for
#' plotting).  If not provided, this defaults to the beginning of the
#' inner groups: For example with
#' \code{concentration~time|Study+Subject/Analyte}, the inner groups
#' start with the first grouping variable before a \code{/},
#' \code{Subject}.  If there is only one grouping variable, it is
#' assumed to be the subject (e.g. \code{concentration~time|Subject}),
#' and if there are multiple grouping variables without a \code{/},
#' subject is assumed to be the last one.  For single-subject data, it
#' is assigned as \code{NULL}.
#' @param labels (optional) Labels for use when plotting.  They are a
#' named list where the names correspond to the names in the data
#' frame and the values are used for xlab and/or ylab as appropriate.
#' @param units (optional) Units for use when plotting and calculating
#' parameters.  Note that unit conversions and simplifications are not
#' done; the text is used as-is.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @export
PKNCAconc <- function(data, formula, subject, labels, units) {
  ## Verify that all the variables in the formula are columns in the
  ## data.
  if (!all(all.vars(formula) %in% names(data))) {
    stop("All of the variables in the formula must be in the data")
  }
  parsedForm <- parseFormula(formula, require.two.sided=TRUE)
  if (length(all.vars(parsedForm$lhs)) != 1)
    stop("The left hand side of the formula must have exactly one variable")
  if (length(all.vars(parsedForm$rhs)) != 1)
    stop("The right hand side of the formula (excluding groups) must have exactly one variable")
  ## Values must be unique (one value per measurement)
  key.cols <- c(all.vars(parsedForm$rhs),
                all.vars(parsedForm$groupFormula))
  if (any(mask.dup <- duplicated(data[,key.cols])))
    stop("Rows that are not unique per group and time (column names: ",
         paste(key.cols, collapse=", "),
         ") found within concentration data.  Row numbers: ",
         paste(seq_along(mask.dup)[mask.dup], collapse=", "))
  ## Assign the subject
  if (missing(subject)) {
    tmp.groups <- all.vars(parsedForm$groupFormula)
    if (length(tmp.groups) == 1) {
      subject <- tmp.groups
    } else {
      subject <- all.vars(findOperator(parsedForm$groupFormula,
                                       "/",
                                       side="left"))
      if (length(subject) == 0) {
        ## There is no / in the group formula, use the last element
        subject <- tmp.groups[length(tmp.groups)]
      } else if (length(subject) == 1) {
        ## There is a subject given; use it as is.
      } else {
        stop("Unknown how to handle subject definition from the fromula")
      }
    }
  } else {
    ## Ensure that the subject is part of the data definition and a
    ## scalar character string.
    if (!is.character(subject))
      stop("subject must be a character string")
    if (!(length(subject) == 1))
      stop("subject must be a scalar")
    if (!(subject %in% names(data)))
      stop("The subject parameter must map to a name in the data")
  }
  ret <- list(data=data,
              formula=formula,
              subject=subject)
  ## check and add labels and units
  if (!missing(labels))
    ret <- set.name.matching(ret, "labels", labels, data)
  if (!missing(units))
    ret <- set.name.matching(ret, "units", units, data)
  class(ret) <- c("PKNCAconc", class(ret))
  ret
}

## Used for setting labels and units
set.name.matching <- function(ret, name, value, data) {
  if (!missing(value)) {
    if (is.null(names(value)))
      stop(paste(name, "must be a named list"))
    if (!(all(names(labels) %in% names(data))))
      stop(paste(name, "names must match data names"))
    ret[[name]] <- value
  }
  ret
}

#' Create a PKNCAdose object
#'
#' @param data A data frame with time and the groups defined in
#'   \code{formula}.
#' @param formula The formula defining the
#'   \code{dose.amount~time|groups} where \code{time} is the time of
#'   the dosing and \code{dose.amount} is the amount administered at
#'   that time.
#' @param labels (optional) Labels for use when plotting.  They are a
#'   named list where the names correspond to the names in the data
#'   frame and the values are used for xlab and/or ylab as
#'   appropriate.
#' @param units (optional) Units for use when plotting and calculating
#'   parameters.  Note that unit conversions and simplifications are
#'   not done; the text is used as-is.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @export
PKNCAdose <- function(data, formula, labels, units) {
  ## Verify that all the variables in the formula are columns in the
  ## data.
  if (!all(all.vars(formula) %in% names(data))) {
    stop("All of the variables in the formula must be in the data")
  }
  parsedForm <- parseFormula(formula, require.two.sided=FALSE)
  if (!(length(all.vars(parsedForm$lhs)) %in% c(0, 1)))
    stop("The left hand side of the formula must have zero or one variable")
  if (length(all.vars(parsedForm$rhs)) != 1)
    stop("The right hand side of the formula (excluding groups) must have exactly one variable")
  ## Values must be unique (one value per measurement)
  key.cols <- c(all.vars(parsedForm$rhs),
                all.vars(parsedForm$groupFormula))
  if (any(mask.dup <- duplicated(data[,key.cols])))
    stop("Rows that are not unique per group and time (column names: ",
         paste(key.cols, collapse=", "),
         ") found within dosing data.  Row numbers: ",
         paste(seq_along(mask.dup)[mask.dup], collapse=", "))
  ret <- list(data=data,
              formula=formula)
  ## check and add labels and units
  if (!missing(labels))
    ret <- set.name.matching(ret, "labels", labels, data)
  if (!missing(units))
    ret <- set.name.matching(ret, "units", units, data)
  class(ret) <- c("PKNCAdose", class(ret))
  ret
}

#' Extract the formula from a PKNCAconc object.
#'
#' @param x The object to extract the formula from.
#' @param \dots Unused
#' @return A formula object
#' @export
formula.PKNCAconc <- function(x, ...)
  x$formula

#' @rdname formula.PKNCAconc
#' @export
formula.PKNCAdose <- formula.PKNCAconc

#' Get the dependent variable (left hand side of the formula) from a
#' PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @return The vector of the dependent variable from the object.
#' @export
getDepVar <- function(x, ...)
  UseMethod("getDepVar", x)

#' @export
getDepVar.PKNCAconc <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$lhs)]
}

#' @export
getDepVar.PKNCAdose <- getDepVar.PKNCAconc

#' Get the independent variable (right hand side of the formula) from
#' a PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @return The vector of the independent variable from the object.
#' @export
getIndepVar <- function(x, ...)
  UseMethod("getIndepVar", x)

#' @export
getIndepVar.PKNCAconc <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$rhs)]
}

#' @export
getIndepVar.PKNCAdose <- getIndepVar.PKNCAconc

#' Get the groups (right hand side after the \code{|} from a PKNCA
#' object.
#'
#' @param object The object to extract the data from
#' @param form The formula to extract the data from (defaults to the
#' formula from \code{object})
#' @param level optional.  If included, this specifies the level(s) of
#' the groups to include.  If a numeric scalar, include the first
#' \code{level} number of groups.  If a numeric vector, include each
#' of the groups specified by the number.  If a character vector,
#' include the named group levels.
#' @param data The data to extract the groups from (defaults to the
#' data from \code{object})
#' @param sep Unused (kept for compatibility with the nlme package)
#' @return A data frame with the (selected) group columns.
#' @export
getGroups.PKNCAconc <- function(object, form=formula(object), level,
                                data=getData(object), sep) {
  grpnames <- all.vars(parseFormula(form)$groups)
  if (!missing(level))
    if (is.factor(level) | is.character(level)) {
      level <- as.character(level)
      if (any(!(level %in% grpnames)))
        stop("Not all levels are listed in the group names.  Missing levels are: ",
             paste(setdiff(level, grpnames), collapse=", "))
      grpnames <- level
    } else if (is.numeric(level)) {
      if (length(level) == 1) {
        grpnames <- grpnames[1:level]
      } else {
        grpnames <- grpnames[level]
      }
    }
  data[, grpnames, drop=FALSE]
}

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAdose <- getGroups.PKNCAconc

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAresults <- function(object,
                                   form=formula(object$data$conc), level,
                                   data=object$result, sep) {
  ## Include the start time as a group; this may be dropped later
  grpnames <- c(all.vars(parseFormula(form)$groups), "start")
  if (!missing(level))
    if (is.factor(level) | is.character(level)) {
      level <- as.character(level)
      if (any(!(level %in% grpnames)))
        stop("Not all levels are listed in the group names.  Missing levels are: ",
             paste(setdiff(level, grpnames), collapse=", "))
      grpnames <- level
    } else if (is.numeric(level)) {
      if (length(level) == 1) {
        grpnames <- grpnames[1:level]
      } else {
        grpnames <- grpnames[level]
      }
    }
  data[, grpnames, drop=FALSE]
}

#' Print and/or summarize a PKNCAconc or PKNCAdose object.
#'
#' @param x The object to print
#' @param object The object to summarize
#' @param n The number of rows of data to show (see
#' \code{\link{head}})
#' @param summarize Summarize the nested number of groups
#' @param \dots Arguments passed to \code{print.formula} and
#' \code{print.data.frame}
#' @export
print.PKNCAconc <- function(x, n=6, summarize=FALSE, ...) {
  cat(sprintf("Formula for concentration:\n "))
  print(stats::formula(x), ...)
  if (is.na(x$subject)) {
    cat("As a single-subject dataset.\n")
  } else {
    cat(sprintf("With %d subjects defined in the '%s' column.\n",
                length(unique(x$data[,x$subject])),
                x$subject))
  }
  if (summarize) {
    cat("\n")
    grp <- getGroups(x)
    if (ncol(grp) > 0) {
      tmp.summary <- data.frame(Group.Name=names(grp),
                                Count=0)
      for (i in 1:ncol(grp))
        tmp.summary$Count[i] <- nrow(unique(grp[,1:i,drop=FALSE]))
      cat("Group summary:\n")
      names(tmp.summary) <- gsub("\\.", " ", names(tmp.summary))
      print.data.frame(tmp.summary, row.names=FALSE)
    } else {
      cat("No groups.\n")
    }
  }
  if (n != 0) {
    if (n >= nrow(x$data)) {
      cat("\nData for concentration:\n")
    } else if (n < 0) {
      cat(sprintf("\nFirst %d rows of concentration data:\n",
                  nrow(x$data)+n))
    } else {
      cat(sprintf("\nFirst %d rows of concentration data:\n",
                  n))
    }
    print.data.frame(utils::head(x$data, n=n), ..., row.names=FALSE)
  }
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAconc <- function(object, n=0, summarize=TRUE, ...)
  print.PKNCAconc(object, n=n, summarize=summarize)

#' @rdname print.PKNCAconc
#' @export
print.PKNCAdose <- function(x, n=6, summarize=FALSE, ...) {
  cat("Formula for dosing:\n ")
  print(stats::formula(x), ...)
  if (summarize) {
    cat("\n")
    grp <- getGroups(x)
    if (ncol(grp) > 0) {
      tmp.summary <- data.frame(Group.Name=names(grp),
                                Count=0)
      for (i in 1:ncol(grp))
        tmp.summary$Count[i] <- nrow(unique(grp[,1:i,drop=FALSE]))
      cat("Group summary:\n")
      names(tmp.summary) <- gsub("\\.", " ", names(tmp.summary))
      print.data.frame(tmp.summary, row.names=FALSE)
    } else {
      cat("No groups.\n")
    }
  }
  if (n != 0) {
    if (n >= nrow(x$data)) {
      cat("\nData for dosing:\n")
    } else if (n < 0) {
      cat(sprintf("\nFirst %d rows of dosing data:\n",
                  nrow(x$data)+n))
    } else {
      cat(sprintf("\nFirst %d rows of dosing data:\n",
                  n))
    }
    print.data.frame(utils::head(x$data, n=n), ..., row.names=FALSE)
  }
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAdose <- summary.PKNCAconc

#' Extract all the original data from a PKNCAconc or PKNCAdose object
#' @param object R object to extract the data from.
#' @export
getData.PKNCAconc <- function(object)
  object$data

#' @rdname getData.PKNCAconc
#' @export
getData.PKNCAdose <- getData.PKNCAconc

#' Extract the columns used in the formula (in order) from a PKNCAconc
#' or PKNCAdose object.
#'
#' @param formula The object to use (parameter name is \code{formula}
#' to use the generic function)
#' @param \dots Unused
#' @return A data frame with the columns from the object in formula
#' order.
#' @export
model.frame.PKNCAconc <- function(formula, ...)
  formula$data[, all.vars(formula$formula), drop=FALSE]

#' @rdname model.frame.PKNCAconc
#' @export
model.frame.PKNCAdose <- model.frame.PKNCAconc

#' Create a PKNCAdata object.
#'
#' \code{PKNCAdata} combines \code{PKNCAconc} and \code{PKNCAdose} and
#' adds in the intervals for PK calculations.
#' 
#' @param data.conc Concentration data as a \code{PKNCAconc} object or
#'   a data frame
#' @param data.dose Dosing data as a \code{PKNCAdose} object
#' @param formula.conc Formula for making a \code{PKNCAconc} object
#'   with \code{data.conc}.  This must be given if \code{data.conc} is
#'   a data.frame, and it must not be given if \code{data.conc} is a
#'   \code{PKNCAconc} object.
#' @param formula.dose Formula for making a \code{PKNCAdose} object
#'   with \code{data.dose}.  This must be given if \code{data.dose} is
#'   a data.frame, and it must not be given if \code{data.dose} is a
#'   \code{PKNCAdose} object.
#' @param intervals A data frame with the AUC interval specifications
#'   as defined in \code{\link{check.interval.specification}}.  If
#'   missing, this will be automatically chosen by
#'   \code{\link{choose.auc.intervals}}.
#' @param options List of changes to the default
#'   \code{\link{PKNCA.options}} for calculations.
#' @param ... arguments passed to \code{PKNCAdata.default}
#' @return A PKNCAdata object with concentration, dose, interval, and
#'   calculation options stored (note that PKNCAdata objects can also
#'   have results after a NCA calculations are done to the data).
#' @seealso \code{\link{PKNCAconc}}, \code{\link{PKNCAdose}},
#'   \code{\link{choose.auc.intervals}}
#' @export
PKNCAdata <- function(data.conc, data.dose, ...)
  UseMethod("PKNCAdata", data.conc)

## Ensure that arguments are reversible
#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAconc <- function(data.conc, data.dose, ...)
  PKNCAdata.default(data.conc=data.conc, data.dose=data.dose, ...)

#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAdose <- function(data.conc, data.dose, ...)
  ## Swap the arguments
  PKNCAdata.default(data.dose=data.conc, data.conc=data.dose, ...)

#' @rdname PKNCAdata
#' @export
PKNCAdata.default <- function(data.conc, data.dose, ...,
                              formula.conc, formula.dose,
                              intervals, options=list()) {
  ret <- list()
  ## Generate the conc element
  if (inherits(data.conc, "PKNCAconc")) {
    if (!missing(formula.conc))
      warning("data.conc was given as a PKNCAconc object.  Ignoring formula.conc")
    ret$conc <- data.conc
  } else {
    ret$conc <- PKNCAconc(data.conc, formula.conc)
  }
  ## Generate the data element
  if (inherits(data.dose, "PKNCAdose")) {
    if (!missing(formula.dose))
      warning("data.dose was given as a PKNCAdose object.  Ignoring formula.dose")
    ret$dose <- data.dose
  } else {
    ret$dose <- PKNCAdose(data.dose, formula.dose)
  }
  ## Check the options
  if (!is.list(options))
    stop("options must be a list.")
  if (length(options) > 0) {
    if (is.null(names(options)))
      stop("options must have names.")
    for (n in names(options)) {
      tmp.opt <- list(options[[n]], TRUE)
      names(tmp.opt) <- c(n, "check")
      do.call(PKNCA.options, tmp.opt)
    }
  }
  ret$options <- options
  ## Check the AUC intervals
  if (missing(intervals)) {
    ## Generate the intervals for each grouping of concentration and
    ## dosing.
    tmp.conc.dose <-
      merge(conc=doBy::splitBy(parseFormula(ret$conc)$groupFormula,
              ret$conc$data),
            dose=doBy::splitBy(parseFormula(ret$dose)$groupFormula,
              ret$dose$data))
    groupid <- attributes(tmp.conc.dose)$groupid
    rownames(groupid) <- NULL
    intervals <- data.frame()
    indep.var.conc <- all.vars(parseFormula(ret$conc)$rhs)
    indep.var.dose <- all.vars(parseFormula(ret$dose)$rhs)
    for (i in seq_len(nrow(groupid))) {
      tmp.group <- groupid[i,,drop=FALSE]
      if (!is.null(tmp.conc.dose[[i]]$conc)) {
        rownames(tmp.group) <- NULL
        new.intervals <-
          cbind(
            tmp.group,
            choose.auc.intervals(tmp.conc.dose[[i]]$conc[,indep.var.conc],
                                 tmp.conc.dose[[i]]$dose[,indep.var.dose],
                                 single.dose.aucs=PKNCA.choose.option("single.dose.aucs",
                                   options)))
        intervals <-
          rbind(intervals, new.intervals)
      } else {
        warning("No intervals generated due to no concentration data for ",
                paste(names(tmp.group),
                      unlist(lapply(tmp.group, as.character)),
                      sep="=", collapse=", "))
      }
    }
  }
  ret$intervals <- check.interval.specification(intervals)
  ## Assign the class and give it all back to the user.
  class(ret) <- c("PKNCAdata", class(ret))
  ret
}

#' Print a PKNCAdata object
#' @param x The object to print
#' @param ... Arguments passed on to \code{\link{print.PKNCAconc}} and
#' \code{\link{print.PKNCAdose}}
#' @export
print.PKNCAdata <- function(x, ...) {
  print.PKNCAconc(x$conc, ...)
  print.PKNCAdose(x$dose, ...)
  cat(sprintf("\nWith %d rows of AUC specifications.\n",
              nrow(x$intervals)))
  if (length(x$options) == 0) {
    cat("No options are set differently than default.\n")
  } else {
    cat("Options changed from default are:\n")
    print(x$options)
  }
}

#' Summarize a PKNCAdata object showing important details about the
#' concentration, dosing, and interval information.
#' @param object The PKNCAdata object to summarize.
#' @param ... arguments passed on to \code{\link{print.PKNCAdata}}
#' @export
summary.PKNCAdata <- function(object, ...)
  print.PKNCAdata(object, summarize=TRUE, ...)

#' Generate a PKNCAresults object
#'
#' This function should not be run directly.  The object is created
#' for summarization and plotting.
#'
#' @param result a data frame with NCA calculation results and groups.
#' Each row is one interval and each column is a group name or the
#' name of an NCA parameter.
#' @param data The PKNCAdata used to generate the result
#' @param provenance Data and calculation provenance
#' @return A PKNCAresults object with each of the above within.
#' @export
PKNCAresults <- function(result, data, provenance) {
  ## Add all the parts into the object
  ret <- list(result=result,
              data=data,
              provenance=provenance)
  class(ret) <- c("PKNCAresults", class(ret))
  ret
}

#' Extract the parameter results from a PKNCAresults and return them
#' as a data frame.
#'
#' @param x The object to extract results from
#' @param ... Ignored (for compatibility with generic
#' \code{\link{as.data.frame}}
#' @return A data frame of results
#' @export
as.data.frame.PKNCAresults <- function(x, ...) {
  x$result
}

#' Count the number of values that are not NA
#'
#' @param x The object to count non-NA values within.
#' @return A scalar count of the non-NA values.
#' @export
count.non.missing <- function(x)
  sum(!is.na(x))

#' During the summarization of PKNCAresults, do the rounding of values
#' based on the instructions given.
#'
#' @param x The values to summarize
#' @param name The NCA parameter name (matching a parameter name in
#' \code{\link{PKNCA.set.summary}})
#' @return A string of the rounded value
#' @export
roundingSummarize <- function(x, name) {
  summaryInstructions <- PKNCA.set.summary()
  if (!(name %in% names(summaryInstructions)))
    stop(name, "is not in the summarization instructions from PKNCA.set.summary")
  roundingInstructions <- summaryInstructions[[name]]$rounding
  if (is.function(roundingInstructions)) {
    ret <- roundingInstructions(x)
  } else if (is.list(roundingInstructions)) {
    if (length(roundingInstructions) != 1)
      stop("Cannot interpret rounding instructions for ", name)
    if ("signif" == names(roundingInstructions)) {
      ret <- signifString(x, roundingInstructions$signif)
    } else if ("round" == names(roundingInstructions)) {
      ret <- roundString(x, roundingInstructions$round)
    } else {
      stop("Invalid rounding instruction list name for ", name)
    }
  }
  if (!is.character(ret))
    ret <- as.character(ret)
  ret
}

#' Summarize PKNCA results
#'
#' @param object The results to summarize
#' @param drop.group Which group(s) should be dropped from the
#'   formula?
#' @param not.requested.string A character string to use when a
#'   parameter summary was not requested for a parameter within an
#'   interval.
#' @param not.calculated.string A character string to use when a
#'   parameter summary was requested, but the point estimate AND
#'   spread calculations (if applicable) returned \code{NA}.
#' @param ... Ignored.
#' @return A data frame of NCA parameter results summarized according
#'   to the summarization settings.
#' @seealso \code{\link{PKNCA.set.summary}}
#' @export
summary.PKNCAresults <- function(object, ...,
                                 drop.group=object$data$conc$subject,
                                 not.requested.string=".",
                                 not.calculated.string="NC") {
  allGroups <- getGroups(object)
  groups <- unique(c("start", "end",
                     setdiff(names(allGroups), drop.group)))
  summaryFormula <- stats::as.formula(paste0("~", paste(groups, collapse="+")))
  summaryInstructions <- PKNCA.set.summary()
  ## Find any parameters that request any summaries
  resultDataCols <- 
    lapply(object$data$intervals[,setdiff(names(object$data$intervals),
                                          c(groups, drop.group,
                                            "start", "end")),
             drop=FALSE],
           FUN=any)
  resultDataCols <- as.data.frame(resultDataCols[unlist(resultDataCols)])
  ret <- cbind(unique(object$result[, groups, drop=FALSE]),
               resultDataCols)
  ret[,setdiff(names(ret), groups)] <- not.requested.string
  ## Loop over every group that needs summarization
  for (i in seq_len(nrow(ret)))
    ## Loop over every column that needs summarziation
    for (n in setdiff(names(ret), groups)) {
      ## Select the rows of the intervals that match the current row
      ## from the return value.
      current.interval <-
        merge(ret[i, groups, drop=FALSE],
              object$data$intervals[,intersect(names(object$data$intervals),
                                               c(groups, n))])
      if (any(current.interval[,n])) {
        currentData <- merge(
          ret[i, groups, drop=FALSE],
          object$result[object$result$PPTESTCD %in% n,,drop=FALSE])
        if (nrow(currentData) == 0) {
          warning("No results to summarize for ", n, " in result row ", i)
        } else {
          ## Calculation is required
          point <- summaryInstructions[[n]]$point(
            currentData$PPORRES)
          na.point <- is.na(point)
          na.spread <- NA
          ## Round the point estimate
          point <- roundingSummarize(point, n)
          current <- point
          if ("spread" %in% names(summaryInstructions[[n]])) {
            spread <- summaryInstructions[[n]]$spread(
              currentData$PPORRES)
            na.spread <- all(is.na(spread))
            if (na.spread) {
              ## The spread couldn't be calculated, so show that
              spread <- not.calculated.string
            } else {
              ## Round the spread
              spread <- roundingSummarize(spread, n)
            }
            ## Collapse the spread into a usable form if it is
            ## longer than one (e.g. a range or a confidence
            ## interval) and put brackets around it.
            spread <- paste0(" [", paste(spread, collapse=", "), "]")
            current <- paste0(current, spread)
          }
          ## Determine if the results were all missing, and if so, give
          ## the not.calculated.string
          if (na.point & (na.spread %in% c(NA, TRUE))) {
            ret[i,n] <- not.calculated.string
          } else {
            ret[i,n] <- current
          }
        }
      }
    }
  ret
}
