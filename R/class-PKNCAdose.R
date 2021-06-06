#' Create a PKNCAdose object
#' 
#' @param data A data frame with time and the groups defined in 
#'   \code{formula}.
#' @param formula The formula defining the 
#'   \code{dose.amount~time|groups} where \code{time} is the time of the
#'   dosing and \code{dose.amount} is the amount administered at that 
#'   time (see Details).
#' @param route Define the route of administration.  The value may be 
#'   either a column name from the \code{data} (checked first) or a 
#'   character string of either \code{"extravascular"} or 
#'   \code{"intravascular"} (checked second).  If given as a column 
#'   name, then every value of the column must be either 
#'   \code{"extravascular"} or \code{"intravascular"}.
#' @param rate,duration (optional) for \code{"intravascular"} dosing, 
#'   the rate or duration of dosing.  If given as a character string, it
#'   is the name of a column from the \code{data}, and if given as a 
#'   number, it is the value for all doses.  Only one may be given, and 
#'   if neither is given, then the dose is assumed to be a bolus 
#'   (\code{duration=0}).  If \code{rate} is given, then the dose amount
#'   must be given (the left hand side of the \code{formula}).
#' @param time.nominal (optional) The name of the nominal time column
#'   (if the main time variable is actual time.  The \code{time.nominal}
#'   is not used during calculations; it is available to assist with
#'   data summary and checking.
#' @param exclude (optional) The name of a column with concentrations to
#'   exclude from calculations and summarization.  If given, the column 
#'   should have values of \code{NA} or \code{""} for concentrations to
#'   include and non-empty text for concentrations to exclude.
#' @param ... Ignored.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @details The \code{formula} for a \code{PKNCAdose} object can be 
#'   given three ways: one-sided (missing left side), one-sided (missing
#'   right side), or two-sided.  Each of the three ways can be given 
#'   with or without groups.  When given one-sided missing the left 
#'   side, the left side can either be omitted or can be given as a 
#'   period (\code{.}): \code{~time|treatment+subject} and 
#'   \code{.~time|treatment+subject} are identical, and dose-related NCA
#'   parameters will all be reported as not calculable (for example, 
#'   clearance).  When given one-sided missing the right side, the right
#'   side must be specified as a period (\code{.}): 
#'   \code{dose~.|treatment+subject}, and only a single row may be given
#'   per group.  When the right side is missing, PKNCA assumes that the 
#'   same dose is given in every interval.  When given as a two-sided 
#'   formula
#' @family PKNCA objects
#' @export
PKNCAdose <- function(data, ...)
  UseMethod("PKNCAdose")

#' @rdname PKNCAdose
#' @export
PKNCAdose.default <- function(data, ...)
  PKNCAdose.data.frame(as.data.frame(data), ...)
#' @rdname PKNCAdose
#' @export
PKNCAdose.tbl_df <- function(data, ...)
  PKNCAdose.data.frame(as.data.frame(data), ...)

#' @rdname PKNCAdose
#' @export
PKNCAdose.data.frame <- function(data, formula, route, rate, duration,
                                 time.nominal, exclude, ...) {
  ## The data must have... data
  if (nrow(data) == 0) {
    stop("data must have at least one row.")
  }
  ## Check inputs
  if (!missing(time.nominal)) {
    if (!(time.nominal %in% names(data))) {
      stop("time.nominal, if given, must be a column name in the input data.")
    }
  }
  ## Verify that all the variables in the formula are columns in the
  ## data.
  parsedForm <- parseFormula(formula, require.two.sided=FALSE)
  ## Check for variable existence and length
  if (!(length(all.vars(parsedForm$lhs)) %in% c(0, 1)))
    stop("The left side of the formula must have zero or one variable")
  if (!(identical(parsedForm$lhs, NA) ||
        all.vars(parsedForm$lhs) %in% c(".", names(data)))) {
    stop("The left side formula must be a variable in the data, empty, or '.'.")
  }
  if (length(all.vars(parsedForm$rhs)) != 1)
    stop("The right side of the formula (excluding groups) must have exactly one variable")
  if (!(all.vars(parsedForm$rhs) %in% c(".", names(data)))) {
    stop("The right side formula must be a variable in the data or '.'.")
  }
  if (!all(all.vars(parsedForm$groups) %in% names(data))) {
    stop("All of the variables in the groups must be in the data")
  }
  ## Values must be unique (one value per measurement)
  key.cols <- c(setdiff(all.vars(parsedForm$rhs), "."),
                all.vars(parsedForm$groupFormula))
  if (any(mask.dup <- duplicated(data[,key.cols])))
    stop("Rows that are not unique per group and time (column names: ",
         paste(key.cols, collapse=", "),
         ") found within dosing data.  Row numbers: ",
         paste(seq_along(mask.dup)[mask.dup], collapse=", "))
  ret <- list(data=data,
              formula=formula)
  class(ret) <- c("PKNCAdose", class(ret))
  if (missing(exclude)) {
    ret <- setExcludeColumn(ret)
  } else {
    ret <- setExcludeColumn(ret, exclude=exclude)
  }
  mask.indep <- is.na(getIndepVar.PKNCAdose(ret))
  if (any(mask.indep) & !all(mask.indep)) {
    stop("Some but not all values are missing for the independent variable, please see the help for PKNCAdose for how to specify the formula and confirm that your data has dose times for all doses.")
  }
  if (missing(route)) {
    ret <- setRoute.PKNCAdose(ret)
  } else {
    ret <- setRoute.PKNCAdose(ret, route)
  }
  ret <- setDuration.PKNCAdose(ret, duration=duration,
                               rate=rate, dose=getDepVar.PKNCAdose(ret))
  if (!missing(time.nominal)) {
    ret <-
      setAttributeColumn(object=ret,
                         attr_name="time.nominal",
                         col_name=time.nominal)
  }
  ret
}

#' Set the dosing route
#' 
#' @param object A PKNCAdose object
#' @param route A character string indicating one of the following:  the
#'   column from the data which indicates the route of administration, a
#'   scalar indicating the route of administration for all subjects, or
#'   a vector indicating the route of administration for each dose in
#'   the dataset.
#' @param ... Arguments passed to another setRoute function
#' @return The object with an updated route
#' @export
setRoute <- function(object, ...) {
  UseMethod("setRoute")
}
#' @rdname setRoute
#' @export
setRoute.PKNCAdose <- function(object, route, ...) {
  if (missing(route)) {
    object <-
      setAttributeColumn(object=object,
                         attr_name="route",
                         default_value="extravascular",
                         message_if_default="Assuming route of administration is extravascular")
  } else {
    object <-
      setAttributeColumn(object=object,
                         attr_name="route",
                         col_or_value=route)
  }
  if (!all(tolower(getAttributeColumn(object=object, attr_name="route")[[1]]) %in%
           c("extravascular", "intravascular"))) {
    stop("route must have values of either 'extravascular' or 'intravascular'.  Please set to one of those values and retry.")
  }
  object
}

#' Set the duration of dosing or measurement
#' 
#' @param object An object to set a duration on
#' @param duration The value to set for the duration or the name of the
#'   column in the data to use for the duration.
#' @param rate (for PKNCAdose objects only) The rate of infusion
#' @param dose (for PKNCAdose objects only) The dose amount
#' @param ... Arguments passed to another setDuration function
#' @return The object with duration set
#' @export
setDuration <- function(object, ...)
  UseMethod("setDuration")
#' @rdname setDuration
#' @export
setDuration.PKNCAdose <- function(object, duration, rate, dose, ...) {
  if (missing(duration) & missing(rate)) {
    object <- setAttributeColumn(object=object, attr_name="duration", default_value=0,
                                 message_if_default="Assuming instant dosing (duration=0)")
                                
  } else if (!missing(duration) & !missing(rate)) {
    stop("Both duration and rate cannot be given at the same time")
    # TODO: A consistency check could be done, but that would get into
    # requiring near-equal checks for floating point error.
  } else if (!missing(duration)) {
    object <- setAttributeColumn(object=object, attr_name="duration", col_or_value=duration)
  } else if (!missing(rate) & !missing(dose)) {
    tmprate <- getColumnValueOrNot(object$data, rate, "rate")
    tmpdose <- getColumnValueOrNot(object$data, dose, "dose")
    duration <- tmpdose$data[[tmpdose$name]]/tmprate$data[[tmprate$name]]
    object <- setAttributeColumn(object=object, attr_name="duration", col_or_value=duration)
  }
  duration.val <- getAttributeColumn(object=object, attr_name="duration")[[1]]
  if (is.numeric(duration.val) &&
      !any(is.na(duration.val)) &&
      !any(is.infinite(duration.val)) &&
      all(duration.val >= 0)) {
    # It passes
  } else {
    stop("duration must be numeric without missing (NA) or infinite values, and all values must be >= 0")
  }
  object
}

#' @rdname formula.PKNCAconc
#' @export
#' @importFrom stats formula
formula.PKNCAdose <-  function(x, ...) {
  x$formula
}

#' @rdname model.frame.PKNCAconc
#' @export
#' @importFrom stats model.frame
model.frame.PKNCAdose <- function(formula, ...) {
  cbind(getDepVar.PKNCAdose(formula),
        getIndepVar.PKNCAdose(formula),
        getGroups.PKNCAdose(formula))
}

#' @export
getDepVar.PKNCAdose <- function(x, ...) {
  parsedForm <- parseFormula(x$formula, require.two.sided=FALSE)
  if (identical(parsedForm$lhs, NA) ||
      identical(all.vars(parsedForm$lhs), ".")) {
    rep(NA_integer_, nrow(x$data))
  } else {
    x$data[, all.vars(parseFormula(x)$lhs)]
  }
}

#' @export
getIndepVar.PKNCAdose <- function(x, ...) {
  parsedForm <- parseFormula(x$formula, require.two.sided=FALSE)
  if (identical(parsedForm$rhs, NA) ||
      identical(all.vars(parsedForm$rhs), ".")) {
    rep(NA_integer_, nrow(x$data))
  } else {
    x$data[, all.vars(parseFormula(x)$rhs)]
  }
}

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAdose <- function(...) {
  getGroups.PKNCAconc(...)
}

#' @rdname getData.PKNCAconc
#' @export
#' @importFrom nlme getData
getData.PKNCAdose <-  function(object)
  object$data

#' @rdname getDataName
getDataName.PKNCAdose <- function(object)
  "data"

#' @rdname print.PKNCAconc
#' @export
#' @importFrom utils head
print.PKNCAdose <- function(x, n=6, summarize=FALSE, ...) {
  cat("Formula for dosing:\n ")
  print(stats::formula(x), showEnv=FALSE, ...)
  if (!is.null(time_nom_data <- getAttributeColumn(x, attr_name="time.nominal", warn_missing=c()))) {
    cat("Nominal time column is: ", names(time_nom_data), "\n", sep="")
  } else {
    cat("Nominal time column is not specified.\n")
  }
  if (summarize) {
    cat("\n")
    grp <- getGroups.PKNCAdose(x)
    if (ncol(grp) > 0) {
      tmp.summary <- as.data.frame(
        lapply(grp, FUN=function(y) length(unique(y))))
      cat("Number unique entries in each group:\n")
      print.data.frame(tmp.summary, row.names=FALSE)
    } else {
      cat("No groups.\n")
    }
  } else {
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
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAdose <- summary.PKNCAconc

#' @rdname split.PKNCAconc
#' @export
split.PKNCAdose <- split.PKNCAconc
