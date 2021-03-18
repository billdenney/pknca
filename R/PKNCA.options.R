## Options for use within the code for setting and getting PKNCA
## default options.

.PKNCA.option.check <- list(
  adj.r.squared.factor=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "The adjusted r^2 for the calculation of lambda.z has this factor",
        "times the number of data points added to it.  It allows for more",
        "data points to be preferred in the calculation of half-life."))
    if (default)
      return(0.0001)
    if (length(x) != 1)
      stop("adj.r.squared.factor must be a scalar")
    if (is.factor(x) |
        !is.numeric(x))
      stop("adj.r.squared.factor must be numeric (and not a factor)")
    ## Must be between 0 and 1, exclusive
    if (x <= 0 | x >= 1)
      stop("adj.r.squared.factor must be between 0 and 1, exclusive")
    if (x > 0.01)
      warning("adj.r.squared.factor is usually <0.01")
    x
  },
  max.missing=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "The maximum fraction of the data that may be missing ('NA') to",
        "calculate summary statistics with the business.* functions."))
    if (default)
      return(0.5)
    if (length(x) != 1)
      stop("max.missing must be a scalar")
    if (is.factor(x) | !is.numeric(x))
      stop("max.missing must be numeric (and not a factor)")
    ## Must be between 0 and 1, inclusive
    if (x < 0 | x >= 1)
      stop("max.missing must be between 0 and 1")
    if (x > 0.5)
      warning("max.missing is usually <= 0.5")
    x
  },
  auc.method=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "The method used to calculate the AUC and related statistics.",
        "Options are 'lin up/log down' and 'linear'."))
    if (default)
      return("lin up/log down")
    match.arg(x, c("lin up/log down", "linear"))
  },
  conc_above=function(x, default=FALSE, description=FALSE) {
    if (description) {
      return(
        "For the time_above parameter, what concentration should the value be above?"
      )
    }
    if (default)
      return(NA_real_)
    stopifnot("conc_above must be a scalar"=length(x) == 1)
    stopifnot("conc_above must be numeric"=is.numeric(x))
    x
  },
  conc.na=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "How should missing ('NA') concentration values be handled?  See",
        "help for 'clean.conc.na' for how to use this option."))
    if (default)
      return("drop")
    if (is.na(x))
      stop("conc.na must not be NA")
    if (is.factor(x)) {
      warning("conc.na may not be a factor; attempting conversion")
      x <- as.character(x)
    }
    if (tolower(x) %in% "drop") {
      x <- tolower(x)
    } else if (is.numeric(x)) {
      if (is.infinite(x)) {
        stop("When a number, conc.na must be finite")
      } else if (x < 0) {
        warning("conc.na is usually not < 0")
      }
    } else {
      stop("conc.na must either be a finite number or the text 'drop'")
    }
    x
  },
  conc.blq=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "How should below the limit of quantification (zero, 0) concentration",
        "values be handled?  See help for 'clean.conc.blq' for how to use",
        "this option."))
    if (default)
      return(list(first="keep",
                  middle="drop",
                  last="keep"))
    check.element <- function(x) {
      if (length(x) != 1)
        stop("conc.blq must be a scalar")
      if (is.na(x))
        stop("conc.blq must not be NA")
      if (is.factor(x)) {
        warning("conc.blq may not be a factor; attempting conversion")
        x <- as.character(x)
      }
      if (tolower(x) %in% c("drop", "keep")) {
        x <- tolower(x)
      } else if (is.numeric(x)) {
        if (is.infinite(x)) {
          stop("When a number, conc.blq must be finite")
        } else if (x < 0) {
          warning("conc.blq is usually not < 0")
        }
      } else {
        stop("conc.blq must either be a finite number or the text 'drop' or 'keep'")
      }
      x
    }
    if (is.list(x)) {
      extra.names <- setdiff(names(x), c("first", "last", "middle"))
      missing.names <- setdiff(c("first", "last", "middle"), names(x))
      if (length(extra.names) != 0)
        stop("When given as a list, conc.blq must only have elements named 'first', 'middle', and 'last'.")
      if (length(missing.names) != 0)
        stop("When given as a list, conc.blq must include elements named 'first', 'middle', and 'last'.")
      ## After the names are confirmed, confirm each value.
      x <- lapply(x, check.element)
    } else {
      x <- check.element(x)
    }
    x
  },
  first.tmax=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "If there is more than one concentration equal to Cmax, which time",
        "should be selected for Tmax?  If 'TRUE', the first will be selected.",
        "If 'FALSE', the last will be selected."))
    if (default)
      return(TRUE)
    if (length(x) != 1)
      stop("first.tmax must be a scalar")
    if (is.na(x))
      stop("first.tmax may not be NA")
    if (!is.logical(x)) {
      x <- as.logical(x)
      if (is.na(x)) {
        stop("Could not convert first.tmax to a logical value")
      } else {
        warning("Converting first.tmax to a logical value: ", x)
      }
    }
    x
  },
  allow.tmax.in.half.life=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "Should the concentration and time at Tmax be allowed in the",
        "half-life calculation?  'TRUE' is yes and 'FALSE' is no."))
    if (default)
      return(FALSE)
    if (length(x) != 1)
      stop("allow.tmax.in.half.life must be a scalar")
    if (is.na(x))
      stop("allow.tmax.in.half.life may not be NA")
    if (!is.logical(x)) {
      x <- as.logical(x)
      if (is.na(x)) {
        stop("Could not convert allow.tmax.in.half.life to a logical value")
      } else {
        warning("Converting allow.tmax.in.half.life to a logical value: ", ret)
      }
    }
    x
  },
  min.hl.points=function(x, default=FALSE, description=FALSE) {
    if (description)
      return("What is the minimum number of points required to calculate half-life?")
    if (default)
      return(3)
    if (length(x) != 1)
      stop("min.hl.points must be a scalar")
    if (is.factor(x))
      stop("min.hl.points cannot be a factor")
    if (!is.numeric(x))
      stop("min.hl.points must be a number")
    if (x < 2)
      stop("min.hl.points must be >=2")
    if (min(x %% 1, 1 - (x %% 1)) >
        100*.Machine$double.eps) {
      warning("Non-integer given for min.hl.points; rounding to nearest integer")
      x <- round(x)
    }
    x
  },
  min.span.ratio=function(x, default=FALSE, description=FALSE) {
    if (description)
      return("What is the minimum span ratio required to consider a half-life valid?")
    if (default)
      return(2)
    if (length(x) != 1)
      stop("min.span.ratio must be a scalar")
    if (is.factor(x))
      stop("min.span.ratio cannot be a factor")
    if (!is.numeric(x))
      stop("min.span.ratio must be a number")
    if (x <= 0)
      stop("min.span.ratio must be > 0")
    if (x < 2)
      warning("min.span.ratio is usually >= 2")
    x
  },
  max.aucinf.pext=function(x, default=FALSE, description=FALSE) {
    if (description)
      return("What is the maximum percent extrapolation to consider an AUCinf valid?")
    if (default)
      return(20)
    if (length(x) != 1)
      stop("max.aucinf.pext must be a scalar")
    if (is.factor(x))
      stop("max.aucinf.pext cannot be a factor")
    if (!is.numeric(x))
      stop("max.aucinf.pext must be a number")
    if (x <= 0)
      stop("max.aucinf.pext must be > 0")
    if (x > 25)
      warning("max.aucinf.pext is usually <=25")
    if (x < 1)
      warning("max.aucinf.pext is on the percent not ratio scale, value given is <1%")
    x
  },
  min.hl.r.squared=function(x, default=FALSE, description=FALSE) {
    if (description)
      return("What is the minimum r-squared value to consider a half-life calculation valid?")
    if (default)
      return(0.9)
    if (length(x) != 1)
      stop("min.hl.r.squared must be a scalar")
    if (is.factor(x))
      stop("min.hl.r.squared cannot be a factor")
    if (!is.numeric(x))
      stop("min.hl.r.squared must be a number")
    if (x <= 0 | x >= 1)
      stop("min.hl.r.squared must be between 0 and 1, exclusive")
    if (x < 0.9)
      warning("min.hl.r.squared is usually >= 0.9")
    x
  },
  tau.choices=function(x, default=FALSE, description=FALSE) {
    if (description)
      return(paste(
        "What values for tau (repeating interdose interval) should be",
        "considered when attempting to automatically determine the intervals",
        "for multiple dosing?  See 'choose.auc.intervals' and 'find.tau' for",
        "more information.  'NA' means automatically look at any potential",
        "interval."))
    if (default)
      return(NA)
    if (is.factor(x))
      stop("tau.choices cannot be a factor")
    if (length(x) > 1 & any(is.na(x)))
      stop("tau.choices may not include NA and be a vector")
    if (!identical(x, NA))
      if (!is.numeric(x))
        stop("tau.choices must be a number")
      if (!is.vector(x)) {
        warning("tau.choices must be a vector, converting")
        x <- as.vector(x)
      }
    x
  },
  single.dose.aucs=function(x, default=FALSE, description=FALSE) {
    if (description)
      return("When data is single-dose, what intervals should be used?")
    if (default) {
      ## It is good to put this through the specification checker in
      ## case they get out of sync during development.  (A free test
      ## case!)
      x <- data.frame(
        start=0,
        end=c(24, Inf),
        auclast=c(TRUE, FALSE),
        aucinf.obs=c(FALSE, TRUE),
        half.life=c(FALSE, TRUE),
        tmax=c(FALSE, TRUE),
        cmax=c(FALSE, TRUE))
    }
    check.interval.specification(x)
  })

#' Set default options for PKNCA functions
#' 
#' This function will set the default PKNCA options.  If given no 
#' inputs, it will provide the current option set.  If given name/value
#' pairs, it will set the option (as in the \code{\link{options}}
#' function).  If given a name, it will return the value for the
#' parameter.  If given the \code{default} option as true, it will
#' provide the default options.
#' 
#' Options are either for calculation or summary functions. Calculation
#' options are required for a calculation function to report a result
#' (otherwise the reported value will be \code{NA}). Summary options are
#' used during summarization and are used for assessing what values are
#' included in the summary.
#' 
#' See the vignette 'Options for Controlling PKNCA' for a current list
#' of options (\code{vignette("Options-for-Controlling-PKNCA", package="PKNCA")}).
#' 
#' @param \dots options to set or get the value for
#' @param default (re)sets all default options
#' @param check check a single option given, but do not set it (for 
#'   validation of the values when used in another function)
#' @param name An option name to use with the \code{value}.
#' @param value An option value (paired with the \code{name}) to set or
#'   check (if \code{NULL}, ).
#' @return If...
#' \describe{
#'   \item{no arguments are given}{returns the current options.}
#'   \item{a value is set (including the defaults)}{returns \code{NULL}}
#'   \item{a single value is requested}{the current value of that option is returned as a scalar} 
#'   \item{multiple values are requested}{the current values of those options are returned as a list}
#' }
#' @family PKNCA calculation and summary settings
#' @seealso \code{\link{PKNCA.options.describe}}
#' @examples
#' 
#' PKNCA.options()
#' PKNCA.options(default=TRUE)
#' PKNCA.options("auc.method")
#' PKNCA.options(name="auc.method")
#' PKNCA.options(auc.method="lin up/log down", min.hl.points=3)
#' @export
PKNCA.options <- function(..., default=FALSE, check=FALSE, name, value) {
  current <- get("options", envir=.PKNCAEnv)
  ## If the options have not been initialized, initialize them and
  ## then proceed.
  if (is.null(current) & !default) {
    PKNCA.options(default=TRUE)
    current <- get("options", envir=.PKNCAEnv)
  }
  args <- list(...)
  ## Put the name/value pair into the args as if they were specified
  ## like another argument.
  if (missing(name)) {
    if (!missing(value))
      stop("Cannot have a value without a name")
  } else {
    if (name %in% names(args))
      stop("Cannot give an option name both with the name argument and as a named argument.")
    if (!missing(value)) {
      args[[name]] <- value
    } else {
      args <- append(args, name)
    }
  }
  if (default & check)
    stop("Cannot request both default and check")
  if (default) {
    if (length(args) > 0)
      stop("Cannot set default and set new options at the same time.")
    ## Extract all the default values
    defaults <- lapply(.PKNCA.option.check,
                       FUN=function(x) x(default=TRUE))
    ## Set the default options
    assign("options", defaults, envir=.PKNCAEnv)
  } else if (check) {
    ## Check an option for accuracy, but don't set it
    if (length(args) != 1)
      stop("Must give exactly one option to check")
    n <- names(args)
    if (!(n %in% names(.PKNCA.option.check)))
      stop(paste("Invalid setting for PKNCA:", n))
    ## Verify the option, and return the sanitized version
    return(.PKNCA.option.check[[n]](args[[n]]))
  } else if (length(args) > 0) {
    if (is.null(names(args))) {
      ## Confirm that the settings exist
      if (length(bad.args <- setdiff(unlist(args), names(current))) > 0)
        stop(sprintf("PKNCA.options does not have value(s) for %s.",
                     paste(bad.args, collapse=", ")))
      ## Get the setting(s)
      if (length(args) == 1) {
        ret <- current[[args[[1]]]]
      } else {
        ret <- list()
        for (i in seq_len(length(args))) {
          ret[[ args[[i]] ]] <- current[[ args[[i]] ]]
        }
      }
      return(ret)
    } else {
      ## Set a value
      ## Verify values are viable and then set them.
      for (n in names(args)) {
        if (!(n %in% names(.PKNCA.option.check)))
          stop(paste("Invalid setting for PKNCA:", n))
        ## Verify and set the option value
        current[[n]] <- .PKNCA.option.check[[n]](args[[n]])
      }
      ## Assign current into the setting environment
      assign("options", current, envir=.PKNCAEnv)
    }
  } else {
    return(current)
  }
}

#' Choose either the value from an option list or the current set value for an
#' option.
#'
#' @param name The option name requested.
#' @param value A value to check for the option (\code{NULL} to choose not to
#'   check the value).
#' @param options The non-default options to choose from.
#' @return The value of the option first from the \code{options} list and if it
#'   is not there then from the current settings.
#' @family PKNCA calculation and summary settings
#' @export
PKNCA.choose.option <- function(name, value=NULL, options=list())
  if (!is.null(value)) {
    PKNCA.options(name=name, value=value, check=TRUE)
  } else if (name %in% names(options)) {
    PKNCA.options(name=name, value=options[[name]], check=TRUE)
  } else {
    PKNCA.options(name)
  }

#' Describe a PKNCA.options option by name.
#' 
#' @param name The option name requested.
#' @return A character string of the description.
#' @seealso \code{\link{PKNCA.options}}
#' @export
PKNCA.options.describe <- function(name) {
  .PKNCA.option.check[[name]](description=TRUE)
}

#' Define how NCA parameters are summarized.
#'
#' @param name The parameter name or a vector of parameter names.  It
#'   must have already been defined (see
#'   \code{\link{add.interval.col}}).
#' @param description A single-line description of the summary
#' @param point The function to calculate the point estimate for the
#'   summary.  The function will be called as \code{point(x)} and must
#'   return a scalar value (typically a number, NA, or a string).
#' @param spread Optional.  The function to calculate the spread (or
#'   variability).  The function will be called as \code{spread(x)} and
#'   must return a scalar or two-long vector (typically a number, NA, or
#'   a string).
#' @param rounding Instructions for how to round the value of point and
#'   spread.  It may either be a list or a function.  If it is a list,
#'   then it must have a single entry with a name of either "signif" or
#'   "round" and a value of the digits to round.  If a function, it is
#'   expected to return a scalar number or character string with the
#'   correct results for an input of either a scalar or a two-long
#'   vector.
#' @param reset Reset all the summary instructions
#' @return All current summary settings (invisibly)
#' @seealso \code{\link{summary.PKNCAresults}}
#' @family PKNCA calculation and summary settings
#' @examples
#' \dontrun{
#' PKNCA.set.summary(
#'   name="half.life",
#'   description="arithmetic mean and standard deviation",
#'   point=business.mean,
#'   spread=business.sd,
#'   rounding=list(signif=3)
#' )
#' }
#' @export
PKNCA.set.summary <- function(name, description, point, spread,
                              rounding=list(signif=3), reset=FALSE) {
  if (reset) {
    current <- list()
  } else {
    current <- get("summary", envir=.PKNCAEnv)
  }
  if (missing(name) & missing(point) & missing(spread)) {
    if (reset)
      assign("summary", current, envir=.PKNCAEnv)
    return(invisible(current))
  }
  # Confirm that the name exists
  if (!all(found_names <- name %in% names(get("interval.cols", envir=.PKNCAEnv)))) {
    stop(paste("You must first define the parameter name with add.interval.col.  Parameters not yet defined are:",
               paste(name[!found_names], collapse=", ")))
  }
  # Reset all names to prep for settings below
  for (current_name in name) {
    current[[current_name]] <- list()
  }
  # Confirm that description is a scalar character string
  if (!is.character(description)) {
    stop("`description` must be a character string.")
  } else if (length(description) != 1) {
    stop("`description` must be a scalar.")
  }
  for (current_name in name) {
    current[[current_name]]$description <- description
  }
  # Confirm that point is a function
  if (!is.function(point)) {
    stop("`point` must be a function")
  }
  for (current_name in name) {
    current[[current_name]]$point <- point
  }
  # Confirm that spread is a function (if given)
  if (!missing(spread)) {
    if (!is.function(spread)) {
      stop("spread must be a function")
    }
    for (current_name in name) {
      current[[current_name]]$spread <- spread
    }
  }
  # Confirm that rounding is either a single-entry list or a function
  if (is.list(rounding)) {
    if (length(rounding) != 1) {
      stop("rounding must have a single value in the list")
    }
    if (!(names(rounding) %in% c("signif", "round"))) {
      stop("When a list, rounding must have a name of either 'signif' or 'round'")
    }
    for (current_name in name) {
      current[[current_name]]$rounding <- rounding
    }
  } else if (is.function(rounding)) {
    for (current_name in name) {
      current[[current_name]]$rounding <- rounding
    }
  } else {
    stop("rounding must be either a list or a function")
  }
  # Set the summary parameters
  assign("summary", current, envir=.PKNCAEnv)
  invisible(current)
}
