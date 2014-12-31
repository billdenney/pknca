## Options for use within the code for setting and getting PKNCA
## default options.

.PKNCA.option.check <- list(
  adj.r.squared.factor=function(x, default=FALSE) {
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
  max.missing=function(x, default=FALSE) {
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
  auc.method=function(x, default=FALSE) {
    if (default)
      return("lin up/log down")
    match.arg(x, c("lin up/log down", "linear"))
  },
  conc.na=function(x, default=FALSE) {
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
  conc.blq=function(x, default=FALSE) {
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
  first.tmax=function(x, default=FALSE) {
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
  min.hl.points=function(x, default=FALSE) {
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
  min.span.ratio=function(x, default=FALSE) {
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
  max.aucinf.pext=function(x, default=FALSE) {
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
  min.hl.r.squared=function(x, default=FALSE) {
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
  tau.choices=function(x, default=FALSE) {
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
  single.dose.aucs=function(x, default=FALSE) {
    if (default) {
      ## It is good to put this through the specification checker in
      ## case they get out of sync during development.  (A free test
      ## case!)
      x <- data.frame(start=0,
                      end=c(24, Inf),
                      auc.type=c("AUClast", "AUCinf"),
                      half.life=c(FALSE, TRUE))
    }
    check.auc.specification(x)
  })

#' Set default options for PKNCA functions
#'
#' This function will set the default PKNCA options.  If given no
#' inputs, it will provide the current option set.  If given
#' name/value pairs, it will set the option (as in the
#' \code{\link{options}} function).  If given a name, it will return
#' the value for the parameter.  If given the \code{default} option as
#' true, it will provide the default options.
#'
#' Options are either for calculation or summary functions.
#' Calculation options are required for a calculation function to
#' report a result (otherwise the reported value will be \code{NA}).
#' Summary options are used during summarization and are used for
#' assessing what values are included in the summary.
#'
#' @param \dots options to set or get the value for
#' @param default (re)sets all default options
#' @param check check a single option given, but do not set it (for
#' validation of the values when used in another function)
#' @export
#' @return
#' If...
#' \describe{
#'   \item{no arguments are given}{returns the current options.}
#'   \item{a value is set (including the defaults)}{returns \code{NULL}}
#'   \item{a single value is requested}{the current value of that option is returned as a scalar}
#'   \item{multiple values are requested}{the current values of those options are returned as a list}
#' }
#' @examples
#'
#' PKNCA.options()
#' PKNCA.options(default=TRUE)
#' PKNCA.options("auc.method")
#' PKNCA.options(auc.method="lin up/log down", min.hl.points=3)
PKNCA.options <- function(..., default=FALSE, check=FALSE) {
  current <- get("options", envir=.PKNCAEnv)
  args <- list(...)
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

## FIXME: Which of the below is needed for a package?
.PKNCAEnv <- new.env()
assign("options",
       lapply(.PKNCA.option.check, FUN=function(x) x(default=TRUE)),
       envir=.PKNCAEnv)

.onLoad <- function(libname, pkgname) {
  PKNCA.options(default=TRUE)
}
