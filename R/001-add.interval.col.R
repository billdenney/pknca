# Setup the default options
.PKNCAEnv <- new.env(parent=emptyenv())
assign("options", NULL, envir=.PKNCAEnv)
assign("summary", list(), envir=.PKNCAEnv)
assign("interval.cols", list(), envir=.PKNCAEnv)

#' Add columns for calculations within PKNCA intervals
#'
#' @param name The column name as a character string
#' @param FUN The function to run (as a character string) or `NA` if the
#'   parameter is automatically calculated when calculating another parameter.
#' @param values Valid values for the column
#' @param depends Character vector of columns that must be run before this
#'   column.
#' @param desc A human-readable description of the parameter (<=40 characters to
#'   comply with SDTM)
#' @param sparse Is the calculation for sparse PK?
#' @param unit_type The type of units to use for assigning and converting units.
#' @param pretty_name The name of the parameter to use for printing in summary
#'   tables with units.  (If an analysis does not include units, then the normal
#'   name is used.)
#' @param formalsmap A named list mapping parameter names in the function call
#'   to NCA parameter names.  See the details for information on use of
#'   `formalsmap`.
#' @param datatype The type of data used for the calculation
#' @returns NULL (Calling this function has a side effect of changing the
#'   available intervals for calculations)
#'
#' @details The `formalsmap` argument enables mapping some alternate formal
#' argument names to parameters.  It is used to generalize functions that may
#' use multiple similar arguments (such as the variants of mean residence time).
#' The names of the list should correspond to function formal parameter names
#' and the values should be one of the following:
#'
#' \itemize{
#'   \item{For the current interval:}
#'   \describe{
#'     \item{character strings of NCA parameter name}{The value of the parameter calculated for the current interval.}
#'     \item{"conc"}{Concentration measurements for the current interval.}
#'     \item{"time"}{Times associated with concentration measurements for the current interval (values start at 0 at the beginning of the current interval).}
#'     \item{"volume"}{Volume associated with concentration measurements for the current interval (typically applies for excretion parameters like urine).}
#'     \item{"duration.conc"}{Durations associated with concentration measurements for the current interval.}
#'     \item{"dose"}{Dose amounts assocuated with the current interval.}
#'     \item{"time.dose"}{Time of dose start associated with the current interval (values start at 0 at the beginning of the current interval).}
#'     \item{"duration.dose"}{Duration of dose (typically infusion duration) for doses in the current interval.}
#'     \item{"route"}{Route of dosing for the current interval.}
#'     \item{"start"}{Time of interval start.}
#'     \item{"end"}{Time of interval end.}
#'     \item{"options"}{PKNCA.options governing calculations.}
#'   }
#'   \item{For the current group:}
#'   \describe{
#'     \item{"conc.group"}{Concentration measurements for the current group.}
#'     \item{"time.group"}{Times associated with concentration measurements for the current group (values start at 0 at the beginning of the current interval).}
#'     \item{"volume.group"}{Volume associated with concentration measurements for the current interval (typically applies for excretion parameters like urine).}
#'     \item{"duration.conc.group"}{Durations assocuated with concentration measurements for the current group.}
#'     \item{"dose.group"}{Dose amounts assocuated with the current group.}
#'     \item{"time.dose.group"}{Time of dose start associated with the current group (values start at 0 at the beginning of the current interval).}
#'     \item{"duration.dose.group"}{Duration of dose (typically infusion duration) for doses in the current group.}
#'     \item{"route.group"}{Route of dosing for the current group.}
#'   }
#' }
#' @examples
#' \dontrun{
#' add.interval.col("cmax",
#'                  FUN="pk.calc.cmax",
#'                  values=c(FALSE, TRUE),
#'                  unit_type="conc",
#'                  pretty_name="Cmax",
#'                  desc="Maximum observed concentration")
#' add.interval.col("cmax.dn",
#'                  FUN="pk.calc.dn",
#'                  values=c(FALSE, TRUE),
#'                  unit_type="conc_dosenorm",
#'                  pretty_name="Cmax (dose-normalized)",
#'                  desc="Maximum observed concentration, dose normalized",
#'                  formalsmap=list(parameter="cmax"),
#'                  depends="cmax")
#' }
#' @family Interval specifications
#' @export
add.interval.col <- function(name,
                             FUN,
                             values=c(FALSE, TRUE),
                             unit_type,
                             pretty_name,
                             depends=NULL,
                             desc="",
                             sparse=FALSE,
                             formalsmap=list(),
                             datatype=c("interval",
                               "individual",
                               "population")) {
  # Check inputs
  if (!is.character(name)) {
    stop("name must be a character string")
  } else if (length(name) != 1) {
    stop("name must have length == 1")
  }
  if (length(FUN) != 1) {
    stop("FUN must have length == 1")
  } else if (!(is.character(FUN) | is.na(FUN))) {
    stop("FUN must be a character string or NA")
  }
  if (!is.null(depends)) {
    if (!is.character(depends)) {
      stop("'depends' must be NULL or a character vector")
    }
  }
  checkmate::assert_logical(sparse, any.missing=FALSE, len=1)
  unit_type <-
    match.arg(
      unit_type,
      choices=c(
        "unitless", "fraction", "%", "count",
        "time", "inverse_time",
        "amount",
        "conc", "conc_dosenorm",
        "dose",
        "volume",
        "auc", "aumc",
        "auc_dosenorm", "aumc_dosenorm",
        "clearance", "renal_clearance"
      )
    )
  stopifnot("pretty_name must be a scalar"=length(pretty_name) == 1)
  stopifnot("pretty_name must be a character"=is.character(pretty_name))
  stopifnot("pretty_name must not be an empty string"=nchar(pretty_name) > 0)
  datatype <- match.arg(datatype)
  if (!(datatype %in% "interval")) {
    stop("Only the 'interval' datatype is currently supported.")
  }
  if (length(desc) != 1) {
    stop("desc must have length == 1")
  } else if (!is.character(desc)) {
    stop("desc must be a character string")
  }
  if (!is.list(formalsmap)) {
    stop("formalsmap must be a list")
  } else if (length(formalsmap) > 0 &
             is.null(names(formalsmap))) {
    stop("formalsmap must be a named list")
  } else if (length(formalsmap) > 0 &
             is.na(FUN)) {
    stop("formalsmap may not be given when FUN is NA.")
  } else if (!all(nchar(names(formalsmap)) > 0)) {
    stop("All formalsmap elements must be named")
  }
  # Ensure that the function exists
  if (!is.na(FUN) &&
      length(utils::getAnywhere(FUN)$objs) == 0) {
    stop("The function named '", FUN, "' is not defined.  Please define the function before calling add.interval.col.")
  }
  if (!is.na(FUN) &
      length(formalsmap) > 0) {
    # Ensure that the formalsmap parameters are all in the list of
    # formal arguments to the function.
    if (!all(names(formalsmap) %in% names(formals(utils::getAnywhere(FUN)$objs[[1]])))) {
      stop("All names for the formalsmap list must be arguments to the function.")
    }
  }
  current <- get("interval.cols", envir=.PKNCAEnv)
  current[[name]] <-
    list(
      FUN=FUN,
      values=values,
      unit_type=unit_type,
      pretty_name=pretty_name,
      desc=desc,
      sparse=sparse,
      formalsmap=formalsmap,
      depends=depends,
      datatype=datatype
    )
  assign("interval.cols", current, envir=.PKNCAEnv)
}

#' Sort the interval columns by dependencies.
#'
#' Columns are always to the right of columns that they depend on.
sort.interval.cols <- function() {
  current <- get("interval.cols", envir=.PKNCAEnv)
  # Only sort if necessary
  sort_order <- get0("interval.cols_sorted", envir=.PKNCAEnv)
  if (identical(sort_order, names(current))) {
    # It is already sorted
    return(sort_order)
  }
  # Build a dependency tree
  myorder <- rep(NA, length(current))
  names(myorder) <- names(current)
  nextnum <- 1
  while (any(is.na(myorder))) {
    for (nextorder in seq_along(myorder)[is.na(myorder)]) {
      if (length(current[[nextorder]]$depends) == 0) {
        # If it doesn't depend on anything then it can go next in order.
        myorder[nextorder] <- nextnum
        nextnum <- nextnum + 1
      } else {
        # If all of its dependencies already have values, then it can be next.
        deps <- unique(unlist(current[[nextorder]]$depends))
        missing_deps <- deps[!(deps %in% names(myorder))]
        if (length(missing_deps) > 0) {
          stop(
            "Invalid dependencies for interval column (please report this as a bug): ",
            names(myorder)[nextorder],
            " The following dependencies are missing: ",
            paste(missing_deps, collapse=", ")
          )
        }
        if (!any(is.na(myorder[deps]))) {
          myorder[nextorder] <- nextnum
          nextnum <- nextnum + 1
        }
      }
    }
  }
  current <- current[names(sort(myorder))]
  assign("interval.cols_sorted", names(current), envir=.PKNCAEnv)
  assign("interval.cols", current, envir=.PKNCAEnv)
  invisible(myorder)
}

#' Get the columns that can be used in an interval specification
#'
#' @param out_format What output format should be provided?
#' @returns If `out_format = "list"`, a list with named elements for each
#'   parameter. Each list element contains the parameter definition. If
#'   `out_format = "sdtm_map"`, a data.frame with the PKNCA PPTESTCD_PKNCA, and
#'   the SDTM PPTEST and PPTESTCD.
#' @seealso [check.interval.specification()] and the vignette "Selection of
#'   Calculation Intervals"
#' @examples
#' get.interval.cols()
#' @family Interval specifications
#' @export
get.interval.cols <- function(out_format = c("list", "sdtm_map")) {
  out_format <- match.arg(out_format)
  sort.interval.cols()
  ret <- get("interval.cols", envir=.PKNCAEnv)
  if (out_format == "sdtm_map") {
    browser()
    stop()
    ret_df <-
      data.frame(
        PPTESTCD_PKNCA = names(ret),
        PPTESTCD = vapply(X = ret, FUN = function(x) x$SDTM_PPTESTCD
      )
  } else if (out_format == "list") {
    # do nothing
  } else {
    stop("Please report a bug: unknown out_format, ", out_format) # nocov
  }
  ret
}

# Add the start and end interval columns
add.interval.col(
  "start",
  FUN = NA,
  values = as.numeric,
  unit_type="time",
  pretty_name="Interval Start",
  desc = "Starting time of the interval"
)
add.interval.col(
  "end",
  FUN = NA,
  values = as.numeric,
  unit_type="time",
  pretty_name="Interval End",
  desc = "Ending time of the interval (potentially infinity)"
)
