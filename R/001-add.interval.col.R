## Setup the default options
.PKNCAEnv <- new.env(parent=emptyenv())
assign("options", NULL, envir=.PKNCAEnv)
assign("summary", list(), envir=.PKNCAEnv)
assign("interval.cols", list(), envir=.PKNCAEnv)

#' Add columns for calculations within PKNCA intervals
#'
#' @param name The column name as a character string
#' @param FUN The function to run (as a character string) or \code{NA}
#'   if the parameter is automatically calculated when calculating
#'   another parameter.
#' @param values Valid values for the column
#' @param depends Character vector of columns that must be run before
#'   this column.
#' @param desc A human-readable description of the parameter (<=40
#'   characters to comply with SDTM)
#' @param formalsmap A named list mapping parameter names in the
#'   function call to NCA parameter names.  See the details for information on use of \code{formalsmap}.
#' @param datatype The type of data used for the calculation
#' @return NULL (Calling this function has a side effect of
#'   changing the available intervals for calculations)
#'
#' @details
#' The \code{formalsmap} argument enables mapping some alternate formal
#' argument names to parameters.  It is used to generalize functions
#' that may use multiple similar arguments (such as the variants of mean
#' residence time). The names of the list should correspond to function
#' formal parameter names and the values should be one of the following:
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
#'                  desc="Maximum observed concentration",
#'                  depends=c())
#' add.interval.col("cmax.dn",
#'                  FUN="pk.calc.dn",
#'                  values=c(FALSE, TRUE),
#'                  desc="Maximum observed concentration, dose normalized",
#'                  formalsmap=list(parameter="cmax"),
#'                  depends=c("cmax"))
#' }
#' @importFrom utils getAnywhere
#' @family Interval specifications
#' @export
add.interval.col <- function(name,
                             FUN,
                             values=c(FALSE, TRUE),
                             depends=c(),
                             desc="",
                             formalsmap=list(),
                             datatype=c("interval",
                               "individual",
                               "population")) {
  ## Check inputs
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
  ## Ensure that the function exists
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
      desc=desc,
      formalsmap=formalsmap,
      depends=depends,
      datatype=datatype
    )
  assign("interval.cols", current, envir=.PKNCAEnv)
}

## Add the start and end interval columns
add.interval.col("start",
  FUN = NA,
  values = as.numeric,
  desc = "Starting time of the interval"
)
add.interval.col("end",
  FUN = NA,
  values = as.numeric,
  desc = "Ending time of the interval (potentially infinity)"
)

#' Sort the interval columns by dependencies.
#'
#' Columns are always to the right of columns that they depend on.
sort.interval.cols <- function() {
  current <- get("interval.cols", envir=.PKNCAEnv)
  ## Build a dependency tree
  myorder <- rep(NA, length(current))
  names(myorder) <- names(current)
  nextnum <- 1
  while (any(is.na(myorder)))
    for (nextorder in (1:length(myorder))[is.na(myorder)]) {
      if (length(current[[nextorder]]$depends) == 0) {
        ## If it doesn't depend on anything then it can go next in
        ## order.
        myorder[nextorder] <- nextnum
        nextnum <- nextnum + 1
      } else {
        ## If all of its dependencies already have values, then it can
        ## be next.
        deps <- unique(unlist(current[[nextorder]]$depends))
        if (!all(deps %in% names(myorder)))
          stop("Invalid dependencies for interval column (please report this as a bug):", # nocov
               names(myorder)[nextorder]) # nocov
        if (!any(is.na(myorder[deps]))) {
          myorder[nextorder] <- nextnum
          nextnum <- nextnum + 1
        }
      }
    }
  current <- current[names(sort(myorder))]
  assign("interval.cols", current, envir=.PKNCAEnv)
  invisible(myorder)
}

#' Get the columns that can be used in an interval specification
#' @return A list with named elements for each parameter.  Each list element
#'   contains the parameter definition.
#' @seealso \code{\link{check.interval.specification}()} and the vignette
#'   "Selection of Calculation Intervals"
#' @examples
#' get.interval.cols()
#' @family Interval specifications
#' @export
get.interval.cols <- function() {
  sort.interval.cols()
  get("interval.cols", envir=.PKNCAEnv)
}
