## Setup the default options
.PKNCAEnv <- new.env(parent=emptyenv())
assign("options", NULL, envir=.PKNCAEnv)
assign("summary", list(), envir=.PKNCAEnv)
assign("interval.cols", list(), envir=.PKNCAEnv)

#' Add columns for calculations within PKNCA intervals
#' 
#' @param name The column name as a character string
#' @param FUN The function to run (as a character string) or \code{NA} if the
#'   parameter is automatically calculated when calculating another parameter.
#' @param values Valid values for the column
#' @param depends Character vector of columns that must be run before this
#'   column.
#' @param desc A human-readable description of the parameter (<=40 characters to
#'   comply with SDTM)
#' @param datatype The type of data used for the calculation
#' @return \code{NULL} (changes the available intervals for calculations
#' @importFrom utils getAnywhere
add.interval.col <- function(name,
                             FUN,
                             values=c(FALSE, TRUE),
                             depends=c(),
                             desc="",
                             datatype=c("interval",
                               "individual",
                               "population")) {
  ## Check inputs
  if (!is.character(name)) {
    stop("name must be a character string")
  } else if (length(name) != 1) {
    stop("name must have length == 1")
  }
  if (!(is.character(FUN) | is.na(FUN))) {
    stop("FUN must be a character string or NA")
  } else if (length(FUN) != 1) {
    stop("FUN must have length == 1")
  }
  datatype <- match.arg(datatype)
  if (!(datatype %in% "interval")) {
    stop("Only the 'interval' datatype is currently supported.")
  }
  if (!is.character(desc)) {
    stop("desc must be a character string")
  } else if (length(desc) != 1) {
    stop("desc must have length == 1")
  }
  current <- get("interval.cols", envir=.PKNCAEnv)
  ## Ensure that the function exists
  if (length(utils::getAnywhere(FUN)) == 0)
    stop("The function named '", FUN, "' is not defined.  Please define the function before calling add.interval.col.")
  current[[name]] <- list(FUN=FUN,
                          values=values,
                          desc=desc,
                          depends=depends,
                          datatype=match.arg(datatype))
  assign("interval.cols", current, envir=.PKNCAEnv)
}

## Add the start and end interval columns
add.interval.col("start",
                 FUN=NA,
                 values=as.numeric,
                 desc="Starting time of the interval")
add.interval.col("end",
                 FUN=NA,
                 values=as.numeric,
                 desc="Ending time of the interval (potentially infinity)")

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
          stop("Invalid dependencies for interval column:",
               names(myorder)[nextorder])
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
#' @examples
#' get.interval.cols()
#' @export
get.interval.cols <- function() {
  sort.interval.cols()
  get("interval.cols", envir=.PKNCAEnv)
}
