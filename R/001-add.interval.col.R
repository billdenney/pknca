## Setup the default options
.PKNCAEnv <- new.env(parent=emptyenv())
assign("options", NULL, envir=.PKNCAEnv)
assign("summary", list(), envir=.PKNCAEnv)
assign("interval.cols", list(), envir=.PKNCAEnv)

#' Add columns for calculations within PKNCA intervals
#'
#' @param name The column name
#' @param FUN The function to run (as a character string)
#' @param values Valid values for the column
#' @param depends Character vector of columns that must be run before
#' this column.
#' @param desc A human-readable description of the parameter (<=40
#' characters to comply with SDTM)
#' @param datatype The type of data 
add.interval.col <- function(name, FUN, values, depends=c(),
                             desc="",
                             datatype=c("interval",
                               "individual",
                               "population")) {
  current <- get("interval.cols", envir=.PKNCAEnv)
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
                 desc="Starting time of the interval",
                 depends=c())
add.interval.col("end",
                 FUN=NA,
                 values=as.numeric,
                 desc="Ending time of the interval (potentially infinity)",
                 depends=c())

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
#' @export
get.interval.cols <- function() {
  sort.interval.cols()
  get("interval.cols", envir=.PKNCAEnv)
}
