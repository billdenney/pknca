# Generics ####

#' Convert an object into a PKNCAconc object
#'
#' @param x The object to convert
#' @param ... Passed to subsequent methods
#' @return A converted object
#' @export
as_PKNCAconc <- function(x, ...) {
  UseMethod("as_PKNCAconc")
}

#' @describeIn as_PKNCAconc Convert an object into a PKNCAdose object
#' @export
as_PKNCAdose <- function(x, ...) {
  UseMethod("as_PKNCAdose")
}

#' @describeIn as_PKNCAconc Convert an object into a PKNCAdata object
#' @export
as_PKNCAdata <- function(x, ...) {
  UseMethod("as_PKNCAdata")
}

#' @describeIn as_PKNCAconc Convert an object into a PKNCAresults object
#' @export
as_PKNCAresults <- function(x, ...) {
  UseMethod("as_PKNCAresults")
}

# as_PKNCAconc ####

#' @export
as_PKNCAconc.PKNCAresults <- function(x, ...) {
  as_PKNCAconc(as_PKNCAdata(x, ...), ...)
}

#' @export
as_PKNCAconc.PKNCAdata <- function(x, ...) {
  x$conc
}

#' @export
as_PKNCAconc.PKNCAconc <- function(x, ...) {
  x
}

# as_PKNCAdose ####

#' @export
as_PKNCAdose.PKNCAresults <- function(x, ...) {
  as_PKNCAdose(as_PKNCAdata(x, ...), ...)
}

#' @export
as_PKNCAdose.PKNCAdata <- function(x, ...) {
  x$dose
}

#' @export
as_PKNCAdose.PKNCAdose <- function(x, ...) {
  x
}

# as_PKNCAdata ####

#' @export
as_PKNCAdata.PKNCAresults <- function(x, ...) {
  x$data
}

#' @export
as_PKNCAdata.PKNCAdata <- function(x, ...) {
  x
}

# as_PKNCAresults ####

#' @export
as_PKNCAresults.PKNCAresults <- function(x, ...) {
  x
}

#' @export
as_PKNCAresults.PKNCAdata <- function(x, ...) {
  # the conversion from PKNCAdata to PKNCAresults is to do the calculation
  pk.nca(x)
}
