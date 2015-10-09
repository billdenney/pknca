#' Plot a PKNCAconc object
#'
#' @param x The object to plot
#' @param groups The grouping variable for the plot (typically the
#' subject column)
#' @param \dots Additional arguments passed to \code{xyplot}
#' @return A trellis object of the plot(s)
#' @export
plot.PKNCAconc <- function(x, ...,
                           groups=x$subject,
                           panel.formula=parseFormula(x)$groupFormula,
                           panel.formula.update) {
  ## Arguments that are set by default only overwrite if the user
  ## doesn't specify others with the same name.
  default.args <- list(lty=1, type="b", auto.key=list(columns=2),
                       scales=list(alternating=FALSE))
  call.args <- list(...)
  set.defaults <- setdiff(names(default.args), names(call.args))
  call.args[set.defaults] <- default.args[set.defaults]
  ## Update the panel formula if applicable
  if (!missing(panel.formula.update)) {
    panel.formula <- update(panel.formula, panel.formula.update)
  }
  conc.formula <- parseFormula(x)
  ## If the groups are given, make sure that they are not in the
  ## panel.formula.
  if (!is.null(groups)) {
    conc.formula$groupFormula <-
      update(panel.formula,
             as.formula(sprintf(".~-%s", groups)))
  }
  call.args[["x"]] <- formula(conc.formula)
  call.args[["data"]] <- x$data
  ## If labels and/or units are given for the x and y variables, use
  ## them.
  xlab <- make.label("rhs", data, conc.formula, x$labels, x$units)
  ylab <- make.label("lhs", data, conc.formula, x$labels, x$units)
  if (!("xlab" %in% names(call.args)))
    call.args[["xlab"]] <- xlab
  if (!("ylab" %in% names(call.args)))
    call.args[["ylab"]] <- ylab
  do.call(lattice::xyplot, call.args)
}

## Make plotting labels from data, a formula, labels, and units.
make.label <- function(side, data, parsed.formula, labels, units) {
  col.text <- as.character(parsed.formula[[side]])
  label <- col.text
  if (!is.null(labels))
    if (col.text %in% names(labels))
      label <- labels[[col.text]]
  if (!is.null(units))
    if (col.text %in% names(units))
      label <- sprintf("%s (%s)", label, units[[col.text]])
  label
}

#' @rdname plot.PKNCAconc
#' @export
plot.PKNCAdata <- function(x, ...)
  plot(x$conc, ...)
