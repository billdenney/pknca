#' Plot a PKNCAconc object
#'
#' @param x The object to plot
#' @param groups The grouping variable for the plot (typically the
#' subject column)
#' @param \dots Additional arguments passed to \code{xyplot}
#' @param panel.formula The formula used for the call to xyplot
#' (defaults to the group formula of \code{x})
#' @param panel.formula.update Updates to the \code{panel.formula} to
#' simplify modifications without having to fully specify the formula.
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
    panel.formula <- stats::update(panel.formula, panel.formula.update)
  }
  conc.formula <- parseFormula(x)
  ## If the groups are given, make sure that they are not in the
  ## panel.formula.
  if (!is.null(groups)) {
    conc.formula$groupFormula <-
      stats::update(panel.formula,
                    stats::as.formula(sprintf(".~-%s", groups)))
  }
  call.args[["x"]] <- stats::formula(conc.formula)
  call.args[["data"]] <- x$data
  ## If labels and/or units are given for the x and y variables, use
  ## them.
  xlab <- make.label("rhs", x$data, conc.formula, x$labels, x$units)
  ylab <- make.label("lhs", x$data, conc.formula, x$labels, x$units)
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
  graphics::plot(x$conc, ...)
