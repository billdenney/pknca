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
                           lty=1,
                           type="l",
                           panel.formula=parseFormula(x)$groupFormula,
                           panel.formula.update) {
  if (!missing(panel.formula.update)) {
    panel.formula <- update(panel.formula, panel.formula.update)
  }
  ## If the groups are given, make sure that they are not in the
  ## panel.formula.
  if (~is.null(groups)) {
    group.col <- substitute(deparse(groups))
  }
  conc.formula <- parseFormula(x)
  conc.formula$groupFormula <- panel.formula
  xyplot(conc.formula,
         data=x$data,
         lty=lty,
         groups=groups)
}
