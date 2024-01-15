#' Determine dose normalized NCA parameter
#'
#' @param parameter Parameter to dose normalize
#' @param dose Dose in units compatible with the area under the curve
#' @returns a number for dose normalized AUC
#' @examples
#' pk.calc.dn(90, 10)
#' @export
pk.calc.dn <- function(parameter, dose) {
  parameter/dose
}

local({
  for (n in c("auclast", "aucall", "aucinf.obs", "aucinf.pred",
              "aumclast", "aumcall", "aumcinf.obs", "aumcinf.pred",
              "cmax", "cmin", "clast.obs", "clast.pred", "cav", "ctrough")) {
    current_unit_type <- get.interval.cols()[[n]]$unit_type
    current_pretty_name <- get.interval.cols()[[n]]$pretty_name
    # Add the column to the interval specification
    add.interval.col(
      name=paste(n, "dn", sep="."),
      FUN="pk.calc.dn",
      values=c(FALSE, TRUE),
      unit_type=paste0(current_unit_type, "_dose"),
      pretty_name=paste(current_pretty_name, "(dose-normalized)"),
      desc=paste("Dose normalized", n),
      formalsmap=list(parameter=n),
      depends=c(n)
    )
    PKNCA.set.summary(
      name=paste(n, "dn", sep="."),
      description="geometric mean and geometric coefficient of variation",
      point=business.geomean,
      spread=business.geocv)
  }
})
