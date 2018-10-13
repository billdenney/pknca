#' Determine dose normalized NCA parameter
#'
#' @param parameter Parameter to dose normalize
#' @param dose Dose in units compatible with the area under the curve
#' @return a number for dose normalized AUC
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
    ## Add the column to the interval specification
    add.interval.col(
      name=paste(n, "dn", sep="."),
      FUN="pk.calc.dn",
      values=c(FALSE, TRUE),
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
