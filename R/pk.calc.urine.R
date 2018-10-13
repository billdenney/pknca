#' Calculate amount excreted (typically in urine or feces)
#'
#' @details ae is \code{sum(conc*volume)}.
#'  
#' @param conc The concentration in the sample
#' @param volume The volume (or mass) of the sample
#' @param check Should the concentration and volume data be checked?
#' @return The amount excreted during the interval
#' @details The units for the concentration and volume should match such
#'   that \code{sum(conc*volume)} has units of mass or moles.
#' @seealso \code{\link{pk.calc.clr}}, \code{\link{pk.calc.fe}}
#' @export
pk.calc.ae <- function(conc, volume, check=TRUE) {
  sum(conc*volume)
}
add.interval.col("ae",
                 FUN="pk.calc.ae",
                 values=c(FALSE, TRUE),
                 desc="The amount excreted (typically into urine or feces)")
PKNCA.set.summary(
  name="ae",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate renal clearance
#' 
#' @details clr is \code{sum(ae)/auc}.
#' 
#' @param ae The amount excreted in urine (as a numeric scalar or
#'   vector)
#' @param auc The area under the curve (as a numeric scalar or vector)
#' @return The renal clearance as a number
#' @details The units for the \code{ae} and \code{auc} should match such
#'   that \code{ae/auc} has units of volume/time.
#' @seealso \code{\link{pk.calc.ae}}, \code{\link{pk.calc.fe}}
#' @export
pk.calc.clr <- function(ae, auc) {
  sum(ae)/auc
}
add.interval.col("clr.last",
                 FUN="pk.calc.clr",
                 values=c(FALSE, TRUE),
                 formalsmap=list(auc="auclast"),
                 desc="The renal clearance calculated using AUClast")
PKNCA.set.summary(
  name="clr.last",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("clr.obs",
                 FUN="pk.calc.clr",
                 values=c(FALSE, TRUE),
                 formalsmap=list(auc="aucinf.obs"),
                 desc="The renal clearance calculated using AUCinf,obs")
PKNCA.set.summary(
  name="clr.obs",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
add.interval.col("clr.pred",
                 FUN="pk.calc.clr",
                 values=c(FALSE, TRUE),
                 formalsmap=list(auc="aucinf.pred"),
                 desc="The renal clearance calculated using AUCinf,pred")
PKNCA.set.summary(
  name="clr.pred",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' Calculate fraction excreted (typically in urine or feces)
#' 
#' @details fe is \code{sum(ae)/dose}
#'
#' @param ae The amount excreted (as a numeric scalar or vector)
#' @param dose The dose (as a numeric scalar or vector)
#' @return The fraction of dose excreted.
#' @details   The units for \code{ae} and \code{dose} should be the same
#'   so that \code{ae/dose} is a unitless fraction.
#' @seealso \code{\link{pk.calc.ae}}, \code{\link{pk.calc.clr}}
#' @export
pk.calc.fe <- function(ae, dose) {
  sum(ae)/dose
}
add.interval.col("fe",
                 FUN="pk.calc.fe",
                 values=c(FALSE, TRUE),
                 desc="The fraction of the dose excreted")
PKNCA.set.summary(
  name="fe",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
