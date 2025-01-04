#' Calculate AUC for intravenous dosing
#'
#' @details The AUC for intravenous (IV) dosing extrapolates the AUC back from
#'   the first measurement to time 0 using c0 and the AUC calculated by another
#'   method (for example the auclast).
#'
#' The calculation method takes the following steps:
#'
#' \itemize{
#'   \item{`time = 0` must be present in the data with a measured concentration.}
#'   \item{The AUC between `time = 0` and the next time point is calculated (`auc_first`).}
#'   \item{The AUC between `time = 0` with `c0` and the next time point is calculated (`auc_second`).}
#'   \item{The final AUC is the initial AUC plus the difference between the two AUCs (`auc_final <- auc + auc_second - auc_first`).}
#' }
#' @inheritParams pk.calc.auxc
#' @inheritParams PKNCA.choose.option
#' @param c0 The concentration at time 0, typically calculated using
#'   `pk.calc.c0()`
#' @param auc The AUC calculated using `conc` and `time` without `c0` (it may be
#'   calculated using any method)
#' @return `pk.calc.auciv`: The AUC calculated using `c0`
#' @export
pk.calc.auciv <- function(conc, time, c0, auc, ..., options = list(), check=TRUE) {
  if (check) {
    assert_conc_time(conc = conc, time = time)
    data <-
      clean.conc.blq(
        conc, time,
        options = options,
        check=FALSE
      )
  } else {
    data <- data.frame(conc = conc, time = time)
  }
  if (!(0 %in% time)) {
    return(structure(NA_real_, exclude="No time 0 in data"))
  } else if (is.na(c0)) {
    return(structure(NA_real_, exclude="c0 is not calculated"))
  }
  auc_first <- pk.calc.auc.last(conc = data$conc[1:2], time = data$time[1:2], ..., check=FALSE)
  auc_second <- pk.calc.auc.last(conc = c(c0, data$conc[2]), time = data$time[1:2], ..., check=FALSE)
  auc_final <- auc + auc_second - auc_first
  auc_final
}

add.interval.col(
  name = "aucivlast",
  FUN = "pk.calc.auciv",
  unit_type = "auc",
  pretty_name = "AUClast (IV dosing)",
  depends = c("auclast", "c0"),
  desc = "The AUClast calculated with back-extrapolation for intravenous dosing using extrapolated C0",
  sparse = FALSE,
  formalsmap = list(auc="auclast")
)

add.interval.col(
  name = "aucivall",
  FUN = "pk.calc.auciv",
  unit_type = "auc",
  pretty_name = "AUCall (IV dosing)",
  depends = c("aucall", "c0"),
  desc = "The AUCall calculated with back-extrapolation for intravenous dosing using extrapolated C0",
  sparse = FALSE,
  formalsmap = list(auc="aucall")
)

add.interval.col(
  name = "aucivint.last",
  FUN = "pk.calc.auciv",
  unit_type = "auc",
  pretty_name = "AUCint,last (IV dosing)",
  depends = c("aucint.last", "c0"),
  desc = "The AUCint,last calculated with back-extrapolation for intravenous dosing using extrapolated C0",
  sparse = FALSE,
  formalsmap = list(auc="aucint.last")
)

add.interval.col(
  name = "aucivint.all",
  FUN = "pk.calc.auciv",
  unit_type = "auc",
  pretty_name = "AUCint,all (IV dosing)",
  depends = c("aucint.all", "c0"),
  desc = "The AUCint,all calculated with back-extrapolation for intravenous dosing using extrapolated C0",
  sparse = FALSE,
  formalsmap = list(auc="aucint.all")
)

add.interval.col(
  name = "aucivinf.obs",
  FUN = "pk.calc.auciv",
  unit_type = "auc",
  pretty_name = "AUCinf,obs (IV dosing)",
  depends = c("aucinf.obs", "c0"),
  desc = "The AUCinf,obs calculated with back-extrapolation for intravenous dosing using extrapolated C0",
  sparse = FALSE,
  formalsmap = list(auc="aucinf.obs")
)

add.interval.col(
  name = "aucivinf.pred",
  FUN = "pk.calc.auciv",
  unit_type = "auc",
  pretty_name = "AUCinf,pred (IV dosing)",
  depends = c("aucinf.pred", "c0"),
  desc = "The  calculated with back-extrapolation for intravenous dosing using extrapolated C0",
  sparse = FALSE,
  formalsmap = list(auc="aucinf.pred")
)

PKNCA.set.summary(
  name=c("aucivlast", "aucivall", "aucivint.last", "aucivint.all", "aucivinf.obs", "aucivinf.pred"),
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)

#' @describeIn pk.calc.auciv Calculate the percent back-extrapolated AUC for IV
#'   administration
#' @details The calculation for back-extrapolation is `100*(1 - auc/auciv)`.
#' @param auciv The AUC calculated using `c0`
#' @returns `pk.calc.auciv_pctbackextrap`: The AUC percent back-extrapolated
#' @export
pk.calc.auciv_pbext <- function(auc, auciv) {
  100*(1 - auc/auciv)
}

add.interval.col(
  name = "aucivpbextlast",
  FUN = "pk.calc.auciv_pbext",
  unit_type = "%",
  pretty_name = "AUCbext (based on AUClast)",
  depends = c("auclast", "aucivlast"),
  desc = "The back-extrapolation percent for intravenous dosing based on AUClast",
  sparse = FALSE,
  formalsmap = list(auc="auclast", auciv="aucivlast")
)

add.interval.col(
  name = "aucivpbextall",
  FUN = "pk.calc.auciv_pbext",
  unit_type = "%",
  pretty_name = "AUCbext (based on AUCall)",
  depends = c("aucall", "aucivall"),
  desc = "The back-extrapolation percent for intravenous dosing based on AUCall",
  sparse = FALSE,
  formalsmap = list(auc="aucall", auciv="aucivall")
)

add.interval.col(
  name = "aucivpbextint.last",
  FUN = "pk.calc.auciv_pbext",
  unit_type = "%",
  pretty_name = "AUCbext (based on AUCint,last)",
  depends = c("aucint.last", "aucivint.last"),
  desc = "The back-extrapolation percent for intravenous dosing based on AUCint,last",
  sparse = FALSE,
  formalsmap = list(auc="aucint.last", auciv="aucivint.last")
)

add.interval.col(
  name = "aucivpbextint.all",
  FUN = "pk.calc.auciv_pbext",
  unit_type = "%",
  pretty_name = "AUCbext (based on AUCint,all)",
  depends = c("aucint.all", "aucivint.all"),
  desc = "The back-extrapolation percent for intravenous dosing based on AUCint,all",
  sparse = FALSE,
  formalsmap = list(auc="aucint.all", auciv="aucivint.all")
)

add.interval.col(
  name = "aucivpbextinf.obs",
  FUN = "pk.calc.auciv_pbext",
  unit_type = "%",
  pretty_name = "AUCbext (based on AUCinf,obs)",
  depends = c("aucinf.obs", "aucivinf.obs"),
  desc = "The back-extrapolation percent for intravenous dosing based on AUCinf,obs",
  sparse = FALSE,
  formalsmap = list(auc="aucinf.obs", auciv="aucivinf.obs")
)

add.interval.col(
  name = "aucivpbextinf.pred",
  FUN = "pk.calc.auciv_pbext",
  unit_type = "%",
  pretty_name = "AUCbext (based on AUCinf,pred)",
  depends = c("aucinf.pred", "aucivinf.pred"),
  desc = "The back-extrapolation percent for intravenous dosing based on AUCinf,pred",
  sparse = FALSE,
  formalsmap = list(auc="aucinf.pred", auciv="aucivinf.pred")
)

PKNCA.set.summary(
  name = c("aucivpbextlast", "aucivpbextall", "aucivpbextint.last", "aucivpbextint.all", "aucivpbextinf.obs", "aucivpbextinf.pred"),
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)
