#' Create a unit assignment and conversion table
#' 
#' This data.frame is typically used for the \code{units} argument for
#' \code{\link{PKNCAdata}()}.
#' 
#' @param concu,doseu,amountu,timeu Units for concentration, dose, amount, and
#'   time
#' @param conversions An optional data.frame with columns of c("PPORRESU",
#'   "PPSTRESU", "conversion_factor") for the original calculation units, the
#'   standardized units, and a conversion factor to multiply the initial value
#'   by to get a standardized value.
#' @return A unit conversion table with columns for "PPTESTCD" and "PPORRESU" if
#'   \code{conversions} is not given, and adding "PPSTRESU" and
#'   "conversion_factor" if \code{conversions} is given.
#' @seealso The \code{units} argument for \code{\link{PKNCAdata}()}
#' @examples
#' pknca_units_table()
#' @export
pknca_units_table <- function(concu, doseu, amountu, timeu, conversions=data.frame()) {
  # The unit conversions are grouped by the type of inputs required
  ret_any <-
    rbind(
    )
  ret <-
    rbind(
      pknca_units_table_unitless(),
      pknca_units_table_time(timeu=timeu),
      pknca_units_table_conc(concu=concu),
      pknca_units_table_amount(amountu=amountu),
      pknca_units_table_conc_dose(concu=concu, doseu=doseu),
      pknca_units_table_conc_time(concu=concu, timeu=timeu),
      pknca_units_table_conc_time_dose(concu=concu, timeu=timeu, doseu=doseu),
      pknca_units_table_conc_time_amount(concu=concu, timeu=timeu, amountu=amountu)
    )
  # You don't have to define parameters for everything for the parameters to be useful
  # missing_cols <- setdiff(names(get.interval.cols()), c(ret$PPTESTCD, "start", "end"))
  # if (length(missing_cols) > 0) {
  #   stop("The following NCA parameters do not have units defined: ", paste(missing_cols, collapse=", "))
  # }
  extra_cols <- setdiff(ret$PPTESTCD, names(PKNCA::get.interval.cols()))
  if (length(extra_cols) > 0) {
    stop("The unknown NCA parameters have units defined: ", paste(extra_cols, collapse=", "))
  }
  if (nrow(conversions) > 0) {
    stopifnot(!duplicated(conversions$PPORRESU))
    stopifnot(!duplicated(conversions$PPSTRESU))
    stopifnot(length(setdiff(names(conversions), c("PPORRESU", "PPSTRESU", "conversion_factor"))) == 0)
    if (!("conversion_factor" %in% names(conversions))) {
      if (!requireNamespace("units", quietly=TRUE)) {
        stop("The units package is required for automatic unit conversion")
      }
      conversions$conversion_factor <- NA_real_
    }
    for (idx in which(is.na(conversions$conversion_factor))) {
      conversions$conversion_factor[idx] <-
        as.numeric(
          units::set_units(
            units::set_units(
              1,
              conversions$PPORRESU[idx], mode="standard"
            ),
            conversions$PPSTRESU[idx], mode="standard"
          )
        )
    }
    ret <- dplyr::left_join(ret, conversions)
    # anything that does not have a specific conversion factor is assumed to be
    # in the correct units, and a unit conversion factor is used to keep the
    # units the same.
    ret$PPSTRESU[is.na(ret$PPSTRESU)] <- ret$PPORRESU[is.na(ret$PPSTRESU)]
    ret$conversion_factor[is.na(ret$conversion_factor)] <- 1
  }
  ret
}

pknca_units_table_unitless <- function() {
  rbind(
    data.frame(
      PPORRESU="unitless",
      PPTESTCD=c("adj.r.squared", "r.squared"),
      stringsAsFactors=FALSE
    ),
    data.frame(
      PPORRESU="fraction",
      PPTESTCD=c("f", "fe", "ptr", "span.ratio"),
      stringsAsFactors=FALSE
    ),
    data.frame(
      PPORRESU="%",
      PPTESTCD=c("aucpext.obs", "aucpext.pred", "deg.fluc", "swing"),
      stringsAsFactors=FALSE
    ),
    data.frame(
      PPORRESU="count",
      PPTESTCD="lambda.z.n.points",
      stringsAsFactors=FALSE
    )
  )
}

pknca_units_table_time <- function(timeu) {
  if (!missing(timeu)) {
    rbind(
      data.frame(
        PPORRESU=timeu,
        PPTESTCD=
          c(
            "half.life", "lambda.z.time.first",
            "mrt.iv.last", "mrt.iv.obs", "mrt.iv.pred", "mrt.last", "mrt.md.obs", "mrt.md.pred", "mrt.obs", "mrt.pred",
            "tfirst", "thalf.eff.iv.last", "thalf.eff.iv.obs", "thalf.eff.iv.pred", 
            "thalf.eff.last", "thalf.eff.obs", "thalf.eff.pred", "time_above", 
            "tlag", "tlast", "tmax"
          ),
        stringsAsFactors=FALSE
      ),
      data.frame(
        PPORRESU=sprintf("1/%s", timeu),
        PPTESTCD=c("kel.iv.last", "kel.iv.obs", "kel.iv.pred", "kel.last", "kel.obs", "kel.pred", "lambda.z"),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc <- function(concu) {
  if (!missing(concu)) {
    data.frame(
      PPORRESU=concu,
      PPTESTCD=c("cav", "ceoi", "clast.obs", "clast.pred", "cmax", "cmin", "ctrough"),
      stringsAsFactors=FALSE
    )
  }
}

pknca_units_table_amount <- function(amountu) {
  if (!missing(amountu)) {
    data.frame(
      PPORRESU=amountu,
      PPTESTCD="ae",
      stringsAsFactors=FALSE
    )
  }
}

pknca_units_table_conc_dose <- function(concu, doseu) {
  if (!missing(concu) & !missing(doseu)) {
    rbind(
      data.frame(
        PPORRESU=sprintf("(%s)/(%s)", concu, doseu),
        PPTESTCD=c("cav.dn", "clast.obs.dn", "clast.pred.dn", "cmax.dn", "cmin.dn", "ctrough.dn"),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # Volume units
        PPORRESU=sprintf("(%s)/(%s)", doseu, concu),
        PPTESTCD=
          c(
            "vd.obs", "vd.pred",
            "vss.iv.last", "vss.iv.obs", "vss.iv.pred", "vss.last",
            "vss.md.obs", "vss.md.pred", 
            "vss.obs", "vss.pred",
            "vz.obs", "vz.pred"
          ),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc_time <- function(concu, timeu) {
  if (!missing(concu) & !missing(timeu)) {
    rbind(
      data.frame(
        # AUC units
        PPORRESU=sprintf("%s*%s", timeu, concu),
        PPTESTCD=
          c(
            "aucall", "aucinf.obs", "aucinf.pred",
            "aucint.all", "aucint.inf.obs", "aucint.inf.pred", "aucint.last",
            "auclast",
            # These are dose-aware, not dose-normalized
            "aucint.all.dose",  "aucint.inf.obs.dose",
            "aucint.inf.pred.dose", "aucint.last.dose"
          ),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # AUMC units
        PPORRESU=sprintf("%s^2*%s", timeu, concu),
        PPTESTCD=c("aumcall", "aumcinf.obs", "aumcinf.pred", "aumclast"),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc_time_dose <- function(concu, timeu, doseu) {
  if (!missing(concu) & !missing(timeu) & !missing(doseu)) {
    rbind(
      data.frame(
        # AUC units, dose-normalized
        PPORRESU=sprintf("(%s*%s)/(%s)", timeu, concu, doseu),
        PPTESTCD=
          c(
            "aucall.dn", "aucinf.obs.dn", "aucinf.pred.dn",
            "auclast.dn"
          ),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # AUMC units, dose-normalized
        PPORRESU=sprintf("(%s^2*%s)/(%s)", timeu, concu, doseu),
        PPTESTCD=c("aumcall.dn", "aumcinf.obs.dn", "aumcinf.pred.dn", "aumclast.dn"),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # Clearance units
        PPORRESU=sprintf("(%s)/(%s*%s)", doseu, timeu, concu),
        PPTESTCD=c("cl.all", "cl.last", "cl.obs", "cl.pred"),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc_time_amount <- function(concu, timeu, amountu) {
  if (!missing(concu) & !missing(timeu) & !missing(amountu)) {
    data.frame(
      # Renal clearance units
      PPORRESU=sprintf("(%s)/(%s*%s)", amountu, timeu, concu),
      PPTESTCD=c("clr.last", "clr.obs", "clr.pred"),
      stringsAsFactors=FALSE
    )
  }
}

#' Find NCA parameters with a given unit type
#' 
#' @param unit_type The type of unit as assigned with \code{add.interval.col}
#' @return A character vector of parameters with a given unit type
#' @keywords Internal
pknca_find_units_param <- function(unit_type) {
  stopifnot(length(unit_type) == 1)
  stopifnot(is.character(unit_type))
  all_intervals <- get.interval.cols()
  ret <- character()
  for (nm in names(all_intervals)) {
    if (all_intervals[[nm]]$unit_type == unit_type) {
      ret <- c(ret, nm)
    }
  }
  if (length(ret) == 0) {
    stop("No parameters found for unit_type=", unit_type)
  }
  ret
}
