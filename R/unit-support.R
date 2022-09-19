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
#' pknca_units_table() # only parameters that are unitless
#' pknca_units_table(
#'   concu="ng/mL", doseu="mg/kg", amountu="mg", timeu="hr"
#' )
#' pknca_units_table(
#'   concu="ng/mL", doseu="mg/kg", amountu="mg", timeu="hr",
#'   # Convert clearance and volume units to more understandable units with
#'   # automatic unit conversion
#'   conversions=data.frame(
#'     PPORRESU=c("(mg/kg)/(hr*ng/mL)", "(mg/kg)/(ng/mL)"),
#'     PPSTRESU=c("mL/hr/kg", "mL/kg")
#'   )
#' )
#' pknca_units_table(
#'   concu="mg/L", doseu="mg/kg", amountu="mg", timeu="hr",
#'   # Convert clearance and volume units to molar units (assuming
#'   conversions=data.frame(
#'     PPORRESU=c("mg/L", "(mg/kg)/(hr*ng/mL)", "(mg/kg)/(ng/mL)"),
#'     PPSTRESU=c("mmol/L", "mL/hr/kg", "mL/kg"),
#'     # Manual conversion of concentration units from ng/mL to mmol/L (assuming
#'     # a molecular weight of 138.121 g/mol)
#'     conversion_factor=c(1/138.121, NA, NA)
#'   )
#' )
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
    stop("Please report a bug.  Unknown NCA parameters have units defined: ", paste(extra_cols, collapse=", ")) # nocov
  }
  if (nrow(conversions) > 0) {
    stopifnot(!duplicated(conversions$PPORRESU))
    stopifnot(!duplicated(conversions$PPSTRESU))
    stopifnot(length(setdiff(names(conversions), c("PPORRESU", "PPSTRESU", "conversion_factor"))) == 0)
    if (!("conversion_factor" %in% names(conversions))) {
      if (!requireNamespace("units", quietly=TRUE)) {
        stop("The units package is required for automatic unit conversion") # nocov
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
    unexpected_conversions <- setdiff(conversions$PPORRESU, ret$PPORRESU)
    if (length(unexpected_conversions) > 0) {
      warning(
        "The following unit conversions were supplied but do not match any units to convert: ",
        paste0("'", unexpected_conversions, "'", collapse=", ")
      )
    }
    ret <-
      dplyr::left_join(
        ret, conversions,
        by=intersect(names(ret), names(conversions))
      )
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
      PPTESTCD=pknca_find_units_param(unit_type="unitless"),
      stringsAsFactors=FALSE
    ),
    data.frame(
      PPORRESU="fraction",
      PPTESTCD=pknca_find_units_param(unit_type="fraction"),
      stringsAsFactors=FALSE
    ),
    data.frame(
      PPORRESU="%",
      PPTESTCD=pknca_find_units_param(unit_type="%"),
      stringsAsFactors=FALSE
    ),
    data.frame(
      PPORRESU="count",
      PPTESTCD=pknca_find_units_param(unit_type="count"),
      stringsAsFactors=FALSE
    )
  )
}

pknca_units_table_time <- function(timeu) {
  if (!missing(timeu) && !is.null(timeu)) {
    rbind(
      data.frame(
        PPORRESU=timeu,
        PPTESTCD=pknca_find_units_param(unit_type="time"),
        stringsAsFactors=FALSE
      ),
      data.frame(
        PPORRESU=sprintf("1/%s", pknca_units_add_paren(timeu)),
        PPTESTCD=pknca_find_units_param(unit_type="inverse_time"),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc <- function(concu) {
  if (!missing(concu) && !is.null(concu)) {
    data.frame(
      PPORRESU=concu,
      PPTESTCD=pknca_find_units_param(unit_type="conc"),
      stringsAsFactors=FALSE
    )
  }
}

pknca_units_table_amount <- function(amountu) {
  if (!missing(amountu) && !is.null(amountu)) {
    data.frame(
      PPORRESU=amountu,
      PPTESTCD=pknca_find_units_param(unit_type="amount"),
      stringsAsFactors=FALSE
    )
  }
}

pknca_units_table_conc_dose <- function(concu, doseu) {
  if (!missing(concu) && !missing(doseu) && !is.null(concu) && !is.null(doseu)) {
    rbind(
      data.frame(
        PPORRESU=sprintf("%s/%s", pknca_units_add_paren(concu), pknca_units_add_paren(doseu)),
        PPTESTCD=pknca_find_units_param(unit_type="conc_dosenorm"),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # Volume units
        PPORRESU=sprintf("%s/%s", pknca_units_add_paren(doseu), pknca_units_add_paren(concu)),
        PPTESTCD=pknca_find_units_param(unit_type="volume"),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc_time <- function(concu, timeu) {
  if (!missing(concu) && !missing(timeu) && !is.null(concu) && !is.null(timeu)) {
    rbind(
      data.frame(
        # AUC units
        PPORRESU=sprintf("%s*%s", timeu, concu),
        PPTESTCD=pknca_find_units_param(unit_type="auc"),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # AUMC units
        PPORRESU=sprintf("%s^2*%s", pknca_units_add_paren(timeu), concu),
        PPTESTCD=pknca_find_units_param(unit_type="aumc"),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc_time_dose <- function(concu, timeu, doseu) {
  if (!missing(concu) && !missing(timeu) && !missing(doseu) &&
      !is.null(concu) && !is.null(timeu) && !is.null(doseu)) {
    rbind(
      data.frame(
        # AUC units, dose-normalized
        PPORRESU=sprintf("(%s*%s)/%s", timeu, concu, pknca_units_add_paren(doseu)),
        PPTESTCD=pknca_find_units_param(unit_type="auc_dosenorm"),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # AUMC units, dose-normalized
        PPORRESU=sprintf("(%s^2*%s)/%s", pknca_units_add_paren(timeu), concu, pknca_units_add_paren(doseu)),
        PPTESTCD=pknca_find_units_param(unit_type="aumc_dosenorm"),
        stringsAsFactors=FALSE
      ),
      data.frame(
        # Clearance units
        PPORRESU=sprintf("%s/(%s*%s)", pknca_units_add_paren(doseu), timeu, concu),
        PPTESTCD=pknca_find_units_param(unit_type="clearance"),
        stringsAsFactors=FALSE
      )
    )
  }
}

pknca_units_table_conc_time_amount <- function(concu, timeu, amountu) {
  if (!missing(concu) && !missing(timeu) && !missing(amountu) &&
      !is.null(concu) && !is.null(timeu) && !is.null(amountu)) {
    data.frame(
      # Renal clearance units
      PPORRESU=sprintf("%s/(%s*%s)", pknca_units_add_paren(amountu), timeu, concu),
      PPTESTCD=pknca_find_units_param(unit_type="renal_clearance"),
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

#' Add parentheses to a unit value, if needed
#'
#' @param unit The text of the unit
#' @return The unit with parentheses around it, if needed
#' @keywords Internal
pknca_units_add_paren <- function(unit) {
  mask_paren <- grepl(x=unit, pattern="[*/]")
  ifelse(mask_paren, yes=paste0("(", unit, ")"), no=unit)
}

#' Perform unit conversion (if possible) on PKNCA results
#'
#' @param result The results data.frame
#' @param units The unit conversion table
#' @return The result table with units converted
#' @keywords Internal
pknca_unit_conversion <- function(result, units) {
  ret <- result
  if (!is.null(units)) {
    ret <-
      dplyr::left_join(
        ret, units,
        by=intersect(names(ret), names(units))
      )
    if ("conversion_factor" %in% names(units)) {
      ret$PPSTRES <- ret$PPORRES * ret$conversion_factor
      # Drop the conversion factor column, since it shouldn't be in the output.
      ret <- ret[, setdiff(names(ret), "conversion_factor")]
    }
  }
  ret
}
