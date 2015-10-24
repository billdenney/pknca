context("Superposition")

## The requirements for superposition are defined below to enable
## test-driven development.

#' Dosing must be single dose
#' * A single dose within the time period
#'   * Time period may be a subset of the data from a multi-dose study
#' * Concentrations must be BLQ at the beginning of the interval
#'   * User can over-ride the requirement of BLQ at the beginning of the interval
#' * Half-life
#'   * If the number of dosing intervals times the duration of the dosing interval is less than the observed data, half-life is not required.  (This will likely be rare.)
#'   * Half-life can be specified on input.
#'   * Half-life can be sufficiently determined from the data given (see the half-life cleaning functions for the definition of "sufficient")
#'   * Mixing methods of half-life determination are supported.
#' * Extrapolation can occur via either Clast,obs or Clast,pred.
#'   * If using Clast,pred, it must be given if the half-life is given.
#' * Dosing interval can be specified as either a simple number for dosing interval (tau) or multiple doses within the interval
#'   * Example: Daily dosing would be: 24
#'   * Example: Breakfast/dinner dosing would be: c(10, 24)
#' * How many times to repeat the dosing interval can either be given as a number of taus to repeat or "steady-state"
#'   * If steady-state, then a relative tolerance will be provided to confirm that no concentration changes by a greater ratio from the previous to the current dose.
#' * Doses amount may be provided
#'   * If doses are not given, they are assumed to be the same as the input dose.
#'   * If doses are given, linear dose-proportionality will be assumed.
#'   * Dose amount can be given either as a scalar which will apply to all intervals or a vector which will apply to each interval.
#'     * If dose amount is given as a vector, it must be the same length as the number of intervals input.
#' * User can specify either integration as either linear or linear-up/log-down
#' * Time points in the output are automatically selected as all observed points modulo the dosing interval.
#'   * More time points can be added by user specification, but points cannot be removed.  (Note that this will ensure the highest possible Cmax in the output data.)
#' * The functions work with either individual or grouped data.
