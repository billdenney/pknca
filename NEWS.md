Note that backward compatibility of function arguments will not be
(near-)guaranteed until version 1.0.  Argument and function changes
will continue until then.  These will be especially noticeable around
the inclusion of IV NCA parameters and additional specifications of
the dosing including dose amount and route.

# PKNCA 0.11.0.9000

## Breaking changes

* PKNCA will now give an error when there are unexpected interval columns.
  The `keep_interval_cols` option can be used to mitigate this error.
* `NA` results from calculating `c0` will now add an exclusion reason.
* AUC for intravenous dosing (all the `auciv*` parameters) now more robustly
  calculate `c0` and does not raise an error when `is.na(c0)` (#353).
* Manual calculation of half.life no longer allows negative half-live values
  (#373).

## New Features

* `PKNCAconc()` and `PKNCAdose()` can now accept unit specifications as either
  column names or units to use (#336).
* PKNCA options can now use `tmax` as a reference for BLQ handling by using new
  names in the `conc.blq` argument (`before.tmax`,`after.tmax`)
* A new parameter `count_conc_measured` was added to enable quality checks,
  typically on AUC measurements. An associated exclusion function,
  `exclude_nca_count_conc_measured()` was also added.
* The `PKNCAconc()` arguments of `include_half.life` and `exclude_half.life` now
  allow `NA` values. If all values are `NA`, then no inclusion or exclusion is
  applied (the interval is treated as-is, like the argument had not been given).
  If some values are `NA` for the interval, those are treated as `FALSE`.

# Minor changes (unlikely to affect PKNCA use)

* PKNCA will now verify the `intervals` data.frame when creating PKNCAdata. The
  checking includes confirming intended column naming and ensuring the correct
  data types.
* PKNCA now contains a `getGroups.PKNCAdata` function to capture grouping columns.
* Duplicate data checks now account for excluded rows.  So, if a row is
  duplicated and all but one of the duplicated rows is excluded, it is not an
  error.  (#298)
* Removed native pipes (`|>`) so that PKNCA will work with older versions of R
  (#304).
* Missing dosing times to `pk.calc.c0()` will not cause an error (#344)

# PKNCA 0.11.0

* PKNCA will now indicate the number of observations included in a summary ("n")
  when it is not the same as the number of subjects included in the summary
  ("N") and the caption will also indicate the definition of "N" and "n".  Note
  that counting of "n" includes all non-missing values that were not excluded
  from summarization; this will included all zeros that are e.g. excluded from
  geometric statistics.
  * If `n == 1`, spread statistics are no longer calculated in the summary.
* A new AUC integration method, "lin-log", has been added using the linear
  method through tmax and log after tmax, with required exceptions for zeros
  (fix #23)
* The parameter `vd` was removed (it was not specific like `vz` or `vss`, and it
  was effectively a duplicate of `vz`).  Use `vz`, instead.
* The `count_conc` NCA parameter was added to assist in data quality checking.
* Subject count ("N") previously counted the number of rows of data, but in
  unusual circumstances, the number of subjects in an NCA result could be fewer
  than the number of rows.  Now the number of subjects is counted (fix #223).
* Extra column in the `intervals` argument to `PKNCAdata()` will no longer cause
  an error (fix #238)
* Many new `assert_*` functions were added to standardize input checking in the
  style of the `checkmate` library.
* Interpolation of zero concentrations in the middle of a set of concentrations
  is now more extensively supported.
* PKNCA has begun the process of deprecating dots in favor of underscores in
  function and parameter names.  Functions with dots instead of underscores
  should continue to work for the foreseeable future (until version 1.0) with
  warnings.
* AUCint will now extrapolate the AUC beyond Tlast using logarithmic
  extrapolation, regardless of the method used (fix #203).
* Imputation will now automatically search for a column named "impute" in the
  interval definition (fix #257).
* Imputation now can look outside the concentration-time of the interval to the
  full concentration-time profile for the group with the `conc.group` and
  `time.group` arguments to the imputation functions.  And,
  `PKNCA_impute_method_start_predose()` imputation performs more reasonably when
  the end of the interval is infinite.
* A progress bar is now available via the `PKNCA.options(progress = )` option
  (fix #193).
* Additional versions of average concentration based on AUCint are now available
  (fix #45).
* A new option "keep_interval_cols" was added to allow keeping a column from the
  intervals in the NCA results.  Note that these are not included in the summary
  groups by default.
* A new argument "filter_requested" for `as.data.frame.PKNCAresults()` allows
  you to filter only to requested results from a PKNCAresults object.
* `pknca_units_table()` now has four new arguments to allow for simplified
  automatic conversion from source units to desired reporting units,
  `concu_pref`, `doseu_pref`, `amountu_pref`, and `timeu_pref` (#197)
* Extraction of PKNCA objects from within other PKNCA objects is now supported
  by various `as_PKNCA*` functions like `as_PKNCAconc()` which can be used to
  extract the concentration data from within a PKNCAdata or PKNCAresults object
  (#278)
* A new "totdose" parameter gives the total dose administered during an interval
* You may exclude parameters from a summary with the new `drop_param` argument
  to `summary()` for PKNCAresults objects.
* The `as.data.frame()` method for `PKNCAresults` objects has a new
  `filter_excluded` argument to remove excluded results from the extracted
  data.frame.  The default behavior is to keep the excluded results with the
  exclude column indicating the reason they were excluded.

## Bugs fixed

* `superpostion()` and the `interp.extrap.conc()` family of functions now
  respect the interpolation and extrapolation types requested rather than using
  default.
* Concentration extrapolation with `extrapolate.conc()` using the "AUCall"
  method now has decreasing instead of increasing concentrations (#249).
* The aucint.inf.obs parameter when calculated with all zero concentrations
  returns zero and aucint.inf.pred returns `NA_real_` (#253)

## Breaking changes

* The arguments `interp.method` and `extrap.method` have been replaced with
  `method` and `auc.type` in the `interp.extrap.conc()` family of functions for
  consistency with the rest of PKNCA (fix #244)
* The AIC.list() function is no longer exported (it was never intended to be an
  external function).
* The `depends` argument to `add.interval.col()` must either be NULL or a
  character vector.
* The names of the `fun.linear`, `fun.log`, and `fun.inf` arguments to
  `pk.calc.auxc` were changed to use underscores.  (If you were using those
  directly, please reach out as they were intended to be internal arguments, and
  I would like to know your use case for changing them.)
* `check.conc.time()` is defunct (it was never intended to be an external
  function).  It has been replaced by `assert_conc()`, `assert_time()` and
  `assert_conc_time()`.
* The clast.obs parameter is now zero when all concentrations are zero (see #253
  for part of the reason).
* (This is not likely to be important for most users.)  The `business...`
  functions (e.g. `business.geomean()`) now include an attribute in non-`NA`
  results with `n`, the number of values included in the statistic.

## Changes under the hood

* Multiple changes were made to speed up calculations.  These will mainly be
  noticed when performing NCA on many subjects (for instance, following
  simulations).  None of these should have external effects that users will
  notice:
  * Adding in dependent parameters required for requested parameters is now more
    efficient (approx 40% time savings)
  * Sorting interval dependencies happens less often (approx 5% time savings)
  * Determining if a parameter is needed for calculation when looking across all
    parameters is more efficient (negligible time savings)
* An internal change was made to make AUC integration and concentration
  interpolation simpler and simplify the ability to create new AUC integration
  or concentration interpolation methods

# PKNCA 0.10.2

* A minor change to `pk.calc.aucpext()` was made so that it now returns
  `NA_real_` instead of `NaN`.
* A minor change was made so that AUC and amount excreted (ae) calculations will
  provide an exclusion reason the result is `NA`.
* A minor new feature makes the specification of imputation easier.  You can
  give the imputation with or without the "PKNCA_impute_method_" prefix.  So,
  "PKNCA_impute_method_start_predose" and "start_predose" are equivalent.

# PKNCA 0.10.1

* A new parameter `aucabove.trough.all` was added to calculate the NCA above the
  trough concentration.
* Testing updates were made to work with dplyr version 1.1.0 (fix #198)
* Internal changes to how columns are identified were made, and the parseFormula
  function was subsequently removed (parseFormula was never intended for
  external use).

# PKNCA 0.10.0

## Bugs Fixed

* When calculating AUC with only a single concentration measurement, NA is now
  returned instead of 0. (fix #176)

## New Features

* Initial support for unit assignment and conversion has been added.  See the
  `units` argument to the `PKNCAdata()` function and the function
  `pknca_units_table()`.
* Initial support for imputation has been added.  See the `impute` argument to
  the `PKNCAdata()` function and the Data Imputation vignette.
* With the addition of units, several outputs now will differ, if units are
  used:
    * `summary()` on a PKNCAresults object shows the units in the column
      heading.
    * When running `as.data.frame()` on a PKNCAresults object with the argument
      `out.format="wide"`, if standardized units values are available, they will
      be used.  And if any unit are available, they will be in the column names.
* Summary tables with units use the "pretty_name" for a parameter which is
  intended for clearer representation in reports.  "pretty_name" use can be
  controlled with the "pretty_names" argument to `summary()`.
    * Note that the pretty names themselves may be modified to help clarify and/or
      shorten the names to make the table heading more useful.  If you intend to
      modify column headers programmatically, set `pretty_names=FALSE` when
      calling the `summary()` function.
* New, IV AUC calculation methods have been added.
* `pk.calc.time_above()` now uses the default AUC calculation method for
  interpolation of time above.  And, it can use 'lin up/log down' interpolation.
* PKNCA can now calculate parameters that require extra information by adding
  the extra information to the intervals data.frame.  For example, add
  `conc_above` as a column to the intervals to allow calculation of
  `time_above`.  With this change, the "conc_above" `PKNCA.options()` value has
  been removed.
* Added dplyr joins, filter, mutate, group_by, and ungroup to allow modification
  of PKNCA objects after creation.  (Note that these functions will make the
  provenance no longer match for PKNCAresults objects.)

## Breaking Changes

* Some functions that were intended to be internal were removed:
    * All `getData()` functions were removed.
    * The `getDataName()` function for PKNCAdata objects was removed.
* `interpolate.conc()` and `interp.extrap.conc()` now give more errors with
  missing (NA) input.  This should not affect typical NCA (where NA values are
  dropped), but it may affect direct calls to the functions themselves.

## Internal Breaking Changes (these should not affect PKNCA users)

* print.parseFormula() was removed from the package.

# PKNCA 0.9.5

* The internals of how PKNCA performs calculations had a significant update. The
  only user-visible change should be that PKNCA does not perform parallel
  computations as of this version. Parallel computation is planned to return in
  the near future.
  * Breaking change:  As part of this change, the split methods for PKNCAconc
    and PKNCAdose objects were removed along with the merge.splitList function.
* Single-subject (ungrouped) analysis works without creating a dummy group (#74)
* PKNCAconc objects are checked earlier for valid data (#154)
* Add time_above parameter to calculate time above a given concentration.
* Fix numeric BLQ replacement when the value is a number and different values
  are given for first, middle, and last (related to #145).  This only affects
  datasets where BLQ is being replaced with a nonzero value (not a common
  scenario).
* Fix issue where intervals could not be tibbles (#141)
* Fix minor issue where only the first exclusion reason would show in the
  exclusion column and other reasons would be ignored (#113). Note that the
  impact of this bug is minimal as the result would have been excluded from
  summaries for the first reason found, but if there were multiple reasons for
  exclusion the subsequent reasons would not be recorded.

# PKNCA 0.9.4

* Additional changes required for compatibility dplyr version 1.0 and CRAN
  checks.  No functionality changed.
* Minor typographical and documentation consistency cleanups throughout.

# PKNCA 0.9.3

* Changes required for compatibility dplyr version 1.0.  No functionality
  changed.

# PKNCA 0.9.2

* New feature: the `time_calc()` function will help convert time values to be
  relative to events (such as calculating time after and before doses)
* Fix issue summarizing results when "start" and "end" are dropped and there are
  multiple interval rows matched for a single group.
* Enable exclusions to be prevented when the input arguments suggest exclusion,
  but the parameter calculating function may be aware of better information about
  exclusion.
* Ensure that exclusions are maintained if an earlier parameter is excluded
  during the initial parameter calculations (Fix #112).
* Two-point half-life calculation works and adjusted r-squared gives a warning
  instead of an error with 2 points (Fix #114).
* Half-life calculation time was decreased by using `.lm.fit()` instead of
  `lm()` decreasing time for a full NCA run by ~30% (and half-life by ~50%).
* For R version 4.0, much more care was taken not to create factors from strings
  unless required (see
  https://developer.r-project.org/Blog/public/2020/02/16/stringsasfactors/index.html)

# PKNCA 0.9.1

* Correct vignette building.

# PKNCA 0.9.0

* Breaking Change: `plot.PKNCAconc()` was moved to the pknca.reporting package
  (https://github.com/billdenney/pknca.reporting)
* Breaking Change: `summary.PKNCAresults()` now provides a caption
  including the summary method for each parameter.  If you change
  summary functions using `PKNCA.set.summary()`, you must now use the
  `description` option to set the description of the summary.
* Breaking Change: ptr now accurately uses ctrough instead of cmin (fix #106)
* Issue fixed where aucint* calculations now respect BLQ and NA rules like other
  calculations. (#104)
* When half.life is not calculated due to insufficient number of points (default
  < 3), an exclusion reason is added. (#102)
* tibbles now work as the interval argument for `PKNCAdata()` (fix #72)
* Issue fixed with summarization of data that has exclusions.
  Exclusions are now correctly handled as missing instead of never
  calculated.
* parseFormula now internally uses NULL for no-group formula definitions.
* signifString and roundString now have sci_range (deprecating si_range) and
  sci_sep arguments.
* Documentation is improved, especially around the selection of
  parameters for intervals.
* Multiple dose data with a single concentration measurement no longer
  generates an error (fixes #84).
* The "start" and "end" columns may now be dropped from the summary of
  `PKNCAresults` objects.
* `PKNCAdata()` is more restrictive on unknown arguments issuing an error
  when unknown arguments are present.
* `intervals` argument to `PKNCAdata()` may now be a tibble (fixes #72).
* Documentation has been extensively updated (fixes #81).
* CRAN changes: Vignettes now better respect not loading suggested
  packages.  Tests are now more permissive in timing.

# PKNCA 0.8.5

* Cleaned AUCint names
* Added dose-count within interval (to warn of multiple doses within an
  interval)
* Various documentation updates
* signifString and roundString now by default use scientific notation
  for values >=1e6 and <=1e-6
* Fix bug in option handling within `pk.nca` (Fix #68)

# PKNCA 0.8.4

* Added AUCint flavors
* Parameter names for NCA parameters will likely be changing in the
  next version; code will still work, but some calculation methods and
  therefore results may be subtly different.  These changes will be
  fully documented.)

# PKNCA 0.8.2

* BACKWARD INCOMPATIBILITY: The function supplied to the exclude
  argument 'FUN' now requires two arguments and operates on the level
  of a single group rather than the full object.  The function can
  also return the reason as a character string instead of a logical
  mask of when to exclude data.
* BACKWARD INCOMPATIBILITY: Added back-end functionality to only
  require one function to handle many NCA parameters that are related
  (e.g. combine pk.calc.aucpext, pk.calc.aucpext.obs,
  pk.calc.aucpext.pred, etc.).  If your current code calls a specific
  function (like pk.calc.aucpext.pred), you must change to using the
  generic function (like pk.calc.aucpext)
* BACKWARD INCOMPATIBILITY: Functions that previously may have
  returned Infinity due to dividing by zero (e.g. when AUC=0
  calculating clearance) now return NA.

* Added Validation vignette.

* Corrected issue where time to steady-state with a single estimate
  may have given more than one estimated time to steady-state.
* Corrected issue with exclude handling where now a blank string is
  also accepted as included (not excluded).
* PKNCAconc now accepts a "volume" argument and pk.nca can now
  calculate urine/feces-related parameters (fe, ae, clr)
* exclude_nca* functions added (Fixes issue #20)
* Add manual half-life point selection (Fixes issue #18)
* Improved summary settings (Fixes issue #54)
* Add parameters for Ceoi and intravenous MRT
* Updated vignettes to improve clarity
* Added dose-normalized PK parameters (Fixes issue #41)
* Added checks to confirm that concentration and time are numeric
  (Fixes feature request #40)
* Improved test coverage

# PKNCA 0.8.1

* A PKNCAdose is no longer required for calculations.
* Data may now be excluded from calculations.

# PKNCA 0.8

This release is not backward compatible.  The switch to observed and
predicted-related NCA parameters (like aucinf.obs and aucinf.pred)
changed the format of the intervals specification.

* Remove dependency on doBy library
* Dose-aware interpolation and extrapolation was added with the interp.extrap.conc.dose function.
* Added Clast.pred related NCA calculations
* Added N to summary of PKNCAresults
* Added parameter selection between Clast,observed and Clast,predicted across all parameters
* Enabled PKNCAdose to be specified with one-sided formula
* Improved error reporting so that the group and time (interval specification) is reported in addition to the error.
* PKNCAdose now allows route of administration and IV infusion parameters of rate/duration to be specified

# PKNCA 0.7.1

* Updated vignettes
* Standardize rounding and significance with missing values in signifString and roundString
* Enable wide data output with as.data.frame(PKNCAresults, out.format="wide")
* Correct calculation of Vz
* Various CRAN-related cleanups

# PKNCA 0.7

* Features added
  * Additional PK parameters to support IV dosing added
  * Fix #11, Intervals can be specified manually, and will apply across appropriate parts of the grouping variables
  * Enable dose and dose.time as parameters to NCA calculations
  * More NCA parameters are calculated, especially related to IV dosing
  * Fix #8, Reporting times for time-based parameters are now within the current interval rather than since first dose (e.g. Tmax on day 14 should be between 0 and 24 not 2*7*24+c(0, 24))
  * Added several vignettes
* Bugs fixed
  * Dosing without concentration is probably placebo; warn and continue
  * Fix #6, make merge.splitByData work with more than one dosing level
  * Export some generic classes that were not previously exported to simplify their use
  * Superposition extensions when lambda.z cannot be calculated
  * Significance rounding into character strings works when the rounding moves up one order of magnitude.
  * Fix #9, summarization of parameters that are not calculated show not calculated instead of missing.

# PKNCA 0.6

First release targeting CRAN
