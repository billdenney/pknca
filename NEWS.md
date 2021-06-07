Note that backward compatibility of function arguments will not be
(near-)guaranteed until version 1.0.  Argument and function changes
will continue until then.  These will be especially noticeable around
the inclusion of IV NCA parameters and additional specifications of
the dosing including dose amount and route.

# PKNCA 0.9.5 (in development)

* Add time_above parameter to calcualte time above a given concentration.
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
