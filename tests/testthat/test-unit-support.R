test_that("all units have been defined", {
  all_units <- pknca_units_table(concu="conc", doseu="dose", amountu="amount", timeu="time")
  expect_equal(
    setdiff(names(get.interval.cols()), c(all_units$PPTESTCD, "start", "end")),
    character()
  )
})

test_that("pknca_find_units_param expected errors", {
  expect_error(
    pknca_find_units_param(unit_type=1:2)
  )
  expect_error(
    pknca_find_units_param(unit_type=1)
  )
  expect_error(
    pknca_find_units_param(unit_type="foo"),
    regexp="No parameters found for unit_type"
  )
})

test_that("pknca_find_units_param", {
  expect_true(
    "cmax" %in% pknca_find_units_param(unit_type="conc")
  )
})

test_that("unit conversion tables are created correctly", {
  expect_true(
    all(
      pknca_units_table(
        concu="ng/mL", doseu="mg/kg", amountu="mg", timeu="hr"
      )$conversion_factor == 1
    )
  )
  automatic_conversion <-
    pknca_units_table(
      concu="ng/mL", doseu="mg/kg", amountu="mg", timeu="hr",
      conversions=data.frame(
        PPORRESU=c("(mg/kg)/(hr*ng/mL)", "(mg/kg)/(ng/mL)"),
        PPSTRESU=c("mL/hr/kg", "mL/kg")
      )
    )
  expect_equal(
    unique(automatic_conversion$conversion_factor[automatic_conversion$PPSTRESU %in% "mL/hr/kg"]),
    1e6,
    tolerance=.Machine$double.eps
  )
  expect_equal(
    unique(automatic_conversion$conversion_factor[automatic_conversion$PPSTRESU %in% "mL/kg"]),
    1e6,
    tolerance=.Machine$double.eps
  )
  expect_equal(
    unique(automatic_conversion$conversion_factor[!(automatic_conversion$PPSTRESU %in% c("mL/hr/kg", "mL/kg"))]),
    1
  )

  mixed_conversion <-
    pknca_units_table(
      concu="mg/L", doseu="mg/kg", amountu="mg", timeu="hr",
      # Convert clearance and volume units to molar units (assuming
      conversions=data.frame(
        PPORRESU=c("mg/L", "(mg/kg)/(hr*mg/L)", "(mg/kg)/(mg/L)"),
        PPSTRESU=c("mmol/L", "mL/hr/kg", "mL/kg"),
        # Manual conversion of concentration units from ng/mL to mmol/L (assuming
        # a molecular weight of 138.121 g/mol)
        conversion_factor=c(1/138.121, NA, NA)
      )
    )
  expect_equal(
    unique(mixed_conversion$conversion_factor[mixed_conversion$PPSTRESU %in% "mL/hr/kg"]),
    1000,
    tolerance=1e-10
  )
  expect_equal(
    unique(mixed_conversion$conversion_factor[mixed_conversion$PPSTRESU %in% "mL/kg"]),
    1000,
    tolerance=1e-10
  )
  expect_equal(
    unique(mixed_conversion$conversion_factor[mixed_conversion$PPSTRESU %in% "mmol/L"]),
    1/138.121,
    tolerance=1e-10
  )
  expect_equal(
    unique(mixed_conversion$conversion_factor[!(mixed_conversion$PPSTRESU %in% c("mL/hr/kg", "mL/kg", "mmol/L"))]),
    1
  )
})

test_that("pknca_units_table", {
  expect_warning(
      pknca_units_table(
        concu="ng/mL", doseu="mg/kg", amountu="mg", timeu="hr",
        conversions=data.frame(PPORRESU="foo", PPSTRESU="bar", conversion_factor=1)
      ),
      regexp="The following unit conversions were supplied but do not match any units to convert"
  )
  # units library errors occur, if units are not convertible
  expect_error(
    pknca_units_table(
      concu="ng/mL", doseu="mg/kg", amountu="mg", timeu="hr",
      conversions=data.frame(PPORRESU="ng/mL", PPSTRESU="mol/L")
    ),
    regexp="cannot convert ng/mL into mol/L"
  )
  expect_error(
    pknca_units_table(
      concu=c("ng/mL", "umol/L"), doseu="mg/kg", amountu="mg", timeu="hr",
      conversions=data.frame(PPORRESU="ng/mL", PPSTRESU="mol/L")
    ),
    regexp = "Only one unit may be provided at a time: ng/mL, umol/L"
  )
})

test_that("pknca_units_add_paren", {
  expect_equal(pknca_units_add_paren("mg"), "mg")
  expect_equal(pknca_units_add_paren("mg/kg"), "(mg/kg)")
  expect_equal(pknca_units_add_paren("mg*kg"), "(mg*kg)")
})

test_that("pknca_units_table treats missing, NULL, and NA the same", {
  expect_equal(
    pknca_units_table(),
    pknca_units_table(concu = NULL, doseu = NULL, amountu = NULL, timeu = NULL)
  )
  expect_equal(
    pknca_units_table(),
    pknca_units_table(concu = NA, doseu = NULL, amountu = NULL, timeu = NULL)
  )
  expect_equal(
    pknca_units_table(),
    pknca_units_table(concu = NA, doseu = NA, amountu = NULL, timeu = NULL)
  )
  expect_equal(
    pknca_units_table(),
    pknca_units_table(concu = NA, doseu = NA, amountu = NA, timeu = NULL)
  )
  expect_equal(
    pknca_units_table(),
    pknca_units_table(concu = NA, doseu = NULL, amountu = NULL, timeu = NA)
  )
  expect_true(all(is.na(
    pknca_units_table(concu = "ng/mL", doseu = "mg", amountu = "umol", timeu = NULL) %>%
      dplyr::filter(PPTESTCD %in% c("start", "lambda.z", "auclast", "aumclast", "auclast.dn", "aumclast.dn", "cl.last", "clr.last")) %>%
      dplyr::pull("PPORRESU")
  )))
  expect_true(all(is.na(
    pknca_units_table(concu = "ng/mL", doseu = "mg", amountu = NULL, timeu = "hr") %>%
      dplyr::filter(PPTESTCD %in% c("ae", "clr.last")) %>%
      dplyr::pull("PPORRESU")
  )))
  expect_true(all(is.na(
    pknca_units_table(concu = "ng/mL", doseu = NULL, amountu = "umol", timeu = "hr") %>%
      dplyr::filter(PPTESTCD %in% c("cmax.dn", "vss.last", "cl.last", "auclast.dn", "aumclast.dn")) %>%
      dplyr::pull("PPORRESU")
  )))
  expect_true(all(is.na(
    pknca_units_table(concu = NULL, doseu = "mg", amountu = "umol", timeu = "hr") %>%
      dplyr::filter(PPTESTCD %in% c("cmax", "cmax.dn", "vss.last", "auclast", "aumclast", "auclast.dn", "aumclast.dn", "cl.last", "clr.last")) %>%
      dplyr::pull("PPORRESU")
  )))
})

test_that("allow duplicate PPSTRESU units", {
  d_conversion <-
    data.frame(
      PPORRESU = c("ng/mL", "(ng/mL)/(mg/kg)", "(mg/kg)/(hr*ng/mL)", "(mg/kg)/(ng/mL)"),
      PPSTRESU = c("mg/mL", "mL/kg", "mL/(h*kg)", "mL/kg")
    )
  # No error for consistent volume conversion
  expect_silent(
    pknca_units_table(concu = "ng/mL", doseu = "mg/kg", timeu = "hr", conversions = d_conversion)
  )
})

test_that("Use preferred units (#197)", {
  prep <-
    pknca_units_table(
      concu = "ng/mL", doseu = "mg/kg", timeu = "hr", amountu = "mg",
      concu_pref = "ug/mL"
    )
  expect_equal(prep$conversion_factor[prep$PPTESTCD == "cmax"], 0.001)
  prep <-
    pknca_units_table(
      concu = "ng/mL", doseu = "mg/kg", timeu = "hr", amountu = "mg",
      doseu_pref = "ug/kg"
    )
  expect_equal(prep$conversion_factor[prep$PPTESTCD == "cmax.dn"], 0.001)
  prep <-
    pknca_units_table(
      concu = "ng/mL", doseu = "mg/kg", timeu = "hr", amountu = "mg",
      timeu_pref = "day"
    )
  expect_equal(prep$conversion_factor[prep$PPTESTCD == "tmax"], 1/24)
  prep <-
    pknca_units_table(
      concu = "ng/mL", doseu = "mg/kg", timeu = "hr", amountu = "mg",
      amountu_pref = "kg"
    )
  expect_equal(prep$conversion_factor[prep$PPTESTCD == "clr.obs"], 1e-6)

  # conversions can override the preferred units parameter
  prep <-
    pknca_units_table(
      concu = "ng/mL", doseu = "mg/kg", timeu = "hr", amountu = "mg",
      timeu_pref = "day",
      conversions = data.frame(PPORRESU = "hr^2*ng/mL", PPSTRESU = "min^2*ng/mL")
    )
  prep_no_conversions <-
    pknca_units_table(
      concu = "ng/mL", doseu = "mg/kg", timeu = "hr", amountu = "mg",
      timeu_pref = "day"
    )
  expect_equal(prep$conversion_factor[prep$PPTESTCD == "tmax"], 1/24)
  expect_equal(prep$conversion_factor[prep$PPTESTCD == "aumcall"], 3600)
  expect_equal(prep_no_conversions$conversion_factor[prep_no_conversions$PPTESTCD == "aumcall"], 24^-2)

  expect_error(
    pknca_units_table(
      concu = "ng/mL", doseu = "mg/kg", timeu = "hr", amountu = "mg",
      timeu_pref = "day",
      conversions = data.frame(PPORRESU = "A", PPSTRESU = "B")
    ),
    regexp = "Cannot find PPORRESU match between conversions and preferred unit conversions.  Check PPORRESU values in 'conversions' argument.",
    fixed = TRUE
  )

  # Not all arguments are required
  expect_silent(
    pknca_units_table(
      concu = "ng/mL",
      doseu = "mg",
      timeu = "hr",
      timeu_pref = "day"
    )
  )
})

test_that("pknca_units_table expected errors", {
  expect_error(
    pknca_units_table(conversions = "A")
  )
  expect_error(
    pknca_units_table(conversions = data.frame(A = 1)),
    # Generate the error to match (in case its text changes)
    regexp =
      attr(
        try(
          checkmate::assert_names("A", subset.of = c("PPORRESU", "PPSTRESU", "conversion_factor"), .var.name = "names(conversions)"),
          silent = TRUE
        ),
        "condition"
      )$message,
    fixed = TRUE
  )
})

test_that("pknca_unit_conversion", {
  results <- data.frame(PPORRES = 1, PPTESTCD = "cmax")

  # No change when no unit conversion occurs
  expect_equal(pknca_unit_conversion(results, units = NULL), results)

  # Adding units with no conversion
  d_u <- data.frame(PPTESTCD = "cmax", PPORRESU = "ng/mL")
  results_u <- data.frame(PPORRES = 1, PPTESTCD = "cmax", PPORRESU = "ng/mL")
  expect_equal(pknca_unit_conversion(results, units = d_u), results_u)

  # Adding units with conversion
  d_u_conv <- data.frame(PPTESTCD = "cmax", PPORRESU = "ng/mL", conversion_factor = 0.001, PPSTRESU = "ug/mL")
  results_u_conv <- data.frame(PPORRES = 1, PPTESTCD = "cmax", PPORRESU = "ng/mL", PPSTRESU = "ug/mL", PPSTRES = 0.001)
  expect_equal(pknca_unit_conversion(results, units = d_u_conv), results_u_conv)

  # Adding units with some units missing gives an error
  d_u_missing <- data.frame(PPTESTCD = "cmax", PPORRESU = NA_character_)
  results_u_missing <- data.frame(PPORRES = 1, PPTESTCD = "cmax", PPORRESU = NA_character_)
  expect_error(
    pknca_unit_conversion(results, units = d_u_missing),
    regexp = "Units are provided for some but not all parameters; missing for: cmax\nThis error can be converted to a warning using `PKNCA.options(allow_partial_missing_units = TRUE)`",
    fixed = TRUE
  )

  # Adding units with some units missing gives an error; that error can be converted to a warning
  d_u_missing <- data.frame(PPTESTCD = "cmax", PPORRESU = NA_character_)
  results_u_missing <- data.frame(PPORRES = 1, PPTESTCD = "cmax", PPORRESU = NA_character_)
  expect_warning(
    pknca_unit_conversion(results, units = d_u_missing, allow_partial_missing_units = TRUE),
    regexp = "Units are provided for some but not all parameters; missing for: cmax",
    fixed = TRUE
  )

  # Fully-integrated test
  d_conc <- as.data.frame(Theoph[Theoph$Subject %in% Theoph$Subject[1], ])
  o_conc <- PKNCAconc(d_conc, conc~Time|Subject)
  o_dose <- PKNCAdose(d_conc[d_conc$Time == 0, ], Dose~Time|Subject)
  # Do not give dose units
  d_units <- pknca_units_table(concu = "mg/L", timeu = "hr")
  d_interval <- data.frame(start = 0, end = Inf, cmax = TRUE, cl.obs = TRUE)
  o_data <- PKNCAdata(o_conc, o_dose, intervals = d_interval, units = d_units)
  expect_error(
    pk.nca(o_data),
    regexp = "Units are provided for some but not all parameters; missing for: cl.obs\nThis error can be converted to a warning using `PKNCA.options(allow_partial_missing_units = TRUE)`",
    fixed = TRUE
  )
  o_data_warn <- PKNCAdata(o_conc, o_dose, intervals = d_interval, units = d_units, options = list(allow_partial_missing_units = TRUE))
  expect_warning(
    pk.nca(o_data_warn),
    regexp = "Units are provided for some but not all parameters; missing for: cl.obs",
    fixed = TRUE
  )
})
