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
})

test_that("pknca_units_add_paren", {
  expect_equal(pknca_units_add_paren("mg"), "mg")
  expect_equal(pknca_units_add_paren("mg/kg"), "(mg/kg)")
  expect_equal(pknca_units_add_paren("mg*kg"), "(mg*kg)")
})
