test_that("units are preserved from input to output", {
  d_conc <-
    tibble::tibble(
      conc=units::set_units(c(8, 4, 2, 1), "ng/mL"),
      time=units::set_units(0:3, "hour")
    )
  o_conc <- PKNCAconc(conc~time, data=d_conc)
  d_intervals <- PKNCA.options("single.dose.aucs")
  d_intervals$start <- units::set_units(d_intervals$start, "hour")
  d_intervals$end <- units::set_units(d_intervals$end, "hour")
  o_data <- PKNCAdata(o_conc, intervals=d_intervals)
  o_nca <- pk.nca(o_data)
})