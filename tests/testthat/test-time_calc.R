context("time_calc()")

test_that("Numeric time_calc() is successful with all scalar combinations of number and NA and numeric order.", {
  expect_equal(
    time_calc(time_event=0, time_obs=0),
    data.frame(
      event_number_before=1,
      event_number_after=1,
      time_after_event=0,
      time_before_event=0,
      time_after_first=0
    )
  )
  expect_equal(
    time_calc(time_event=1, time_obs=0),
    data.frame(
      event_number_before=NA_integer_,
      event_number_after=1,
      time_after_event=NA_real_,
      time_before_event=-1,
      time_after_first=-1
    )
  )
  expect_equal(
    time_calc(time_event=0, time_obs=1),
    data.frame(
      event_number_before=1,
      event_number_after=NA_integer_,
      time_after_event=1,
      time_before_event=NA_real_,
      time_after_first=1
    )
  )
  expect_equal(
    time_calc(time_event=NA_real_, time_obs=1),
    data.frame(
      event_number_before=NA_integer_,
      event_number_after=NA_integer_,
      time_after_event=NA_real_,
      time_before_event=NA_real_,
      time_after_first=NA_real_
    )
  )
  expect_equal(
    time_calc(time_event=0, time_obs=NA_real_),
    data.frame(
      event_number_before=NA_integer_,
      event_number_after=NA_integer_,
      time_after_event=NA_real_,
      time_before_event=NA_real_,
      time_after_first=NA_real_
    )
  )
  expect_equal(
    time_calc(time_event=NA_real_, time_obs=NA_real_),
    data.frame(
      event_number_before=NA_integer_,
      event_number_after=NA_integer_,
      time_after_event=NA_real_,
      time_before_event=NA_real_,
      time_after_first=NA_real_
    )
  )
  expect_error(
    time_calc(
      time_event=5,
      time_obs=as.difftime(5, units="hours")
    ),
    regexp="Both `time_event` and `time_obs` must be the same class (numeric).",
    fixed=TRUE
  )
  expect_equal(
    expect_warning(
      time_calc(time_event=numeric(0), time_obs=1),
      regexp="No events provided",
      fixed=TRUE
    ),
    time_calc(time_event=NA_real_, time_obs=1)
  )
})

test_that("Numeric time_calc() is successful with all vector combinations of number and NA and numeric order.", {
  expect_equal(
    time_calc(time_event=c(0, 10, 20), time_obs=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)),
    data.frame(
      event_number_before=c(NA_integer_, 1, 1, 1, 2, 2, 2, 3, 3),
      event_number_after=c(1, 1, 2, 2, 2, 3, 3, 3, NA_integer_),
      time_after_event=c(NA_real_, 0, 1, 9, 0, 1, 9, 0, 1),
      time_before_event=c(-1, 0, -9, -1, 0, -9, -1, 0, NA_real_),
      time_after_first=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)
    )
  )
  expect_equal(
    time_calc(time_event=c(0, 10, 10, 20), time_obs=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)),
    data.frame(
      event_number_before=c(NA_integer_, 1, 1, 1, 3, 3, 3, 4, 4),
      event_number_after=c(1, 1, 2, 2, 2, 4, 4, 4, NA_integer_),
      time_after_event=c(NA_real_, 0, 1, 9, 0, 1, 9, 0, 1),
      time_before_event=c(-1, 0, -9, -1, 0, -9, -1, 0, NA_real_),
      time_after_first=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)
    ),
    info="Duplicated events are counted as the first for event_number_before and the last for event_number_after."
  )
  expect_equal(
    time_calc(time_event=c(0, 10, 10, 20, NA_real_), time_obs=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)),
    data.frame(
      event_number_before=c(NA_integer_, 1, 1, 1, 3, 3, 3, 4, 4),
      event_number_after=c(1, 1, 2, 2, 2, 4, 4, 4, NA_integer_),
      time_after_event=c(NA_real_, 0, 1, 9, 0, 1, 9, 0, 1),
      time_before_event=c(-1, 0, -9, -1, 0, -9, -1, 0, NA_real_),
      time_after_first=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)
    ),
    info="NA time_event is ignored."
  )
  expect_equal(
    time_calc(time_event=c(0, 10, NA_real_, 10, 20), time_obs=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)),
    data.frame(
      event_number_before=c(NA_integer_, 1, 1, 1, 4, 4, 4, 5, 5),
      event_number_after=c(1, 1, 2, 2, 2, 5, 5, 5, NA_integer_),
      time_after_event=c(NA_real_, 0, 1, 9, 0, 1, 9, 0, 1),
      time_before_event=c(-1, 0, -9, -1, 0, -9, -1, 0, NA_real_),
      time_after_first=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)
    ),
    info="NA time_event is counted correctly in the middle."
  )
  expect_error(
    time_calc(time_event=c(10, 0), time_obs=c(-1, 0, 1, 9, 10, 11, 19, 20, 21)),
    regexp="`time_event` must be sorted."
  )
})

test_that("POSIXt objects work", {
  expect_equal(
    time_calc(
      time_event=as.POSIXlt("2019-07-30 06:16"),
      time_obs=as.POSIXct("2019-07-30 06:00"),
      units="hours"
    ),
    data.frame(
      event_number_before=NA_integer_,
      event_number_after=1L,
      time_after_event=NA_real_,
      time_before_event=-16/60,
      time_after_first=-16/60
    )
  )
  expect_equal(
    time_calc(
      time_event=as.POSIXlt("2019-07-30 06:16"),
      time_obs=as.POSIXct("2019-07-30 06:00"),
      units="days"
    ),
    data.frame(
      event_number_before=NA_integer_,
      event_number_after=1L,
      time_after_event=NA_real_,
      time_before_event=-16/1440,
      time_after_first=-16/1440
    ),
    info="Units are respected"
  )
  expect_equal(
    time_calc(
      time_event=as.POSIXlt("2019-07-30 06:16"),
      time_obs=as.POSIXct("2019-07-30 06:00"),
      units="hours"
    ),
    data.frame(
      event_number_before=NA_integer_,
      event_number_after=1L,
      time_after_event=NA_real_,
      time_before_event=-16/60,
      time_after_first=-16/60
    ),
    info="Mixed POSIXlt and POSIXct inputs work"
  )
  expect_error(
    time_calc(
      time_event=as.POSIXlt("2019-07-30 06:16"),
      time_obs=5,
      units="hours"
    ),
    regexp="Both `time_event` and `time_obs` must be the same class (POSIXt).",
    fixed=TRUE
  )
  expect_error(
    time_calc(
      time_event=as.POSIXlt("2019-07-30 06:16"),
      time_obs=as.POSIXct("2019-07-30 06:00")
    ),
    regexp="`units` must be provided.",
    fixed=TRUE
  )
})

test_that("difftime", {
  expect_error(
    time_calc(
      time_event=as.difftime(5, units="hours"),
      time_obs=as.difftime(5, units="hours")
    ),
    regexp="`units` must be provided.",
    fixed=TRUE
  )
  expect_error(
    time_calc(
      time_event=as.difftime(5, units="hours"),
      time_obs=5,
      units="hours"
    ),
    regexp="Both `time_event` and `time_obs` must be the same class (difftime).",
    fixed=TRUE
  )
  expect_equal(
    time_calc(
      time_event=as.difftime(5, units="hours"),
      time_obs=as.difftime(5, units="hours"),
      units="hours"
    ),
    time_calc(
      time_event=5,
      time_obs=5
    )
  )
  expect_equal(
    time_calc(
      time_event=as.difftime(5, units="hours"),
      time_obs=as.difftime(5, units="hours"),
      units="days"
    ),
    time_calc(
      time_event=5/24,
      time_obs=5/24
    ),
    info="Units are not inherited from the difftime object."
    # because they could differ, and the default automatic selection may not
    # match the user's desire.
  )
})
