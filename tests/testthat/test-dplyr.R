# based on purrr:::capture_output
msg_grabber <- function(code) {
  messages <- character()
  handler <- function(m) {
    messages <<- c(messages, m$message)
    invokeRestart("muffleMessage")
  }
  withCallingHandlers(code, message = handler)
  messages
}

test_that("dplyr filter", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  filtered <- filter(myresult, PPTESTCD == "auclast")
  filtered_manual <- myresult
  filtered_manual$result <- filtered_manual$result[filtered_manual$result$PPTESTCD == "auclast", ]
  expect_equal(filtered, filtered_manual)

  filtered <- filter(myconc, ID == 1)
  filtered_manual <- myconc
  filtered_manual$data <- myconc$data[myconc$data$ID == 1, ]
  expect_equal(filtered, filtered_manual)
})

test_that("dplyr left_join", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  joindf <- data.frame(ID=1, foo="bar")
  msg_join_id <- msg_grabber(left_join(data.frame(ID = 1), data.frame(ID = 1)))
  expect_message(
    joined <- left_join(myresult, joindf),
    msg_join_id,
    fixed = TRUE
  )
  joined_manual <- myresult
  expect_message(
    joined_manual$result <- dplyr::left_join(joined_manual$result, joindf),
    msg_join_id,
    fixed = TRUE
  )
  expect_equal(joined, joined_manual)

  expect_message(
    joined <- left_join(myconc, joindf),
    msg_join_id,
    fixed = TRUE
  )
  joined_manual <- myconc
  expect_message(
    joined_manual$data <- left_join(joined_manual$data, joindf),
    msg_join_id,
    fixed = TRUE
  )
  expect_equal(joined, joined_manual)
})

test_that("dplyr mutate", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  mutated <- mutate(myresult, foo="bar")
  mutated_manual <- myresult
  mutated_manual$result <- mutate(mutated_manual$result, foo="bar")
  expect_equal(mutated, mutated_manual)

  mutated <- mutate(myconc, foo="bar")
  mutated_manual <- myconc
  mutated_manual$data <- mutate(mutated_manual$data, foo="bar")
  expect_equal(mutated, mutated_manual)
})

test_that("dplyr group_by and ungroup", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  grouped <- group_by(myconc, treatment)
  expect_s3_class(grouped$data, "grouped_df")
  ungrouped <- ungroup(grouped)
  expect_false("grouped_df" %in% class(ungrouped$data))
  expect_s3_class(ungrouped$data, "data.frame")
})
