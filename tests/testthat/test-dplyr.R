source("generate.data.R")

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
  joined <- left_join(myresult, joindf)
  joined_manual <- myresult
  joined_manual$result <- left_join(joined_manual$result, joindf)
  expect_equal(joined, joined_manual)
  
  joined <- left_join(myconc, joindf)
  joined_manual <- joined
  joined_manual$data <- left_join(joined_manual$data, joindf)
  expect_equal(joined, joined_manual)
})
