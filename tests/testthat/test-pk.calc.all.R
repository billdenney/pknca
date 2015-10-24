context("All NCA calculations")

generate.data <- function(nsub, ntreat, time.points, resid=0.1) {
  one.cmt.oral <- function(time, dose, v, ka, kel, f=1)
    f*dose*ka/(v*(ka-kel))*(exp(-kel*time) - exp(-ka*time))
  ## one.cmt.iv <- function(time, dose, v, kel)
  ##   dose/v*exp(-kel*time)
  set.seed(5)
  ret <- data.frame()
  for (i in 1:ntreat)
    for (j in 1:nsub)
      ret <- rbind(ret,
                   data.frame(treatment=paste("Trt", i),
                              ID=j,
                              time=time.points,
                              conc=exp(
                                rnorm(length(time.points), mean=0, sd=resid))*
                                one.cmt.oral(time.points,
                                             dose=i,
                                             v=1,
                                             ka=1,
                                             kel=0.05),
                              stringsAsFactors=FALSE))
  ret
}

test_that("pk.nca.interval", {
  tmpdata <- generate.data(nsub=5, ntreat=2,
                           time.points=c(0, 1, 2, 4, 6, 8, 12, 24))
})
