test_that("the calibrated baseline hazard reproduces the incidence rates", {
  set.seed(7)
  rates <- data.frame(age = 20:89, rate = seq(1e-4, 4e-3, length.out = 70))
  lp <- stats::rnorm(500, sd = 0.6)
  weights <- stats::runif(500)
  baseline <- calibrate_baseline_hazard(rates, lp, weights / sum(weights),
                                        tolerance = 1e-12)
  # the survivor-weighted mean hazard at each age equals the incidence rate
  rr <- exp(lp)
  cumulative <- c(0, cumsum(baseline$hazard))[seq_len(70)]
  implied <- vapply(seq_len(70), function(k) {
    s <- exp(-rr * cumulative[k])
    baseline$hazard[k] * sum(weights * s * rr) / sum(weights * s)
  }, numeric(1))
  expect_equal(implied, rates$rate, tolerance = 1e-10)
  expect_equal(baseline$hazard, reference_baseline_hazard(rates, lp, weights),
               tolerance = 1e-10)
})

test_that("the default tolerance stops where earlier versions stopped", {
  rates <- data.frame(age = 0:85, rate = seq(0, 5e-3, length.out = 86))
  lp <- seq(-1, 1, length.out = 101)
  baseline <- calibrate_baseline_hazard(rates, lp, rep(1 / 101, 101))
  expect_gte(baseline$iterations, 2L)
  expect_equal(baseline$hazard, reference_baseline_hazard(rates, lp, rep(1, 101)),
               tolerance = 1e-3)
})

test_that("integer weights are equivalent to replicated reference rows", {
  set.seed(8)
  rates <- data.frame(age = 30:80, rate = seq(5e-4, 3e-3, length.out = 51))
  lp <- stats::rnorm(60, sd = 0.8)
  counts <- sample(1:5, 60, replace = TRUE)
  weighted <- calibrate_baseline_hazard(rates, lp, counts / sum(counts))
  replicated <- calibrate_baseline_hazard(rates, rep(lp, counts),
                                          rep(1, sum(counts)) / sum(counts))
  expect_equal(weighted$hazard, replicated$hazard, tolerance = 1e-12)
  expect_identical(weighted$iterations, replicated$iterations)
})

test_that("a population without risk factors has the incidence rates as baseline", {
  rates <- data.frame(age = 40:60, rate = seq(1e-3, 3e-3, length.out = 21))
  baseline <- calibrate_baseline_hazard(rates, rep(0, 10), rep(0.1, 10))
  expect_equal(baseline$hazard, rates$rate)
  expect_error(calibrate_baseline_hazard(rates, c(0, 1000), c(0.5, 0.5)), "not finite")
})
