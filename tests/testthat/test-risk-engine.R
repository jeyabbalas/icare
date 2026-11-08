engine_setup <- function() {
  baseline <- list(age = 20:89, hazard = seq(2e-4, 4e-3, length.out = 70))
  competing <- data.frame(age = 20:89, rate = seq(1e-3, 6e-2, length.out = 70))
  list(baseline = baseline, competing = competing)
}

test_that("the vectorized engine equals the defining sum for whole and fractional ages", {
  setup <- engine_setup()
  set.seed(9)
  n <- 300
  lp <- stats::rnorm(n, sd = 0.7)
  start <- c(sample(20:60, n / 2, TRUE), stats::runif(n / 2, 20, 60))
  length <- c(sample(0:25, n / 2, TRUE), stats::runif(n / 2, 0, 25))
  result <- compute_absolute_risk_engine(lp, start, length, setup$baseline,
                                         setup$competing, survival = TRUE)
  expected <- t(vapply(seq_len(n), function(i) {
    reference_absolute_risk(lp[i], start[i], length[i], setup$baseline, setup$competing)
  }, numeric(2)))
  expect_equal(result$risk, expected[, "risk"], tolerance = 1e-12)
  expect_equal(result$survival, expected[, "survival"], tolerance = 1e-12)
  expect_equal(result$risk[length == 0], rep(0, sum(length == 0)))
})

test_that("risks over adjacent windows chain through the survival", {
  setup <- engine_setup()
  lp <- c(-0.5, 0, 0.8)
  whole <- compute_absolute_risk_engine(lp, 40.5, 20.5, setup$baseline, setup$competing)
  # exactly at a whole-age cut
  first <- compute_absolute_risk_engine(lp, 40.5, 6.5, setup$baseline, setup$competing,
                                        survival = TRUE)
  second <- compute_absolute_risk_engine(lp, 47, 14, setup$baseline, setup$competing)
  expect_equal(whole, first$risk + first$survival * second, tolerance = 1e-14)
  # a cut inside a year of age is not exactly additive in the yearly model
  first <- compute_absolute_risk_engine(lp, 40.5, 6.75, setup$baseline, setup$competing,
                                        survival = TRUE)
  second <- compute_absolute_risk_engine(lp, 47.25, 13.75, setup$baseline,
                                         setup$competing)
  chained <- first$risk + first$survival * second
  expect_false(isTRUE(all.equal(whole, chained, tolerance = 1e-8)))
  expect_equal(whole, chained, tolerance = 1e-3)
})

test_that("risk increases with the relative risk and the window, within [0, 1]", {
  setup <- engine_setup()
  risk <- compute_absolute_risk_engine(c(-1, 0, 1, 3), 50, 30, setup$baseline,
                                       setup$competing)
  expect_true(all(diff(risk) > 0))
  expect_true(all(risk > 0 & risk < 1))
  by_length <- compute_absolute_risk_engine(rep(0, 4), 50, c(1, 5, 10, 30),
                                            setup$baseline, setup$competing)
  expect_true(all(diff(by_length) > 0))
})
