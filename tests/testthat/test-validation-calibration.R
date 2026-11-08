calibration_data <- function(n = 3000, seed = 21) {
  set.seed(seed)
  lp <- stats::rnorm(n)
  predicted <- stats::plogis(lp - 3)
  outcome <- stats::rbinom(n, 1, pmin(1, 1.1 * predicted))
  list(lp = lp, predicted = predicted, outcome = outcome)
}

test_that("the overall E/O interval uses the binomial variance in a cohort", {
  d <- calibration_data()
  result <- overall_expected_observed(d$predicted, d$outcome)
  e <- mean(d$predicted)
  o <- mean(d$outcome)
  n <- length(d$outcome)
  expect_equal(c(result$expected, result$observed, result$ratio), c(e, o, e / o))
  expect_equal(result$ci, exp(log(e / o) + c(-1, 1) * 1.96 * sqrt((1 - o) / (n * o))))
})

test_that("in sampled data the E/O interval adds the sampling variance of E and O", {
  d <- calibration_data()
  set.seed(22)
  prob <- ifelse(d$outcome == 1, 1, stats::runif(length(d$outcome), 0.1, 0.6))
  keep <- stats::runif(length(prob)) < prob
  p <- d$predicted[keep]
  y <- d$outcome[keep]
  w <- 1 / prob[keep]
  result <- overall_expected_observed(p, y, w)
  size <- sum(w)
  e <- sum(w * p) / size
  o <- sum(w * y) / size
  u <- (p / e - y / o) / size
  variance <- (1 - o) / (size * o) + sum(w * (w - 1) * u^2)
  expect_equal(result$ratio, e / o)
  expect_equal(result$ci, exp(log(e / o) + c(-1, 1) * 1.96 * sqrt(variance)))
  expect_equal(overall_expected_observed(p, y, rep(1, length(y))),
               overall_expected_observed(p, y))
})

test_that("group statistics and both Hosmer-Lemeshow variances follow their formulas", {
  d <- calibration_data()
  groups <- risk_score_groups(d$lp, NULL, 5)
  g <- as.integer(groups)
  size <- tabulate(g)
  e <- as.vector(tapply(d$predicted, g, mean))
  o <- as.vector(tapply(d$outcome, g, mean))
  result <- risk_group_calibration(d$predicted, d$outcome, groups)
  expect_equal(result$expected, e)
  expect_equal(result$observed, o)
  expect_equal(result$subjects, size)
  expect_equal(result$cases, as.vector(tapply(d$outcome, g, sum)))
  expect_equal(result$observed_lower, o - 1.96 * sqrt(o * (1 - o) / size))
  expect_equal(result$eo_ratio, e / o)
  expect_equal(result$eo_upper, exp(log(e / o) + 1.96 * sqrt((1 - o) / (size * o))))
  hl <- sum((o - e)^2 / (e * (1 - e) / size))
  expect_equal(result$hl$statistic, hl)
  expect_equal(result$hl$df, 5L)
  expect_equal(result$hl$p_value, stats::pchisq(hl, 5, lower.tail = FALSE))
  observed_hl <- risk_group_calibration(d$predicted, d$outcome, groups,
                                        hl_variance = "observed")$hl
  expect_equal(observed_hl$statistic, sum((o - e)^2 / (o * (1 - o) / size)))
})

test_that("the relative-risk test is the delta-method Wald test on G - 1 groups", {
  d <- calibration_data()
  groups <- risk_score_groups(d$lp, NULL, 4)
  g <- as.integer(groups)
  size <- tabulate(g)
  e <- as.vector(tapply(d$predicted, g, mean))
  o <- as.vector(tapply(d$outcome, g, mean))
  result <- risk_group_calibration(d$predicted, d$outcome, groups)
  # Jacobian of log(O_g / mean(O)) with respect to O_k
  jacobian <- diag(1 / o) - 1 / (4 * mean(o))
  covariance <- jacobian %*% diag(o * (1 - o) / size) %*% t(jacobian)
  expect_equal(result$variance_log_rr, covariance)
  difference <- log(o / mean(o)) - log(e / mean(e))
  statistic <- drop(t(difference[1:3]) %*% solve(covariance[1:3, 1:3]) %*% difference[1:3])
  expect_equal(result$rr$statistic, statistic)
  expect_equal(result$rr$df, 3L)
  expect_equal(result$observed_rr_lower,
               exp(log(o / mean(o)) - 1.96 * sqrt(diag(covariance))))
})

test_that("sampled data add the sampling variance to each group", {
  d <- calibration_data()
  set.seed(23)
  w <- ifelse(d$outcome == 1, 1, 1 / stats::runif(length(d$outcome), 0.1, 0.6))
  groups <- risk_score_groups(d$lp, w, 4)
  g <- as.integer(groups)
  group_sum <- function(x) as.vector(tapply(x, g, sum))
  size <- group_sum(w)
  e <- group_sum(w * d$predicted) / size
  o <- group_sum(w * d$outcome) / size
  result <- risk_group_calibration(d$predicted, d$outcome, groups, w)
  expect_equal(result$size, size)
  expected_variance <- e * (1 - e) / size +
    group_sum(w * (w - 1) * (d$outcome - d$predicted)^2) / size^2
  expect_equal(result$hl$statistic, sum((o - e)^2 / expected_variance))
  observed_variance <- o * (1 - o) / size +
    group_sum(w * (w - 1) * (d$outcome - o[g])^2) / size^2
  expect_equal(diag(result$variance_matrix), observed_variance)
})

test_that("groups without cases give NA intervals and an uncomputed relative-risk test", {
  predicted <- rep(c(0.01, 0.02, 0.05, 0.1), each = 50)
  outcome <- c(rep(0, 50), rep(c(1, 0, 0, 0, 0), 30))
  groups <- factor(rep(1:4, each = 50))
  expect_warning(result <- risk_group_calibration(predicted, outcome, groups),
                 "relative-risk test needs at least one case")
  expect_true(is.finite(result$hl$statistic))
  expect_true(is.na(result$rr$statistic))
  expect_equal(result$eo_ratio[1], Inf)
  expect_true(is.na(result$eo_lower[1]) && is.na(result$eo_upper[1]))
  expect_warning(expect_warning(
    observed <- risk_group_calibration(predicted, outcome, groups, hl_variance = "observed"),
    "Hosmer-Lemeshow statistic is infinite"), "relative-risk test")
  expect_equal(observed$hl$p_value, 0)
})
