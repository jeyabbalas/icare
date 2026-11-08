test_that("weighted quantiles interpolate the weighted empirical distribution", {
  set.seed(5)
  x <- stats::rnorm(200)
  probs <- seq(0, 1, 0.1)
  # unit weights without ties: linear interpolation of the ECDF (type 4)
  expect_equal(weighted_quantile_interpolated(x, rep(1, 200), probs),
               stats::quantile(x, probs, type = 4, names = FALSE))
  w <- stats::runif(200)
  expect_equal(weighted_quantile_interpolated(x, 7 * w, probs),
               weighted_quantile_interpolated(x, w, probs))
  # zero weights drop their values
  expect_equal(weighted_quantile_interpolated(c(x, 100), c(w, 0), probs),
               weighted_quantile_interpolated(x, w, probs))
  # ties: the ECDF steps at the distinct values 1, 2, 3 to 0.5, 0.75, 1
  expect_equal(weighted_quantile_interpolated(c(1, 1, 2, 3), rep(1, 4),
                                              c(0, 0.25, 0.6, 0.875, 1)),
               c(1, 1, 1.4, 2.5, 3))
})

test_that("risk-score groups equal those of the fork, with and without ties", {
  fork <- readRDS(test_path("fixtures", "fork-outputs.rds"))
  cases <- group_cases()
  for (i in seq_along(cases)) {
    groups <- droplevels(quantile_groups(cases[[i]]$x, cases[[i]]$weights, cases[[i]]$q))
    expect_identical(as.integer(groups), fork$groups[[i]]$code, label = paste("case", i))
    expect_identical(levels(groups), fork$groups[[i]]$labels, label = paste("case", i))
  }
})

test_that("untied scores are cut at their quantiles", {
  set.seed(6)
  x <- stats::rnorm(500)
  expect_identical(risk_score_groups(x, NULL, 5),
                   cut(x, stats::quantile(x, seq(0, 1, 0.2)), include.lowest = TRUE))
  groups <- risk_score_groups(x, stats::runif(500, 1, 20), 10)
  expect_equal(nlevels(groups), 10L)
  expect_false(anyNA(groups))
})

test_that("tied scores form their own groups and every subject keeps a group", {
  x <- c(rep(0, 60), 1:40)
  groups <- risk_score_groups(x, NULL, 10)
  expect_false(anyNA(groups))
  expect_equal(levels(groups)[1], "0")
  expect_true(all(groups[x == 0] == "0"))
  expect_equal(sum(table(groups)), 100L)
})

test_that("groups at cut-offs ignore cut-offs outside the scores", {
  x <- seq(-1, 1, length.out = 101)
  groups <- risk_score_groups(x, cutoffs = c(-0.5, 0, 0.5))
  expect_equal(levels(groups), c("[-1,-0.5]", "(-0.5,0]", "(0,0.5]", "(0.5,1]"))
  expect_equal(as.vector(table(groups)), c(26L, 25L, 25L, 25L))
  outside <- risk_score_groups(x, cutoffs = c(-5, -0.5, 0, 0.5, 1, 5))
  expect_identical(outside, groups)
  expect_error(risk_score_groups(x, cutoffs = c(-5, 5)), "1 non-empty group")
  expect_error(risk_score_groups(rep(0.3, 10)), "1 non-empty group")
  expect_error(risk_score_groups(x, cutoffs = c(0, NA)), "no missing values")
})
