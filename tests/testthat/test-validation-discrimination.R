# AUC and DeLong variance by enumerating every case-control pair
pairwise_auc <- function(score, outcome, weights = rep(1, length(score))) {
  case <- outcome == 1
  kernel <- outer(score[case], score[!case], function(a, b) (a > b) + 0.5 * (a == b))
  wc <- weights[case]
  wn <- weights[!case]
  auc <- sum(outer(wc, wn) * kernel) / (sum(wc) * sum(wn))
  list(auc = auc, case_placement = as.vector(kernel %*% wn) / sum(wn),
       control_placement = as.vector(wc %*% kernel) / sum(wc))
}

test_that("weighted shares below count ties as one half", {
  reference <- c(1, 2, 2, 3)
  weights <- c(1, 2, 3, 4)
  expect_equal(weighted_share_below(c(0, 1, 2, 2.5, 3, 9), reference, weights),
               c(0, 0.5, 3.5, 6, 8, 10) / 10)
})

test_that("the full-cohort AUC and its DeLong variance equal the pairwise definitions", {
  set.seed(10)
  for (scores in list(stats::rnorm(400), round(stats::rnorm(400), 1))) {
    outcome <- stats::rbinom(400, 1, stats::plogis(scores - 1.5))
    direct <- pairwise_auc(scores, outcome)
    result <- compute_auc(scores, outcome)
    expect_equal(result$auc, direct$auc, tolerance = 1e-12)
    expect_equal(result$variance,
                 stats::var(direct$case_placement) / sum(outcome) +
                   stats::var(direct$control_placement) / sum(1 - outcome),
                 tolerance = 1e-12)
    expect_equal(result$ci, result$auc + c(-1, 1) * 1.96 * sqrt(result$variance))
  }
})

test_that("the weighted AUC is the weighted Mann-Whitney statistic", {
  set.seed(11)
  scores <- round(stats::rnorm(300), 1)
  outcome <- stats::rbinom(300, 1, 0.3)
  prob <- ifelse(outcome == 1, 1, stats::runif(300, 0.05, 0.5))
  direct <- pairwise_auc(scores, outcome, 1 / prob)
  result <- compute_auc(scores, outcome, prob)
  expect_equal(result$auc, direct$auc, tolerance = 1e-12)
  # Variance as in earlier versions: Hmisc::wtd.var() of the placements plus
  # the (1 - pi) / pi sampling term
  wtd_var <- function(x, w) sum(w * (x - sum(w * x) / sum(w))^2) / (sum(w) - 1)
  part <- function(placement, p) {
    w <- 1 / p
    (wtd_var(placement, w) + sum(w * (placement - direct$auc)^2 * (1 - p) / p) / sum(w)) /
      sum(w)
  }
  case <- outcome == 1
  expect_equal(result$variance, part(direct$case_placement, prob[case]) +
                 part(direct$control_placement, prob[!case]), tolerance = 1e-12)
  # integer weights equal replicated subjects
  counts <- sample(1:3, 300, replace = TRUE)
  expect_equal(compute_auc(scores, outcome, 1 / counts)$auc,
               compute_auc(rep(scores, counts), rep(outcome, counts))$auc,
               tolerance = 1e-12)
})

test_that("the AUC is invariant to increasing transformations and 1/2 without information", {
  set.seed(12)
  scores <- stats::rnorm(200)
  outcome <- stats::rbinom(200, 1, 0.2)
  expect_equal(compute_auc(exp(3 * scores), outcome)$auc, compute_auc(scores, outcome)$auc)
  expect_equal(compute_auc(rep(1, 200), outcome)$auc, 0.5)
  expect_equal(compute_auc(-scores, outcome)$auc, 1 - compute_auc(scores, outcome)$auc)
})
