# Calibration: expected and observed risks
# ========================================
#
# Notation: p_i is the predicted (expected) absolute risk of subject i, y_i
# the observed outcome over the same interval, and w_i = 1/pi_i the inverse
# inclusion probability (1 in a full cohort). Within a set of subjects g
# (everyone, or a risk-score group) with weighted size N_g = sum w_i,
#
#   E_g = sum w_i p_i / N_g     (expected risk)
#   O_g = sum w_i y_i / N_g     (observed risk).
#
# Variances. Inference is conditional on the risk factors: in a full cohort
# E_g is known and only O_g varies, with binomial variance O_g (1 - O_g) / N_g.
# In a sample from the cohort the sampling adds the variance of a
# Horvitz-Thompson total under Poisson sampling, sum (1 - pi_i)/pi_i^2 u_i^2 =
# sum w_i (w_i - 1) u_i^2, of each estimator's linearization u_i. Since both
# E_g and O_g are estimated from the sample, the linearization of
# log(E_g / O_g) is u_i = (p_i / E_g - y_i / O_g) / N_g. Earlier versions
# ignored the sampling variability of E, which made the E/O intervals of
# sampled data too narrow; in a full cohort the intervals are unchanged.
#
# Tests over G risk-score groups:
# * Hosmer-Lemeshow: HL = sum_g (O_g - E_g)^2 / V_g on G degrees of freedom.
#   By default V_g is the variance under the null hypothesis that the model is
#   calibrated, E_g (1 - E_g) / N_g plus the sampling term with
#   u_i = (y_i - p_i) / N_g (as in the classical Hosmer-Lemeshow test).
#   `hosmer.lemeshow.variance = "observed"` uses the variance of O_g (the
#   statistic of earlier versions), which is anti-conservative when groups
#   have few cases and infinite when a group has none.
# * Relative risk: the observed relative risks O_g / mean(O) are compared with
#   E_g / mean(E) by a Wald test of the log relative risks of G - 1 groups,
#   with the delta-method covariance from the variance of O_g (as in earlier
#   versions). The test needs at least one case in every group.

#' Overall expected-to-observed ratio with its 95% confidence interval
#' @param weights Inverse inclusion probabilities, or NULL for a full cohort.
#' @noRd
overall_expected_observed <- function(predicted, outcome, weights = NULL) {
  stats_eo <- expected_observed_by_group(predicted, outcome, weights,
                                         factor(rep(1L, length(outcome))))
  list(expected = stats_eo$expected, observed = stats_eo$observed,
       ratio = stats_eo$ratio, ci = c(stats_eo$lower, stats_eo$upper))
}

#' Expected and observed risks, and E/O with its interval, per group
#' @return A list of vectors, one element per level of `groups`.
#' @noRd
expected_observed_by_group <- function(predicted, outcome, weights, groups) {
  sampled <- !is.null(weights)
  if (!sampled) weights <- rep(1, length(outcome))
  index <- as.integer(groups)
  group_sum <- function(x) as.vector(rowsum(x, index, reorder = TRUE))
  size <- group_sum(weights)
  expected <- group_sum(weights * predicted) / size
  observed <- group_sum(weights * outcome) / size
  ratio <- expected / observed

  variance_log <- (1 - observed) / (size * observed)
  if (sampled) {
    u <- predicted / expected[index] - outcome / observed[index]
    variance_log <- variance_log + group_sum(weights * (weights - 1) * u^2) / size^2
  }
  half_width <- 1.96 * sqrt(variance_log)
  lower <- exp(log(ratio) - half_width)
  upper <- exp(log(ratio) + half_width)
  # Without cases the ratio is infinite and has no interval.
  no_cases <- observed == 0
  lower[no_cases] <- NA_real_
  upper[no_cases] <- NA_real_
  list(size = size, expected = expected, observed = observed, ratio = ratio,
       lower = lower, upper = upper)
}

#' Calibration within risk-score groups
#'
#' @param groups Factor of risk-score groups (all levels non-empty).
#' @param weights Inverse inclusion probabilities, or NULL for a full cohort.
#' @param hl_variance "expected" or "observed" (see the header).
#' @return A list with the per-group results, the variance matrices and the
#'   two tests.
#' @noRd
risk_group_calibration <- function(predicted, outcome, groups, weights = NULL,
                                   hl_variance = "expected") {
  sampled <- !is.null(weights)
  phase_weights <- if (sampled) weights else rep(1, length(outcome))
  index <- as.integer(groups)
  n_groups <- nlevels(groups)
  group_sum <- function(x) as.vector(rowsum(x, index, reorder = TRUE))

  eo <- expected_observed_by_group(predicted, outcome, weights, groups)
  size <- eo$size
  expected <- eo$expected
  observed <- eo$observed
  # Sampling variance term sum w (w - 1) r^2 / N^2 of a residual r per group;
  # zero in a full cohort.
  sampling_term <- function(residual) {
    group_sum(phase_weights * (phase_weights - 1) * residual^2) / size^2
  }

  # Variance of the observed risk of each group (for its confidence interval
  # and the relative-risk test).
  variance_observed <- observed * (1 - observed) / size +
    sampling_term(outcome - observed[index])
  hl_group_variance <- switch(
    hl_variance,
    expected = expected * (1 - expected) / size + sampling_term(outcome - predicted),
    observed = variance_observed
  )
  hl_statistic <- sum((observed - expected)^2 / hl_group_variance)
  if (!is.finite(hl_statistic)) {
    warn("The Hosmer-Lemeshow statistic is infinite because some risk-score groups ",
         "have no cases; consider fewer groups or `hosmer.lemeshow.variance = ",
         "\"expected\"`.")
  }
  hl_p_value <- stats::pchisq(hl_statistic, df = n_groups, lower.tail = FALSE)

  # Relative risks with respect to the mean risk over the groups
  mean_observed <- mean(observed)
  observed_rr <- observed / mean_observed
  expected_rr <- expected / mean(expected)
  variance_matrix <- diag(variance_observed, nrow = n_groups)
  off_diagonal <- -1 / (n_groups * mean_observed)
  derivative <- matrix(off_diagonal, n_groups, n_groups)
  diag(derivative) <- 1 / observed + off_diagonal
  variance_log_rr <- derivative %*% variance_matrix %*% derivative
  sd_log_rr <- sqrt(diag(variance_log_rr))

  rr_test <- relative_risk_test(observed_rr, expected_rr, derivative, variance_matrix)

  list(
    size = size,
    subjects = tabulate(index, n_groups),
    cases = as.vector(rowsum(outcome, index, reorder = TRUE)),
    expected = expected, observed = observed,
    observed_lower = observed - 1.96 * sqrt(variance_observed),
    observed_upper = observed + 1.96 * sqrt(variance_observed),
    expected_rr = expected_rr, observed_rr = observed_rr,
    observed_rr_lower = exp(log(observed_rr) - 1.96 * sd_log_rr),
    observed_rr_upper = exp(log(observed_rr) + 1.96 * sd_log_rr),
    eo_ratio = eo$ratio, eo_lower = eo$lower, eo_upper = eo$upper,
    variance_matrix = variance_matrix, variance_log_rr = variance_log_rr,
    hl = list(statistic = hl_statistic, df = n_groups, p_value = hl_p_value),
    rr = rr_test
  )
}

#' Wald test of the log relative risks of G - 1 groups
#' @noRd
relative_risk_test <- function(observed_rr, expected_rr, derivative, variance_matrix) {
  n_groups <- length(observed_rr)
  df <- n_groups - 1L
  if (any(observed_rr <= 0)) {
    warn("The relative-risk test needs at least one case in every risk-score group; ",
         sum(observed_rr <= 0), " group(s) have none, so the test is not computed. ",
         "Consider fewer groups.")
    return(list(statistic = NA_real_, df = df, p_value = NA_real_))
  }
  keep <- -n_groups
  contrast <- derivative[keep, , drop = FALSE]
  difference <- (log(observed_rr) - log(expected_rr))[keep]
  covariance <- contrast %*% variance_matrix %*% t(contrast)
  statistic <- tryCatch(
    as.numeric(crossprod(difference, solve(covariance, difference))),
    error = function(e) NA_real_
  )
  if (is.na(statistic)) {
    warn("The relative-risk test could not be computed: the covariance matrix of the ",
         "log relative risks is singular.")
    return(list(statistic = NA_real_, df = df, p_value = NA_real_))
  }
  list(statistic = statistic, df = df,
       p_value = stats::pchisq(statistic, df = df, lower.tail = FALSE))
}
