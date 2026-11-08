# Discrimination: area under the ROC curve
# ========================================
#
# The AUC is the probability that a case has a higher risk score than a
# non-case, counting ties as 1/2 (the Mann-Whitney statistic):
#
#   AUC = sum_ij w_i v_j K(s_i, s_j) / (sum_i w_i sum_j v_j),
#   K(s_case, s_noncase) = 1 if s_case > s_noncase, 1/2 if tied, 0 otherwise,
#
# with weights w (cases) and v (non-cases) equal to 1 in a full cohort and to
# the inverse inclusion probabilities 1/pi in a sample from the cohort (e.g. a
# nested case-control study). Earlier versions counted ties as 0, which biases
# the AUC down when risk scores are tied.
#
# The variance comes from the placements (DeLong et al., 1988): a case's
# placement is the weighted share of non-cases it outranks, and a non-case's
# placement the weighted share of cases that outrank it. Placements are
# computed from sorted scores and cumulative weights in O(n log n) time and
# O(n) memory; earlier versions built the full cases x non-cases matrix.
# The weighted variance adds a (1 - pi)/pi term for the sampling, as in
# earlier versions.

#' Weighted share of `reference` values below each `x`, ties counting 1/2
#' @noRd
weighted_share_below <- function(x, reference, weights) {
  ord <- order(reference)
  sorted <- reference[ord]
  cumulative <- c(0, cumsum(weights[ord]))
  below <- cumulative[findInterval(x, sorted, left.open = TRUE) + 1L]
  at_or_below <- cumulative[findInterval(x, sorted) + 1L]
  (below + at_or_below) / (2 * cumulative[length(cumulative)])
}

#' Weighted variance with the "unbiased" denominator sum(w) - 1 (the
#' definition of Hmisc::wtd.var() used by earlier versions)
#' @noRd
weighted_variance <- function(x, w) {
  mean_x <- sum(w * x) / sum(w)
  sum(w * (x - mean_x)^2) / (sum(w) - 1)
}

#' AUC of a risk score with its variance and 95% confidence interval
#'
#' @param score Risk scores.
#' @param outcome 0/1 outcomes; there must be cases and non-cases.
#' @param sampling_prob Inclusion probabilities, or NULL for a full cohort.
#' @return `list(auc, variance, ci)`.
#' @noRd
compute_auc <- function(score, outcome, sampling_prob = NULL) {
  cases <- outcome == 1
  weights <- if (is.null(sampling_prob)) rep(1, length(score)) else 1 / sampling_prob
  case_score <- score[cases]
  control_score <- score[!cases]
  case_weight <- weights[cases]
  control_weight <- weights[!cases]

  # Share of non-cases each case outranks, and share of cases outranking each
  # non-case.
  case_placement <- weighted_share_below(case_score, control_score, control_weight)
  control_placement <- 1 - weighted_share_below(control_score, case_score, case_weight)
  auc <- sum(case_weight * case_placement) / sum(case_weight)

  variance <- if (is.null(sampling_prob)) {
    stats::var(case_placement) / length(case_placement) +
      stats::var(control_placement) / length(control_placement)
  } else {
    case_prob <- sampling_prob[cases]
    control_prob <- sampling_prob[!cases]
    case_term <- weighted_variance(case_placement, case_weight) +
      sum(case_weight * (case_placement - auc)^2 * (1 - case_prob) / case_prob) /
      sum(case_weight)
    control_term <- weighted_variance(control_placement, control_weight) +
      sum(control_weight * (control_placement - auc)^2 * (1 - control_prob) /
            control_prob) / sum(control_weight)
    case_term / sum(case_weight) + control_term / sum(control_weight)
  }
  list(auc = auc, variance = variance,
       ci = auc + c(-1, 1) * 1.96 * sqrt(variance))
}
