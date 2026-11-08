# Calibration of the baseline hazard
# ==================================
#
# Under the proportional-hazards model lambda(t | z) = lambda0(t) RR(z), the
# disease incidence rate of the population at age t is the mean hazard among
# the people still free of disease at t:
#
#   lambda(t) = lambda0(t) E[RR(Z) | T >= t],
#   E[RR(Z) | T >= t] = sum_i w_i S_i(t) RR_i / sum_i w_i S_i(t),
#   S_i(t) = exp(-RR_i Lambda0(t)),   Lambda0(t) = sum_{u < t} lambda0(u),
#
# where the sums run over the reference population with weights w_i and the
# cumulative hazard starts at the youngest age of the incidence table
# (competing mortality is assumed independent of the risk factors and cancels).
# lambda0 = lambda / E[RR | T >= t] is found by fixed-point iteration,
# starting from the rare-disease approximation E[RR | T >= t] = E[RR], and
# stopping when the summed absolute change of E[RR | T >= t] over ages is at
# most `tolerance` (the rule of earlier versions of iCARE).

#' Calibrate the baseline hazard to the population incidence rates
#'
#' @param incidence Yearly incidence table (`age`, `rate`), consecutive ages.
#' @param population_lp Linear predictor of each reference-population row.
#' @param weights Non-negative weights of the rows (any positive scale).
#' @param tolerance,max_iterations Convergence controls.
#' @return `list(age, hazard, expected_rr, iterations)`: the calibrated
#'   baseline hazard lambda0 and the survivor-weighted mean relative risk at
#'   each age of the incidence table.
#' @noRd
calibrate_baseline_hazard <- function(incidence, population_lp, weights,
                                      tolerance = 1e-3, max_iterations = 1000L) {
  rr <- exp(population_lp)
  if (any(!is.finite(rr))) {
    abort("The relative risks of the reference population are not finite; ",
          "check `model.log.RR` and `model.ref.dataset`.")
  }
  lambda <- incidence$rate
  n_ages <- length(lambda)
  weighted_rr <- weights * rr

  previous <- weighted_average(rr, weights) - 1
  expected <- weighted_average(rr, weights)
  iterations <- 0L
  while (sum(abs(expected - previous)) > tolerance) {
    iterations <- iterations + 1L
    if (iterations > max_iterations) {
      abort("The calibration of the baseline hazard did not converge in ",
            max_iterations, " iterations.")
    }
    previous <- expected
    lambda0 <- lambda / previous
    # Cumulative baseline hazard up to, but not including, each age.
    cumulative <- c(0, cumsum(lambda0))[seq_len(n_ages)]
    expected <- numeric(n_ages)
    for (k in seq_len(n_ages)) {
      survival <- exp(-rr * cumulative[k])
      expected[k] <- sum(weighted_rr * survival) / sum(weights * survival)
    }
  }
  list(age = incidence$age, hazard = lambda / expected, expected_rr = expected,
       iterations = iterations)
}
