# Direct (slow) implementations of the model, used as references in the tests

#' Absolute risk of one person by the defining sum, year by year
#'
#' AR = sum_t e_t lambda0(t) RR exp(-sum_{u < t} e_u (lambda0(u) RR + m(u))),
#' with e_t the share of year t inside the window.
reference_absolute_risk <- function(lp, start, length, baseline, competing) {
  end <- start + length
  if (length <= 0) {
    return(c(risk = 0, survival = 1))
  }
  rr <- exp(lp)
  risk <- 0
  cumulative <- 0
  for (t in seq(floor(start), ceiling(end) - 1)) {
    exposure <- min(end, t + 1) - max(start, t)
    lambda0 <- baseline$hazard[baseline$age == t]
    mortality <- competing$rate[competing$age == t]
    if (length(mortality) == 0L) mortality <- 0
    risk <- risk + exposure * lambda0 * rr * exp(-cumulative)
    cumulative <- cumulative + exposure * (lambda0 * rr + mortality)
  }
  c(risk = risk, survival = exp(-cumulative))
}

#' Baseline hazard solving the calibration equations to high precision
reference_baseline_hazard <- function(rates, lp, weights) {
  rr <- exp(lp)
  lambda0 <- rates$rate / sum(weights * rr) * sum(weights)
  for (iteration in 1:200) {
    cumulative <- c(0, cumsum(lambda0))[seq_along(lambda0)]
    expected <- vapply(cumulative, function(h) {
      s <- exp(-rr * h)
      sum(weights * s * rr) / sum(weights * s)
    }, numeric(1))
    lambda0 <- rates$rate / expected
  }
  lambda0
}

#' Donors of a profile with missing risk factors, by direct search
#'
#' Population rows (reference rows, with their simulated SNPs) in the
#' profile's percentile bin of the linear predictor of its observed design
#' columns, widened while the bin has no positive weight.
reference_donors <- function(model, design_row) {
  population <- model$population
  weights <- population$weights
  present <- which(!is.na(design_row))
  design <- population_design(population)
  ref_lp <- as.vector(design[, present, drop = FALSE] %*% model$beta[present])
  observed_lp <- sum(design_row[present] * model$beta[present])
  cuts <- get_cutpoints(ref_lp, seq(0, 1, 0.01))
  bin <- get_miss_perc(cuts, observed_lp)
  in_range <- function(low, high) ref_lp >= low & ref_lp < high
  has_donors <- function(low, high) sum(weights[in_range(low, high)]) > 0
  range <- cuts[c(bin, bin + 1)]
  if (!has_donors(range[1], range[2])) {
    range <- get_new_endpoints(cuts, bin, has_donors)
  }
  which(in_range(range[1], range[2]) & weights > 0)
}

#' Weighted mean over donors of the risks of the completed linear predictors
reference_average_risk <- function(model, donors, lp, start, length) {
  risks <- vapply(lp, function(l) {
    reference_absolute_risk(l, start, length, model$baseline, model$competing)[["risk"]]
  }, numeric(1))
  w <- model$population$weights[donors]
  c(risk = sum(w * risks) / sum(w), lp = sum(w * lp) / sum(w))
}

#' Risk of one covariate profile with missing values, by brute force
#'
#' The profile's observed covariates are completed with each donor's values of
#' the missing ones, and the completed design rows are built by model.matrix().
reference_missing_covariate_risk <- function(model, formula, profile, start, length) {
  design_of <- function(data) {
    terms <- stats::delete.response(stats::terms(formula))
    frame <- stats::model.frame(terms, data, na.action = stats::na.pass)
    stats::model.matrix(terms, frame)[, -1, drop = FALSE]
  }
  design_row <- design_of(profile)[1, names(model$beta)]
  donors <- reference_donors(model, design_row)
  missing <- names(profile)[vapply(profile, function(x) is.na(x[1]), logical(1))]
  completed <- profile[rep(1, length(donors)), , drop = FALSE]
  rows <- model$population$reference_row[donors]
  for (v in missing) completed[[v]] <- model$covariates$data[[v]][rows]
  lp <- as.vector(design_of(completed)[, names(model$beta), drop = FALSE] %*% model$beta)
  reference_average_risk(model, donors, lp, start, length)
}

#' Risk of one profile whose missing design columns are each fully missing
#' (no interaction of an observed and a missing risk factor), by brute force
reference_missing_column_risk <- function(model, design_row, start, length) {
  donors <- reference_donors(model, design_row)
  missing <- which(is.na(design_row))
  completed <- matrix(design_row, nrow = length(donors), ncol = length(design_row),
                      byrow = TRUE)
  completed[, missing] <- population_design(model$population)[donors, missing]
  reference_average_risk(model, donors, as.vector(completed %*% model$beta), start, length)
}
