# Incidence rates of the study by attained age
# ============================================
#
# The study's age-specific incidence rate is estimated by splitting each
# subject's follow-up by single years of attained age (a Lexis diagram): for
# each age a, the number of onsets at ages in [a, a + 1) divided by the
# person-years lived at ages in [a, a + 1) between entry and the earlier of
# exit and onset. These are annual hazards on the scale of the population
# incidence rates (`model.disease.incidence.rates`), so the two can be
# compared directly. In a sample from the cohort, onsets and person-years are
# weighted by the inverse inclusion probabilities.
#
# Earlier versions counted, for each age a, onsets in [a, a + 1) among
# subjects who entered at least a year before a, divided by the number of
# those subjects still under follow-up at age a, without accounting for
# follow-up that ends within the year.
#
# Confidence intervals: exact Poisson intervals in a full cohort; in a sample
# from the cohort, log-scale Wald intervals whose variance adds the sampling
# term sum w (w - 1) (d_i - rate * t_i)^2 / T^2 to the Poisson variance D / T^2
# (d_i and t_i are subject i's onsets and person-years at that age).

#' Person-time incidence rates by single year of attained age
#'
#' @param entry,exit Ages at entry and exit.
#' @param onset Time from entry to onset (Inf without onset during follow-up).
#' @param weights Inverse inclusion probabilities, or NULL for a full cohort.
#' @return A data frame with columns `ages`, `study_incidence` (onsets per
#'   person-year), `lower`, `upper` (95% confidence interval), `cases` and
#'   `person_years` (weighted in a sample from the cohort).
#' @noRd
study_incidence_by_age <- function(entry, exit, onset, weights = NULL) {
  sampled <- !is.null(weights)
  if (!sampled) weights <- rep(1, length(entry))
  onset_age <- entry + onset
  end <- pmin(exit, onset_age)
  has_onset <- is.finite(onset_age)
  # An onset at a whole age a + 1 ends the person-time of year a, so it is
  # counted in that year.
  onset_year <- pmax(ceiling(onset_age) - 1, floor(entry))

  ages <- seq(floor(min(entry)), ceiling(max(end)) - 1)
  n_ages <- length(ages)
  cases <- person_years <- variance <- numeric(n_ages)
  for (k in seq_len(n_ages)) {
    a <- ages[k]
    time_at_age <- pmax(0, pmin(end, a + 1) - pmax(entry, a))
    onset_at_age <- has_onset & onset_year == a
    cases[k] <- sum(weights[onset_at_age])
    person_years[k] <- sum(weights * time_at_age)
    if (sampled && person_years[k] > 0) {
      rate <- cases[k] / person_years[k]
      residual <- onset_at_age - rate * time_at_age
      variance[k] <- (cases[k] + sum(weights * (weights - 1) * residual^2)) /
        person_years[k]^2
    }
  }
  observed <- person_years > 0
  ages <- ages[observed]
  cases <- cases[observed]
  person_years <- person_years[observed]
  variance <- variance[observed]
  rate <- cases / person_years

  if (sampled) {
    half_width <- 1.96 * sqrt(variance) / rate
    lower <- rate * exp(-half_width)
    upper <- rate * exp(half_width)
    lower[cases == 0] <- NA_real_
    upper[cases == 0] <- NA_real_
  } else {
    lower <- ifelse(cases > 0, stats::qchisq(0.025, 2 * cases) / 2, 0) / person_years
    upper <- stats::qchisq(0.975, 2 * (cases + 1)) / 2 / person_years
  }
  data.frame(ages = ages, study_incidence = rate, lower = lower, upper = upper,
             cases = cases, person_years = person_years)
}
