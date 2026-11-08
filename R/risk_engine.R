# Absolute risk engine
# ====================
#
# The absolute risk of the disease over the window [a, a + tau) for a person
# with relative risk RR = exp(lp), baseline hazard lambda0 and competing
# mortality m is iCARE's discrete-time formula over years of age t:
#
#   AR = sum_t e_t lambda0(t) RR exp(-H(t)),
#   H(t) = sum_{u < t} e_u (lambda0(u) RR + m(u)),
#
# where e_t in [0, 1] is the exposure of the window to the year of age
# [t, t + 1) (see ages.R); for whole-number ages e_t is 1 inside the window
# and 0 outside, which is exactly the formula of earlier versions. The
# probability of being alive and free of the disease at the end of the window
# is exp(-H(a + tau)); it lets risks over adjacent windows be chained:
# AR(a, c) + exp(-H(c)) AR(c, b) = AR(a, b).
#
# The computation is one pass over the years of age, vectorized over people,
# so its cost is proportional to (number of people) x (number of years).

# Upper bound on people per vectorized pass, to cap memory use.
engine_chunk_size <- 250000L

#' Absolute risks for a set of linear predictors and windows
#'
#' @param lp Linear predictors (log relative risks), one per person.
#' @param start,length Risk windows, recycled to `length(lp)`.
#' @param baseline Calibrated baseline hazard: `list(age, hazard)`.
#' @param competing Competing mortality table (`age`, `rate`).
#' @param survival Also return the event-free survival at the end of each window.
#' @return A numeric vector of risks, or `list(risk, survival)`.
#' @noRd
compute_absolute_risk_engine <- function(lp, start, length, baseline, competing,
                                         survival = FALSE) {
  n <- base::length(lp)
  if (base::length(start) == 1L) start <- rep(start, n)
  if (base::length(length) == 1L) length <- rep(length, n)
  if (n > engine_chunk_size) {
    chunks <- split(seq_len(n), ceiling(seq_len(n) / engine_chunk_size))
    parts <- lapply(chunks, function(i) {
      compute_absolute_risk_engine(lp[i], start[i], length[i], baseline, competing,
                                   survival = TRUE)
    })
    risk <- unlist(lapply(parts, `[[`, "risk"), use.names = FALSE)
    if (!survival) {
      return(risk)
    }
    return(list(risk = risk,
                survival = unlist(lapply(parts, `[[`, "survival"), use.names = FALSE)))
  }

  risk <- numeric(n)
  cumulative_hazard <- numeric(n)
  open <- length > 0
  if (any(open)) {
    end <- start + length
    years <- seq(min(floor(start[open])), max(ceiling(end[open])) - 1)
    # Rates are only looked up for years some window is exposed to; other years
    # of the grid get 0 and contribute nothing (their exposure is 0).
    lambda0 <- baseline$hazard[match(years, baseline$age)]
    mortality <- competing$rate[match(years, competing$age)]
    lambda0[is.na(lambda0)] <- 0
    mortality[is.na(mortality)] <- 0
    rr <- exp(lp)
    for (k in seq_along(years)) {
      t <- years[k]
      exposure <- pmin(end, t + 1) - pmax(start, t)
      exposed <- which(exposure > 0)
      if (base::length(exposed) == 0L) next
      e <- exposure[exposed]
      hazard <- lambda0[k] * rr[exposed]
      h <- cumulative_hazard[exposed]
      risk[exposed] <- risk[exposed] + e * hazard * exp(-h)
      cumulative_hazard[exposed] <- h + e * (hazard + mortality[k])
    }
  }
  if (!survival) {
    return(risk)
  }
  list(risk = risk, survival = exp(-cumulative_hazard))
}
