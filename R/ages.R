# Risk windows
# ============
#
# A risk window is the age interval [start, start + length). Absolute risk is
# computed year of age by year of age: the exposure of a window to the year of
# age [t, t + 1) is the overlap of the two intervals, between 0 and 1 years.
# For whole-number start ages and lengths every exposure is 0 or 1 and the
# formula is iCARE's classical yearly sum; fractional ages and lengths (for
# example the follow-up of a validation cohort) contribute the corresponding
# fraction of a year instead of being rounded up to whole years.

# Inputs closer than this (in years, about 0.3 seconds) to a whole number are
# treated as that whole number, so that sums such as 50.1 + 4.9 do not create
# a spurious exposure of 1e-15 years to an extra year of age.
age_snap_tolerance <- 1e-8

snap_to_whole <- function(x) {
  nearest <- round(x)
  ifelse(abs(x - nearest) < age_snap_tolerance, nearest, x)
}

#' Validate and recycle the start ages and lengths of the risk windows
#'
#' @param start,length `apply.age.start` and `apply.age.interval.length`.
#' @param n Number of profiles.
#' @param profile_name Name of the profile argument, for messages.
#' @return `list(start, length)`, each of length `n`.
#' @noRd
normalize_risk_windows <- function(start, length, n, profile_name = "apply.cov.profile") {
  if (is.null(start) || is.null(length)) {
    abort("`apply.age.start` and `apply.age.interval.length` must be supplied.")
  }
  if (is.data.frame(start) || is.data.frame(length)) {
    abort("`apply.age.start` and `apply.age.interval.length` must be numeric ",
          "vectors, not data frames.")
  }
  start <- as.vector(start)
  length <- as.vector(length)
  # A bare NA is logical; report it as a missing value below.
  if (is.logical(start) && all(is.na(start))) start <- as.numeric(start)
  if (is.logical(length) && all(is.na(length))) length <- as.numeric(length)
  if (!is.numeric(start) || !is.numeric(length)) {
    abort("`apply.age.start` and `apply.age.interval.length` must be numeric.")
  }
  if (base::length(start) == 1L) start <- rep(start, n)
  if (base::length(length) == 1L) length <- rep(length, n)
  if (base::length(start) != base::length(length)) {
    abort("`apply.age.start` and `apply.age.interval.length` must have the ",
          "same length (or length 1).")
  }
  if (base::length(start) != n) {
    abort("`apply.age.start` and `apply.age.interval.length` must have length ",
          "1 or one value per row of `", profile_name, "` (", n, ").")
  }
  if (anyNA(start) || anyNA(length)) {
    abort("`apply.age.start` and `apply.age.interval.length` may not contain ",
          "missing values.")
  }
  if (any(!is.finite(start)) || any(!is.finite(length))) {
    abort("`apply.age.start` and `apply.age.interval.length` must be finite.")
  }
  if (any(start < 0) || any(length < 0)) {
    abort("`apply.age.start` and `apply.age.interval.length` cannot be negative.")
  }
  start <- snap_to_whole(as.numeric(start))
  end <- snap_to_whole(start + as.numeric(length))
  list(start = start, length = end - start)
}

#' Integer ages with positive exposure in at least one window
#'
#' The window from start to start + length (right-open) is exposed to the years of age
#' floor(start), ..., ceiling(start + length) - 1. Windows of length zero are
#' exposed to none.
#' @noRd
exposed_ages <- function(start, length) {
  open <- length > 0
  if (!any(open)) {
    return(numeric(0))
  }
  first <- floor(start[open])
  last <- ceiling(start[open] + length[open]) - 1
  origin <- min(first)
  nbins <- max(last) - origin + 1
  # Difference array: +1 at each window's first year, -1 after its last.
  delta <- tabulate(first - origin + 1, nbins) -
    c(0L, tabulate(last - origin + 1, nbins)[-nbins])
  origin - 1 + which(cumsum(delta) > 0)
}
