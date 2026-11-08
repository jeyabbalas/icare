# Person-time incidence by brute force, subject by subject and year by year
brute_force_incidence <- function(entry, exit, onset, weights = rep(1, length(entry))) {
  onset_age <- entry + onset
  end <- pmin(exit, onset_age)
  ages <- seq(floor(min(entry)), ceiling(max(end)) - 1)
  cases <- person_years <- numeric(length(ages))
  for (i in seq_along(entry)) {
    for (k in seq_along(ages)) {
      a <- ages[k]
      person_years[k] <- person_years[k] +
        weights[i] * max(0, min(end[i], a + 1) - max(entry[i], a))
    }
    if (is.finite(onset_age[i])) {
      # onsets in (a, a + 1] end the person-time of year a
      year <- if (onset_age[i] == entry[i]) floor(entry[i]) else ceiling(onset_age[i]) - 1
      cases[ages == year] <- cases[ages == year] + weights[i]
    }
  }
  keep <- person_years > 0
  data.frame(ages = ages[keep], cases = cases[keep], person_years = person_years[keep])
}

incidence_data <- function(n = 400, seed = 30) {
  set.seed(seed)
  entry <- round(stats::runif(n, 40, 60), 2)
  exit <- entry + round(stats::runif(n, 0.5, 12), 2)
  onset <- ifelse(stats::runif(n) < 0.15, stats::runif(n) * (exit - entry), Inf)
  onset[1:5] <- ceiling(entry[1:5]) + 1 - entry[1:5]  # onsets at whole ages
  exit[1:5] <- pmax(exit[1:5], entry[1:5] + onset[1:5])
  list(entry = entry, exit = exit, onset = onset)
}

test_that("full-cohort rates are onsets per person-year with exact Poisson intervals", {
  d <- incidence_data()
  result <- study_incidence_by_age(d$entry, d$exit, d$onset)
  direct <- brute_force_incidence(d$entry, d$exit, d$onset)
  expect_equal(result$ages, direct$ages)
  expect_equal(result$cases, direct$cases)
  expect_equal(result$person_years, direct$person_years, tolerance = 1e-12)
  expect_equal(result$study_incidence, direct$cases / direct$person_years, tolerance = 1e-12)
  for (k in seq_len(nrow(result))) {
    exact <- stats::poisson.test(result$cases[k], result$person_years[k])$conf.int
    expect_equal(c(result$lower[k], result$upper[k]), as.vector(exact), tolerance = 1e-8)
  }
})

test_that("sampled data weight onsets and person-time and widen the intervals", {
  d <- incidence_data()
  set.seed(31)
  w <- ifelse(is.finite(d$onset), 1, 1 / stats::runif(length(d$entry), 0.2, 0.8))
  result <- study_incidence_by_age(d$entry, d$exit, d$onset, w)
  direct <- brute_force_incidence(d$entry, d$exit, d$onset, w)
  expect_equal(result$cases, direct$cases, tolerance = 1e-12)
  expect_equal(result$person_years, direct$person_years, tolerance = 1e-12)
  k <- which.max(result$cases)
  a <- result$ages[k]
  rate <- result$study_incidence[k]
  end <- pmin(d$exit, d$entry + d$onset)
  time <- pmax(0, pmin(end, a + 1) - pmax(d$entry, a))
  onset_age <- d$entry + d$onset
  event <- is.finite(onset_age) & pmax(ceiling(onset_age) - 1, floor(d$entry)) == a
  variance <- (result$cases[k] + sum(w * (w - 1) * (event - rate * time)^2)) /
    result$person_years[k]^2
  expect_equal(result$upper[k], rate * exp(1.96 * sqrt(variance) / rate))
  expect_true(all(is.na(result$lower[result$cases == 0])))
})

test_that("an onset at a whole age is counted in the year that ends there", {
  result <- study_incidence_by_age(entry = c(50, 50), exit = c(52, 52), onset = c(1, Inf))
  expect_equal(result$ages, c(50, 51))
  expect_equal(result$cases, c(1, 0))
  expect_equal(result$person_years, c(2, 1))
})

test_that("onsets are placed at exit and at whole ages up to rounding error", {
  # Onset times computed separately from the ages (for example from dates) can
  # put the onset a rounding error past exit or past a whole age.
  entry <- c(70.5, 69.3, 70)
  exit <- c(71, 72, 72)
  onset <- c(0.5 + 1e-12, 71 - 69.3 + 1e-12, Inf)
  expect_gt(entry[1] + onset[1], exit[1])
  result <- study_incidence_by_age(entry, exit, onset)
  expect_equal(result$ages, c(69, 70, 71))
  # Both onsets are at age 71, which ends year 70.
  expect_equal(result$cases, c(0, 2, 0))
  expect_equal(result$person_years, c(0.7, 2.5, 1), tolerance = 1e-9)
})
