# Monte Carlo study of the intervals and tests of ModelValidation()
# =================================================================
#
# Simulates validation studies of a model that is exactly calibrated and
# records how often the 95% confidence intervals contain the true value and
# how often the tests reject at the 5% level. The results are summarized in
# the "Statistical details" section of vignette("model-validation").
#
# Usage, from the root of the repository with this distribution installed:
#
#   Rscript tools/simulations/validation-simulation.R [scale] [cores] [output.rds]
#
# `scale` multiplies the default numbers of replicates (5,000 per full-cohort
# setting, 2,000 per sampled design, 1,500 per incidence setting); `cores`
# defaults to all cores but one (1 on Windows); the replicate-level results
# are saved to `output.rds` when it is given. Every replicate has its own
# random-number stream, so the results do not depend on the number of cores.
# The default run takes about 9 minutes of computing time (under a minute on
# 16 cores).
#
# Data-generating model. The risk score is s ~ N(0, 0.5^2) and the 5-year
# risk is p = c exp(s), with c giving the chosen mean risk; the outcome is
# y ~ Bernoulli(p). Entry ages are uniform on [50, 70) with 5 years of
# follow-up, and onsets are uniform over follow-up. The predicted risks are
# therefore exactly calibrated, and the true AUC of the score has the closed
# form (Phi(sd / sqrt(2)) - mean / 2) / (1 - mean), because the scores of
# cases are N(sd^2, sd^2).

suppressPackageStartupMessages(library(iCARE))
options(width = 160)

arguments <- commandArgs(trailingOnly = TRUE)
scale <- if (length(arguments) >= 1) as.numeric(arguments[1]) else 1
cores <- if (length(arguments) >= 2) as.integer(arguments[2]) else
  max(1L, parallel::detectCores() - 1L)
output_file <- if (length(arguments) >= 3) arguments[3] else NULL
if (.Platform$OS.type == "windows") cores <- 1L

score_sd <- 0.5
follow_up <- 5

simulate_cohort <- function(size, mean_risk) {
  score <- stats::rnorm(size, 0, score_sd)
  risk <- pmin(mean_risk / exp(score_sd^2 / 2) * exp(score), 1)
  outcome <- stats::rbinom(size, 1, risk)
  entry <- stats::runif(size, 50, 70)
  data.frame(
    observed.outcome = outcome,
    study.entry.age = entry,
    study.exit.age = entry + follow_up,
    time.of.onset = ifelse(outcome == 1, stats::runif(size, 0, follow_up), Inf),
    score = score,
    risk = risk
  )
}

true_auc <- function(mean_risk) {
  (stats::pnorm(score_sd / sqrt(2)) - mean_risk / 2) / (1 - mean_risk)
}

validate <- function(study, hl_variance = "expected") {
  columns <- intersect(c("observed.outcome", "study.entry.age", "study.exit.age",
                         "time.of.onset", "sampling.weights"), names(study))
  suppressWarnings(ModelValidation(
    study.data = study[, columns],
    predicted.risk.interval = follow_up,
    predicted.risk = study$risk,
    linear.predictor = study$score,
    hosmer.lemeshow.variance = hl_variance
  ))
}

covers <- function(interval, value) interval[1] <= value && value <= interval[2]

# One row of results per replicate
summarize_validation <- function(result, mean_risk, observed_hl = NULL) {
  c(
    cases = sum(result$Category_Results$Number_of_Cases),
    eo_covers = covers(result$CI_Overall_Expected_to_Observed_Ratio, 1),
    auc_covers = covers(result$CI_AUC, true_auc(mean_risk)),
    hl_rejects = result$HL_pvalue < 0.05,
    hl_observed_rejects = if (is.null(observed_hl)) NA else observed_hl$HL_pvalue < 0.05,
    rr_rejects = result$RR_test_pvalue < 0.05   # NA when a group has no case
  )
}

# Replicates in parallel, each with its own L'Ecuyer-CMRG stream
run_replicates <- function(count, replicate, seed) {
  count <- max(1L, round(count * scale))
  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1]), add = TRUE)
  set.seed(seed)
  streams <- vector("list", count)
  stream <- .Random.seed
  for (i in seq_len(count)) {
    streams[[i]] <- stream
    stream <- parallel::nextRNGStream(stream)
  }
  rows <- parallel::mclapply(seq_len(count), function(i) {
    assign(".Random.seed", streams[[i]], envir = globalenv())
    replicate()
  }, mc.cores = cores)
  failed <- vapply(rows, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("Replicates failed: ", rows[[which(failed)[1]]])
  as.data.frame(do.call(rbind, rows))
}

# Percentage with its Monte Carlo standard error
percent <- function(x) {
  x <- x[!is.na(x)]
  p <- mean(x)
  sprintf("%.1f (%.1f)", 100 * p, 100 * sqrt(p * (1 - p) / length(x)))
}

# ---------------------------------------------------------------------------
# 1. Full cohorts of 5,000 with 10 risk-score groups, about 10, 30 and 50
#    cases per group on average
# ---------------------------------------------------------------------------

full_settings <- data.frame(mean_risk = c(0.02, 0.06, 0.10))
full_results <- lapply(seq_len(nrow(full_settings)), function(k) {
  mean_risk <- full_settings$mean_risk[k]
  run_replicates(5000, function() {
    study <- simulate_cohort(5000, mean_risk)
    summarize_validation(validate(study), mean_risk,
                         observed_hl = validate(study, "observed"))
  }, seed = 1000 + k)
})

full_table <- data.frame(
  `Cases per group` = vapply(full_results, function(r) mean(r$cases) / 10, numeric(1)),
  Replicates = vapply(full_results, nrow, integer(1)),
  `E/O CI coverage` = vapply(full_results, function(r) percent(r$eo_covers), ""),
  `AUC CI coverage` = vapply(full_results, function(r) percent(r$auc_covers), ""),
  `HL size` = vapply(full_results, function(r) percent(r$hl_rejects), ""),
  `HL size, observed variance` = vapply(full_results, function(r) percent(r$hl_observed_rejects), ""),
  `RR test size` = vapply(full_results, function(r) percent(r$rr_rejects), ""),
  `RR test not computed` = vapply(full_results, function(r) percent(is.na(r$rr_rejects)), ""),
  check.names = FALSE
)

# ---------------------------------------------------------------------------
# 2. Samples from cohorts of 50,000: all cases and non-cases included with
#    probability 0.1, or with a probability rising or falling with the score
#    (known inclusion probabilities, passed as sampling.weights)
# ---------------------------------------------------------------------------

inclusion_designs <- list(
  `10% of non-cases` = function(score) rep(0.1, length(score)),
  `high scores oversampled` = function(score) stats::plogis(stats::qlogis(0.1) + score),
  `low scores oversampled` = function(score) stats::plogis(stats::qlogis(0.1) - score)
)
sampled_settings <- expand.grid(design = names(inclusion_designs),
                                mean_risk = c(0.01, 0.05), stringsAsFactors = FALSE)
sampled_results <- lapply(seq_len(nrow(sampled_settings)), function(k) {
  mean_risk <- sampled_settings$mean_risk[k]
  inclusion <- inclusion_designs[[sampled_settings$design[k]]]
  run_replicates(2000, function() {
    cohort <- simulate_cohort(50000, mean_risk)
    probability <- ifelse(cohort$observed.outcome == 1, 1, inclusion(cohort$score))
    included <- stats::runif(nrow(cohort)) < probability
    study <- cohort[included, ]
    study$sampling.weights <- probability[included]
    summarize_validation(validate(study), mean_risk)
  }, seed = 2000 + k)
})

sampled_table <- data.frame(
  Design = sampled_settings$design,
  `Mean risk` = sampled_settings$mean_risk,
  `Cases per group` = vapply(sampled_results, function(r) mean(r$cases) / 10, numeric(1)),
  Replicates = vapply(sampled_results, nrow, integer(1)),
  `E/O CI coverage` = vapply(sampled_results, function(r) percent(r$eo_covers), ""),
  `AUC CI coverage` = vapply(sampled_results, function(r) percent(r$auc_covers), ""),
  `HL size` = vapply(sampled_results, function(r) percent(r$hl_rejects), ""),
  `RR test size` = vapply(sampled_results, function(r) percent(r$rr_rejects), ""),
  check.names = FALSE
)

# ---------------------------------------------------------------------------
# 3. Study incidence by attained age, ages 50 to 74: a constant hazard, a
#    hazard rising 6% a year (constant within each year of age), and a
#    constant hazard in a sample of all cases and 10% of non-cases
# ---------------------------------------------------------------------------

ages <- 50:74
rising_hazard <- 0.005 * exp(0.06 * (ages - 50))   # hazard on [a, a + 1)
cumulative_hazard <- function(age) {
  stats::approx(c(ages, 75), c(0, cumsum(rising_hazard)), xout = age, rule = 2)$y
}
age_at_cumulative_hazard <- function(hazard) {
  stats::approx(c(0, cumsum(rising_hazard)), c(ages, 75), xout = hazard, rule = 2)$y
}

incidence_settings <- list(
  `constant hazard 0.01, cohort of 5,000` =
    list(size = 5000, hazard = 0.01, sampled = FALSE),
  `hazard rising 6% a year, cohort of 10,000` =
    list(size = 10000, hazard = NULL, sampled = FALSE),
  `constant hazard 0.002, sample from 50,000` =
    list(size = 50000, hazard = 0.002, sampled = TRUE)
)
incidence_results <- lapply(seq_along(incidence_settings), function(k) {
  setting <- incidence_settings[[k]]
  truth <- if (is.null(setting$hazard)) rising_hazard else rep(setting$hazard, length(ages))
  run_replicates(1500, function() {
    entry <- stats::runif(setting$size, 50, 70)
    if (is.null(setting$hazard)) {
      onset <- age_at_cumulative_hazard(cumulative_hazard(entry) +
                                          stats::rexp(setting$size)) - entry
      risk <- 1 - exp(-(cumulative_hazard(entry + follow_up) - cumulative_hazard(entry)))
    } else {
      onset <- stats::rexp(setting$size, setting$hazard)
      risk <- rep(1 - exp(-follow_up * setting$hazard), setting$size)
    }
    outcome <- as.numeric(onset <= follow_up)
    study <- data.frame(observed.outcome = outcome, study.entry.age = entry,
                        study.exit.age = entry + follow_up,
                        time.of.onset = ifelse(outcome == 1, onset, Inf),
                        score = stats::rnorm(setting$size), risk = risk)
    if (setting$sampled) {
      probability <- ifelse(outcome == 1, 1, 0.1)
      included <- stats::runif(setting$size) < probability
      study <- study[included, ]
      study$sampling.weights <- probability[included]
    }
    rates <- validate(study)$Study_Incidence_Rate
    rates <- rates[match(ages, rates$ages), ]
    c(rate_ratio = rates$study_incidence / truth,
      covers = rates$lower <= truth & truth <= rates$upper)
  }, seed = 3000 + k)
})

incidence_table <- data.frame(
  Setting = names(incidence_settings),
  Replicates = vapply(incidence_results, nrow, integer(1)),
  `Mean estimate / truth, by age` = vapply(incidence_results, function(r) {
    ratio <- colMeans(r[, grep("^rate_ratio", names(r))], na.rm = TRUE)
    sprintf("%.3f to %.3f", min(ratio), max(ratio))
  }, ""),
  `CI coverage, by age` = vapply(incidence_results, function(r) {
    coverage <- 100 * colMeans(r[, grep("^covers", names(r))], na.rm = TRUE)
    sprintf("%.1f to %.1f", min(coverage), max(coverage))
  }, ""),
  `CI coverage, all ages` = vapply(incidence_results, function(r) {
    percent(unlist(r[, grep("^covers", names(r))]))
  }, ""),
  check.names = FALSE
)

# ---------------------------------------------------------------------------

cat("ModelValidation() under exact calibration: iCARE",
    as.character(utils::packageVersion("iCARE")), "\n")
cat("Percentages with their Monte Carlo standard errors; sizes at the 5% level.\n")
cat("RR test size is among the replicates in which the test was computed.\n\n")
cat("1. Full cohorts of 5,000, 10 risk-score groups\n\n")
print(full_table, row.names = FALSE, digits = 3)
cat("\n2. Samples from cohorts of 50,000 (all cases)\n\n")
print(sampled_table, row.names = FALSE, digits = 3)
cat("\n3. Study incidence by attained age, ages 50-74\n\n")
print(incidence_table, row.names = FALSE)

if (!is.null(output_file)) {
  saveRDS(list(full = full_results, sampled = sampled_results,
               incidence = incidence_results,
               settings = list(full = full_settings, sampled = sampled_settings,
                               incidence = incidence_settings)),
          output_file)
}
