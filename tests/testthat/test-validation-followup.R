toy_study <- function() {
  data.frame(observed.outcome = c(1, 0, 1, 0, 1),
             study.entry.age = c(40, 45, 50, 55.5, 60),
             study.exit.age = c(48, 52, 53, 57, 70),
             time.of.onset = c(3, Inf, 2.5, NA, 7))
}

test_that("study data are validated", {
  study <- validate_study_data(toy_study())
  expect_equal(study$onset, c(3, Inf, 2.5, Inf, 7))
  expect_null(study$sampling_prob)
  data <- toy_study()
  expect_error(validate_study_data(as.matrix(data)[, 1:3]), "must have the columns")
  expect_error(validate_study_data(data[0, ]), "no rows")
  bad <- data
  bad$observed.outcome[2] <- 2
  expect_error(validate_study_data(bad), "1 for cases and 0 for non-cases")
  bad <- data
  bad$study.exit.age[3] <- 50
  expect_error(validate_study_data(bad), "first being row 3")
  bad <- data
  bad$time.of.onset[1] <- NA
  expect_error(validate_study_data(bad), "Every case")
  bad <- data
  bad$time.of.onset[5] <- -1
  expect_error(validate_study_data(bad), "non-negative")
  bad <- data
  bad$study.entry.age[2] <- NA
  expect_error(validate_study_data(bad), "no missing values")
  bad <- data
  bad$sampling.weights <- c(0.5, 1, 0.2, 0, 1)
  expect_error(validate_study_data(bad), "probability of inclusion")
  data$sampling.weights <- c(1, 0.1, 1, 0.2, 1)
  expect_equal(validate_study_data(data)$sampling_prob, c(1, 0.1, 1, 0.2, 1))
})

test_that("follow-up and outcomes are restricted to the prediction interval", {
  study <- validate_study_data(toy_study())
  five <- adjust_validation_followup(study, rep(5, 5))
  expect_equal(five$followup, c(5, 5, 3, 1.5, 5))
  expect_equal(five$outcome, c(1, 0, 1, 0, 0))
  observed <- adjust_validation_followup(study, NULL)
  expect_equal(observed$followup, c(8, 7, 3, 1.5, 10))
  expect_equal(observed$outcome, c(1, 0, 1, 0, 1))
})

test_that("onsets after exit are censored with a warning (fork patch 5b)", {
  data <- toy_study()
  data$time.of.onset[3] <- 4
  study <- validate_study_data(data)
  expect_warning(adjusted <- adjust_validation_followup(study, NULL),
                 "1 case\\(s\\) have `time.of.onset` after `study.exit.age`")
  expect_equal(adjusted$outcome, c(1, 0, 0, 0, 1))
  expect_equal(adjusted$onset, c(3, Inf, Inf, Inf, 7))
})

test_that("cases diagnosed at exit are kept despite rounding error in the ages", {
  # 70.94 - 61.71 is 9.2299999999999969 in floating point, below the onset
  # time 9.23; earlier versions turned such cases into non-cases.
  data <- data.frame(observed.outcome = c(1, 1, 0), study.entry.age = c(61.71, 50.1, 50),
                     study.exit.age = c(70.94, 55.2, 60), time.of.onset = c(9.23, 5.1, Inf))
  expect_lt(data$study.exit.age[1] - data$study.entry.age[1], data$time.of.onset[1])
  study <- validate_study_data(data)
  expect_no_warning(observed <- adjust_validation_followup(study, NULL))
  expect_equal(observed$outcome, c(1, 1, 0))
  expect_equal(observed$onset[1:2], c(70.94 - 61.71, 55.2 - 50.1))
  expect_equal(observed$onset[3], Inf)
  # An onset at the end of the prediction interval, up to rounding error, counts.
  interval <- adjust_validation_followup(study, c(70.94 - 61.71, 5.1 - 1e-12, 5))
  expect_equal(interval$outcome, c(1, 1, 0))
  # Onsets later than exit by more than rounding error are still censored.
  data$time.of.onset[1] <- 9.24
  expect_warning(adjust_validation_followup(validate_study_data(data), NULL),
                 "1 case\\(s\\) have `time.of.onset` after")
})

test_that("ModelValidation() keeps cases whose onset is recorded at exit", {
  set.seed(14)
  n <- 400
  entry <- round(stats::runif(n, 50, 70), 2)
  followup <- round(stats::runif(n, 1, 10), 2)
  case <- stats::rbinom(n, 1, 0.2)
  data <- data.frame(observed.outcome = case, study.entry.age = entry,
                     study.exit.age = entry + followup,
                     time.of.onset = ifelse(case == 1, followup, Inf))
  expect_no_warning(result <- ModelValidation(
    data, total.followup.validation = TRUE,
    predicted.risk = stats::runif(n, 0.05, 0.3), linear.predictor = stats::rnorm(n)))
  expect_equal(result$Subject_Specific_Observed_Outcome, case)
  expect_equal(sum(result$Study_Incidence_Rate$cases), sum(case))
})

test_that("the prediction interval is validated and described", {
  expect_null(normalize_prediction_interval(NULL, TRUE, 3))
  expect_null(normalize_prediction_interval(5, TRUE, 3))
  expect_error(normalize_prediction_interval(NULL, FALSE, 3),
               "Supply `predicted.risk.interval`")
  expect_equal(normalize_prediction_interval(5, FALSE, 3), c(5, 5, 5))
  expect_equal(normalize_prediction_interval(c(1, 2, 3), FALSE, 3), c(1, 2, 3))
  expect_error(normalize_prediction_interval(c(1, 2), FALSE, 3), "length 1 or one value")
  expect_error(normalize_prediction_interval(0, FALSE, 3), "positive")
  expect_error(normalize_prediction_interval(NA, FALSE, 3), "positive")
  expect_equal(describe_prediction_interval(NULL), "Observed Followup")
  expect_equal(describe_prediction_interval(c(5, 5)), "5 years")
  expect_equal(describe_prediction_interval(c(2, 5, 3)), "Subject-specific (2 to 5 years)")
})
