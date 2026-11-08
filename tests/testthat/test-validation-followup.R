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
