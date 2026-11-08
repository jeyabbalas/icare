bc <- load_example_data()
inputs <- scenario_inputs(bc)
fork <- readRDS(test_path("fixtures", "fork-outputs.rds"))$validation
official <- readRDS(test_path("fixtures", "official-outputs.rds"))$validation
cohort <- inputs$cohort
ncc <- inputs$ncc

risks_of <- function(data, interval) {
  observed <- data$study.exit.age - data$study.entry.age
  followup <- if (is.null(interval)) observed else pmin(interval, observed)
  result <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = data$study.entry.age, apply.age.interval.length = followup,
    apply.cov.profile = data[, inputs$covariates], return.lp = TRUE)))
  list(risk = as.numeric(result$risk), lp = as.numeric(result$lps))
}
cohort_risk <- risks_of(cohort, 5)
ncc_risk <- risks_of(ncc, NULL)

# The quantities recorded in the fixtures (tools/fixtures/reference-outputs.R)
summarize_validation <- function(result) {
  groups <- result$Category_Results
  list(
    auc = result$AUC, variance_auc = result$Variance_AUC,
    eo = result$Overall_Expected_to_Observed_Ratio,
    eo_ci = result$CI_Overall_Expected_to_Observed_Ratio,
    hl = unname(result$Hosmer_Lemeshow_Results$statistic),
    hl_df = unname(result$Hosmer_Lemeshow_Results$parameter), hl_p = result$HL_pvalue,
    rr = unname(result$RR_test_result$statistic),
    rr_df = unname(result$RR_test_result$parameter), rr_p = result$RR_test_pvalue,
    categories = groups$Categories,
    observed = groups$Observed_Absolute_Risk, expected = groups$Predicted_Absolute_Risk,
    observed_lower = groups$CI_Absolute_Risk_Lower,
    observed_upper = groups$CI_Absolute_Risk_Upper,
    observed_rr = groups$Observed_Relative_Risk, expected_rr = groups$Predicted_Relative_Risk,
    observed_rr_lower = groups$CI_Relative_Risk_Lower,
    observed_rr_upper = groups$CI_Relative_Risk_Upper,
    variance_log_rr = result$Variance_Matrix_LogRelative_Risk,
    predicted_risk = thin(result$Subject_Specific_Predicted_Absolute_Risk),
    cases = sum(result$Subject_Specific_Observed_Outcome),
    followup_sum = sum(result$Adjusted_Followup)
  )
}

expect_fixture <- function(result, expected, skip = character(0), tolerance = 1e-10) {
  summary <- summarize_validation(result)
  for (field in setdiff(names(expected), skip)) {
    expect_equal(summary[[field]], expected[[field]], tolerance = tolerance,
                 label = field, ignore_attr = TRUE)
  }
}

validate <- function(...) ModelValidation(..., hosmer.lemeshow.variance = "observed")

test_that("full-cohort validation equals the fork and Bioconductor", {
  result <- validate(study.data = cohort, predicted.risk.interval = 5,
                     predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp)
  expect_fixture(result, fork$cohort)
  expect_fixture(result, official$cohort)
  # The default Hosmer-Lemeshow variance is the variance under calibration
  default <- ModelValidation(study.data = cohort, predicted.risk.interval = 5,
                             predicted.risk = cohort_risk$risk,
                             linear.predictor = cohort_risk$lp)
  expect_fixture(default, fork$cohort, skip = c("hl", "hl_p"))
  expect_false(isTRUE(all.equal(default$HL_pvalue, fork$cohort$hl_p)))
})

test_that("nested case-control validation equals the fork except for the E/O interval", {
  result <- validate(study.data = ncc, total.followup.validation = TRUE,
                     predicted.risk = ncc_risk$risk, linear.predictor = ncc_risk$lp)
  expect_fixture(result, fork$ncc, skip = "eo_ci")
  # E is no longer treated as fixed
  expect_equal(result$CI_Overall_Expected_to_Observed_Ratio, fork$ncc$eo_ci,
               tolerance = 1e-2)
  expect_false(isTRUE(all.equal(result$CI_Overall_Expected_to_Observed_Ratio,
                                fork$ncc$eo_ci, tolerance = 1e-6)))
})

test_that("cut-off groups use the number of groups formed (fork patch 5c)", {
  result <- validate(study.data = cohort, predicted.risk.interval = 5,
                     predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp,
                     linear.predictor.cutoffs = c(-0.4, -0.2, 0, 0.2, 0.4))
  expect_fixture(result, fork$cutoffs)
  expect_equal(nrow(result$Category_Results), 6L)
  expect_equal(unname(result$Hosmer_Lemeshow_Results$parameter), 6)
  expect_false(official$cutoffs$hl_df == 6)
})

test_that("tied risk scores count as 1/2 in the AUC (fork patch 1)", {
  tied <- round(cohort_risk$lp, 1)
  result <- validate(study.data = cohort, predicted.risk.interval = 5,
                     predicted.risk = cohort_risk$risk, linear.predictor = tied)
  expect_fixture(result, fork$tied_cohort)
  expect_gt(result$AUC, official$tied_cohort$auc)
  tied_ncc <- validate(study.data = ncc, total.followup.validation = TRUE,
                       predicted.risk = ncc_risk$risk,
                       linear.predictor = round(ncc_risk$lp, 1))
  expect_fixture(tied_ncc, fork$tied_ncc, skip = "eo_ci")
  # Bioconductor's version stops on these data
  expect_type(official$tied_ncc, "character")
})

test_that("risks computed from the model object equal those of earlier versions", {
  object <- c(covariate_model_args(bc), list(apply.cov.profile = cohort[, inputs$covariates]))
  result <- validate(study.data = cohort, predicted.risk.interval = 5,
                     iCARE.model.object = object)
  expect_fixture(result, fork$internal)
  expect_fixture(result, official$internal)
  expect_equal(result$Subject_Specific_Predicted_Absolute_Risk, cohort_risk$risk)
  expect_equal(result$Subject_Specific_Risk_Score, cohort_risk$lp)
  # settings that earlier versions needed are accepted and ignored
  expect_silent(validate(study.data = cohort, predicted.risk.interval = 5,
                         iCARE.model.object = c(object, list(use.c.code = 1, return.lp = TRUE,
                                                             return.refs.risk = TRUE))))
})

test_that("a supplied risk or risk score is kept (fork patch 5d)", {
  object <- c(covariate_model_args(bc), list(apply.cov.profile = cohort[, inputs$covariates]))
  supplied <- pmin(1, cohort_risk$risk * 1.2)
  result <- validate(study.data = cohort, predicted.risk.interval = 5,
                     predicted.risk = supplied, iCARE.model.object = object)
  expect_equal(result$Subject_Specific_Predicted_Absolute_Risk, supplied)
  expect_equal(result$Subject_Specific_Risk_Score, cohort_risk$lp)
  scores <- -cohort_risk$lp
  result <- validate(study.data = cohort, predicted.risk.interval = 5,
                     linear.predictor = scores, iCARE.model.object = object)
  expect_equal(result$Subject_Specific_Risk_Score, scores)
  expect_equal(result$Subject_Specific_Predicted_Absolute_Risk, cohort_risk$risk)
})

test_that("the model object passes the reference weights and n.imp through", {
  weights <- inputs$weights
  profiles <- cohort[1:400, inputs$covariates]
  weighted <- do.call(computeAbsoluteRisk, c(
    covariate_model_args(bc, model.ref.dataset.weights = weights),
    list(apply.age.start = cohort$study.entry.age[1:400],
         apply.age.interval.length = pmin(5, cohort$study.exit.age[1:400] -
                                             cohort$study.entry.age[1:400]),
         apply.cov.profile = profiles)))
  object <- c(covariate_model_args(bc, model.ref.dataset.weights = weights),
              list(apply.cov.profile = profiles))
  study <- cohort[1:400, ]
  study$observed.outcome[1:5] <- 1
  study$time.of.onset[1:5] <- 0.5
  # (few cases: the tests by risk group are not the point here)
  result <- suppressWarnings(validate(study.data = study, predicted.risk.interval = 5,
                                      iCARE.model.object = object))
  expect_equal(result$Subject_Specific_Predicted_Absolute_Risk, as.numeric(weighted$risk))
})

test_that("reference risks are computed over the requested windows", {
  object <- covariate_model_args(bc)
  result <- validate(study.data = cohort, predicted.risk.interval = 5,
                     predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp,
                     iCARE.model.object = object, reference.entry.age = 50,
                     reference.exit.age = 55)
  direct <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 50, apply.age.interval.length = 5,
    apply.cov.profile = bc$ref_cov_dat[1, ], return.refs.risk = TRUE)))
  expect_equal(result$Reference_Absolute_Risk, direct$refs.risk)
  design <- stats::model.matrix(bc$bc_model_formula, cbind(diagnosis = 0, bc$ref_cov_dat))
  expect_equal(result$Reference_Risk_Score,
               as.vector(design[, -1] %*% bc$bc_model_log_or[colnames(design)[-1]]))
  # one window per reference row
  entry <- rep(c(45, 55), length.out = nrow(bc$ref_cov_dat))
  per_row <- validate(study.data = cohort, predicted.risk.interval = 5,
                      predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp,
                      iCARE.model.object = object, reference.entry.age = entry,
                      reference.exit.age = entry + 10)
  odd <- seq(1, nrow(bc$ref_cov_dat), by = 2)
  at_45 <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 45, apply.age.interval.length = 10,
    apply.cov.profile = bc$ref_cov_dat[1, ], return.refs.risk = TRUE)))$refs.risk
  expect_equal(per_row$Reference_Absolute_Risk[odd], at_45[odd])
  # supplied reference values are returned as given
  given <- validate(study.data = cohort, predicted.risk.interval = 5,
                    predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp,
                    predicted.risk.ref = c(0.1, 0.2), linear.predictor.ref = c(-1, 1))
  expect_equal(given$Reference_Absolute_Risk, c(0.1, 0.2))
  expect_equal(given$Reference_Risk_Score, c(-1, 1))
})

test_that("the result has the documented structure", {
  result <- ModelValidation(study.data = cohort, predicted.risk.interval = 5,
                            predicted.risk = cohort_risk$risk,
                            linear.predictor = cohort_risk$lp, dataset = "Cohort",
                            model.name = "Model")
  expect_s3_class(result, "icareValid")
  expect_named(result, c(
    "Subject_Specific_Observed_Outcome", "Risk_Prediction_Interval", "Adjusted_Followup",
    "Subject_Specific_Predicted_Absolute_Risk", "Reference_Absolute_Risk",
    "Subject_Specific_Risk_Score", "Reference_Risk_Score", "Population_Incidence_Rate",
    "Study_Incidence_Rate", "Category_Results",
    "Category_Specific_Observed_Absolute_Risk", "Category_Specific_Predicted_Absolute_Risk",
    "Category_Specific_Observed_Relative_Risk", "Category_Specific_Predicted_Relative_Risk",
    "Category_Specific_Expected_to_Observed_Ratio", "Variance_Matrix_Absolute_Risk",
    "Variance_Matrix_LogRelative_Risk", "Hosmer_Lemeshow_Results", "HL_pvalue",
    "RR_test_result", "RR_test_pvalue", "AUC", "Variance_AUC", "CI_AUC",
    "Overall_Expected_to_Observed_Ratio", "CI_Overall_Expected_to_Observed_Ratio",
    "input.args"))
  expect_equal(result$Risk_Prediction_Interval, "5 years")
  groups <- result$Category_Results
  expect_true(all(vapply(groups[-1], is.numeric, logical(1))))
  expect_equal(sum(groups$Number_of_Subjects), nrow(cohort))
  expect_equal(sum(groups$Number_of_Cases), sum(result$Subject_Specific_Observed_Outcome))
  expect_equal(groups$Expected_to_Observed_Ratio,
               groups$Predicted_Absolute_Risk / groups$Observed_Absolute_Risk)
  expect_s3_class(result$Hosmer_Lemeshow_Results, "htest")
  expect_match(result$Hosmer_Lemeshow_Results$data.name,
               "10 risk-score groups; variance from the expected risks")
  expect_named(result$Study_Incidence_Rate,
               c("ages", "study_incidence", "lower", "upper", "cases", "person_years"))
  expect_equal(result$input.args$dataset, "Cohort")
  expect_null(result$Population_Incidence_Rate)
})

test_that("intervals per subject and the observed follow-up are supported", {
  interval <- rep(c(3, 5), length.out = nrow(cohort))
  result <- validate(study.data = cohort, predicted.risk.interval = interval,
                     predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp)
  expect_equal(result$Risk_Prediction_Interval, "Subject-specific (3 to 5 years)")
  observed <- cohort$study.exit.age - cohort$study.entry.age
  expect_equal(result$Adjusted_Followup, pmin(interval, observed))
  total <- validate(study.data = cohort, total.followup.validation = TRUE,
                    predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp)
  expect_equal(total$Risk_Prediction_Interval, "Observed Followup")
  expect_equal(total$Adjusted_Followup, observed)
  expect_equal(sum(total$Subject_Specific_Observed_Outcome), sum(cohort$observed.outcome))
})

test_that("invalid inputs give informative errors and warnings", {
  expect_error(ModelValidation(study.data = cohort, predicted.risk.interval = 5),
               "must be computed from `iCARE.model.object`")
  expect_error(ModelValidation(study.data = cohort, predicted.risk.interval = 5,
                               linear.predictor = cohort_risk$lp,
                               iCARE.model.object = covariate_model_args(bc)),
               "`predicted.risk` was not supplied")
  expect_warning(
    ModelValidation(study.data = cohort, predicted.risk.interval = 5,
                    predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp,
                    iCARE.model.object = list(model.formula = bc$bc_model_formula,
                                              model.cov.inf = bc$bc_model_cov_info)),
    "Ignoring unknown element\\(s\\) of `iCARE.model.object`: model.cov.inf")
  expect_error(ModelValidation(study.data = cohort, predicted.risk.interval = 5,
                               predicted.risk = cohort_risk$risk[-1],
                               linear.predictor = cohort_risk$lp), "one value per row")
  expect_error(ModelValidation(study.data = cohort, predicted.risk.interval = 5,
                               predicted.risk = cohort_risk$risk + 1,
                               linear.predictor = cohort_risk$lp), "between 0 and 1")
  no_cases <- cohort
  no_cases$observed.outcome <- 0
  expect_error(ModelValidation(study.data = no_cases, predicted.risk.interval = 5,
                               predicted.risk = cohort_risk$risk,
                               linear.predictor = cohort_risk$lp), "both cases and non-cases")
  expect_error(ModelValidation(study.data = cohort, predicted.risk = cohort_risk$risk,
                               linear.predictor = cohort_risk$lp),
               "Supply `predicted.risk.interval`")
  expect_error(ModelValidation(study.data = cohort, predicted.risk.interval = 5,
                               predicted.risk = cohort_risk$risk,
                               linear.predictor = cohort_risk$lp,
                               hosmer.lemeshow.variance = "pooled"), "should be one of")
})
