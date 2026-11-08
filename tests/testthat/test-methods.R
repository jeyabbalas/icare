bc <- load_example_data()

test_that("absolute risks print as a compact table", {
  risks <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 50, apply.age.interval.length = 10,
    apply.cov.profile = bc$ref_cov_dat[1:12, ])))
  output <- capture.output(print(risks))
  expect_equal(output[1], "Absolute risks of 12 profiles")
  expect_true(any(grepl("Int_Start +Int_End +Risk_Estimate", output)))
  expect_true(any(grepl("... 2 more profiles", output, fixed = TRUE)))
  capture.output(expect_invisible(print(risks, n = 2)))
  expect_equal(summary(risks), summary(risks$risk))
  split <- do.call(computeAbsoluteRiskSplitInterval, c(covariate_model_args(bc), list(
    apply.age.start = 45, apply.age.interval.length = 10, cut.time = 50,
    apply.cov.profile = bc$ref_cov_dat[1, ])))
  output <- capture.output(print(split))
  expect_equal(output[1], "Absolute risks of 1 profile")
  expect_true(any(grepl("Int_Start +Cut_Age +Int_End +Risk_Estimate", output)))
})

test_that("validation results print a report", {
  validations <- example_validations()
  output <- capture.output(print(validations$ncc))
  expect_equal(output[1], "iCARE model validation")
  for (pattern in c("Dataset: +Example nested case-control study", "Weighted subjects",
                    "Risk prediction interval: +Observed Followup", "E/O: +[0-9.]+ \\(95% CI",
                    "Hosmer-Lemeshow test: +chi-square", "^ +\\(variance from the expected risks\\)$",
                    "AUC: +0\\.[0-9]+ \\(95% CI", "E/O \\(95% CI\\)")) {
    expect_true(any(grepl(pattern, output)), label = pattern)
  }
  capture.output(expect_invisible(print(validations$cohort)))
})

test_that("the validation report fits in 80 columns and wraps the formula after '+'", {
  validations <- example_validations()
  output <- capture.output(print(validations$ncc))
  expect_true(all(nchar(output) <= 80))
  formula_text <- paste(deparse(bc$bc_model_formula, width.cutoff = 500L), collapse = " ")
  wrapped <- wrap_formula(formula_text, 53)
  expect_true(all(nchar(wrapped) <= 53))
  expect_true(all(grepl(" \\+$", wrapped[-length(wrapped)])))
  expect_equal(paste(wrapped, collapse = " "), formula_text)
  expect_equal(wrap_formula("y ~ x", 53), "y ~ x")
  long <- wrap_formula(paste0("y ~ ", strrep("x", 40), " + b"), 30)
  expect_equal(long, c("y ~", paste(strrep("x", 40), "+"), "b"))
})
