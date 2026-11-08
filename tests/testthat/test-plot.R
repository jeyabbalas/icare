validations <- example_validations()

draw_quietly <- function(figure) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  print(figure)
  invisible(TRUE)
}

test_that("the figure is a patchwork of square panels with the statistics", {
  figure <- plotModelValidation(validation.results = validations$cohort, draw = FALSE)
  expect_s3_class(figure, "patchwork")
  expect_length(figure, 4L)
  for (i in seq_along(figure)) {
    expect_equal(figure[[i]]$theme$aspect.ratio, 1)
  }
  expect_equal(figure$patches$annotation$title, "Breast cancer model")
  expect_match(figure[[1]]$labels$subtitle, "E/O")
  expect_match(figure[[3]]$labels$subtitle, "AUC")
  # calibration panels have equal axes
  for (i in 1:2) {
    ranges <- ggplot2::ggplot_build(figure[[i]])$layout$panel_params[[1]]
    expect_equal(ranges$x.range, ranges$y.range)
  }
})

test_that("cohort and sampled results draw without warnings, with every panel", {
  panels <- c("absolute", "relative", "discrimination", "incidence", "eo")
  for (result in validations) {
    figure <- plotModelValidation(validation.results = result, panels = panels, draw = FALSE)
    expect_length(figure, 5L)
    expect_no_warning(draw_quietly(figure))
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plotModelValidation(validation.results = validations$ncc))
  expect_s3_class(plot(validations$cohort, panels = "absolute", draw = FALSE), "patchwork")
})

test_that("panels, labels and limits can be chosen", {
  figure <- plotModelValidation(validation.results = validations$cohort,
                                panels = c("eo", "absolute"), dataset = "Other data",
                                model.name = "Other model", x.lim.absrisk = c(0, 5),
                                y.lim.absrisk = c(0, 5), x.lab.absrisk = "Expected",
                                draw = FALSE)
  expect_length(figure, 2L)
  expect_equal(figure$patches$annotation$title, "Other model")
  expect_match(figure$patches$annotation$subtitle, "Other data")
  expect_equal(figure[[2]]$labels$x, "Expected")
  ranges <- ggplot2::ggplot_build(figure[[2]])$layout$panel_params[[1]]
  expect_equal(ranges$x.range, c(0, 5))
})

test_that("invalid arguments are rejected", {
  result <- validations$cohort
  expect_error(plotModelValidation(validation.results = list()), "result of ModelValidation")
  expect_error(plotModelValidation(validation.results = result, panels = "roc"),
               "should be one of")
  expect_error(plotModelValidation(validation.results = result, x.lim.RR = c(0, 2)),
               "two increasing positive numbers")
  expect_error(plotModelValidation(validation.results = result, y.lim.absrisk = c(5, 1)),
               "two increasing numbers")
  expect_error(plotModelValidation(validation.results = result, base.size = 2),
               "greater than 4")
  expect_error(plotModelValidation(study.data = result$input.args$study.data[1:10, ],
                                   validation.results = result), "must be the data given")
})

test_that("absolute risks plot as a histogram", {
  bc <- load_example_data()
  risks <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 50, apply.age.interval.length = 10,
    apply.cov.profile = bc$ref_cov_dat[1:500, ])))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_warning(figure <- plot(risks))
  expect_s3_class(figure, "ggplot")
  expect_equal(figure$labels$x, "Absolute risk (%)")
})

test_that("sampled studies report weighted means in the figure and the report", {
  ncc_figure <- plotModelValidation(validation.results = validations$ncc, draw = FALSE)
  expect_match(ncc_figure$patches$annotation$subtitle, "Weighted mean follow-up")
  expect_match(ncc_figure$patches$annotation$subtitle, "Risk interval: observed follow-up")
  cohort_figure <- plotModelValidation(validation.results = validations$cohort, draw = FALSE)
  expect_match(cohort_figure$patches$annotation$subtitle, "\\| +Mean follow-up")
  expect_match(cohort_figure$patches$annotation$subtitle, "Risk interval: 5 years")
  ncc_report <- capture.output(print(validations$ncc))
  expect_true(any(grepl("^Follow-up \\(years\\): +weighted mean", ncc_report)))
  cohort_report <- capture.output(print(validations$cohort))
  expect_true(any(grepl("^Follow-up \\(years\\): +mean", cohort_report)))
})

test_that("a study with a single case can be plotted", {
  set.seed(9)
  n <- 300
  entry <- stats::runif(n, 50, 70)
  outcome <- c(1, rep(0, n - 1))
  study <- data.frame(observed.outcome = outcome, study.entry.age = entry,
                      study.exit.age = entry + 5,
                      time.of.onset = ifelse(outcome == 1, 2, Inf))
  validation <- suppressWarnings(ModelValidation(
    study, predicted.risk.interval = 5, predicted.risk = stats::runif(n, 0.01, 0.05),
    linear.predictor = stats::rnorm(n), number.of.percentiles = 2))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plotModelValidation(validation.results = validation, draw = FALSE),
                  "patchwork")
  expect_equal(density_bandwidth(1.5, 1, "nrd0"), 0.15)
  expect_equal(density_bandwidth(-0.2, 1, "SJ"), 0.1)
})
