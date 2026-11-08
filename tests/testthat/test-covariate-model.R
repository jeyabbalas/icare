bc <- load_example_data()

toy_cov_info <- list(
  list(name = "x", type = "continuous"),
  list(name = "g", type = "factor", levels = c("a", "b", "c"), ref = "b")
)

test_that("covariate information is matched by name and ordered like the formula", {
  info <- normalize_cov_info(rev(toy_cov_info), y ~ x + g)
  expect_named(info, c("x", "g"))
  # the reference level comes first
  expect_equal(info$g$levels, c("b", "a", "c"))
  # 'name' may be abbreviated, as in the example data
  expect_named(normalize_cov_info(bc$bc_model_cov_info, bc$bc_model_formula),
               all.vars(bc$bc_model_formula)[-1])
})

test_that("invalid covariate information is rejected", {
  expect_error(normalize_cov_info(toy_cov_info[1], y ~ x + g), "Not described: g")
  expect_error(normalize_cov_info(list(list(name = "x", type = "binary")), y ~ x),
               "\"continuous\" or \"factor\"")
  expect_error(normalize_cov_info(list(list(name = "g", type = "factor")), y ~ g),
               "needs a 'levels' field")
  expect_error(normalize_cov_info(list(list(name = "g", type = "factor",
                                            levels = c(1, 1))), y ~ g), "distinct")
  expect_error(normalize_cov_info(list(list(name = "g", type = "factor",
                                            levels = 1:2, ref = 3)), y ~ g),
               "'ref'.*one of its levels")
  expect_warning(normalize_cov_info(list(list(name = "x", type = "continuous",
                                              levels = 1:2)), y ~ x), "not used")
  expect_error(check_model_formula(y ~ x - 1), "intercept")
  expect_error(check_model_formula(y ~ 1), "no covariates")
})

test_that("design matrices code factors with the declared levels and reference", {
  info <- normalize_cov_info(toy_cov_info, y ~ x + g)
  reference <- prepare_covariate_data(data.frame(x = 1:6, g = rep(c("a", "b", "c"), 2)),
                                      info, "ref", allow_missing = FALSE)
  setup <- build_design_spec(y ~ x + g, reference)
  expect_equal(setup$spec$columns, c("x", "ga", "gc"))
  # profiles with a subset of the levels, and a single row, keep every column
  profile <- prepare_covariate_data(data.frame(x = 2.5, g = "c"), info, "profile",
                                    allow_missing = TRUE)
  expect_equal(build_design_matrix(setup$spec, profile),
               matrix(c(2.5, 0, 1), 1, dimnames = list(NULL, c("x", "ga", "gc"))))
  # missing covariates give NA in their columns only
  partial <- prepare_covariate_data(data.frame(x = NA, g = "a"), info, "profile",
                                    allow_missing = TRUE)
  expect_equal(unname(build_design_matrix(setup$spec, partial)[1, ]), c(NA, 1, 0))
  expect_error(prepare_covariate_data(data.frame(x = 1, g = "z"), info, "profile", TRUE),
               "outside the levels")
  expect_error(prepare_covariate_data(data.frame(x = "1", g = "a"), info, "profile", TRUE),
               "must be numeric")
  expect_error(prepare_covariate_data(data.frame(x = 1), info, "profile", TRUE),
               "lacks the covariate column")
  expect_error(prepare_covariate_data(data.frame(x = NA, g = "a"), info, "ref", FALSE),
               "cannot contain missing values")
})

test_that("the example design matrix equals model.matrix()", {
  info <- normalize_cov_info(bc$bc_model_cov_info, bc$bc_model_formula)
  reference <- prepare_covariate_data(bc$ref_cov_dat, info, "ref", allow_missing = FALSE)
  setup <- build_design_spec(bc$bc_model_formula, reference)
  expected <- stats::model.matrix(bc$bc_model_formula,
                                  cbind(diagnosis = 0, bc$ref_cov_dat))[, -1]
  expect_equal(unname(setup$design), unname(expected))
  expect_length(setup$spec$columns, 77L)
  # the interaction columns depend on both of their covariates
  interaction <- grep(":", setup$spec$columns)[1]
  expect_setequal(setup$spec$column_covariates[[interaction]], c("rd_menohrt", "bmi_dec"))
})

test_that("log relative risks are matched to the design columns by name", {
  columns <- c("x", "ga", "gc")
  beta <- c(gc = 0.3, x = 0.1, ga = -0.2)
  expect_equal(match_log_rr(beta, columns), c(x = 0.1, ga = -0.2, gc = 0.3))
  expect_equal(match_log_rr(matrix(beta, dimnames = list(names(beta), NULL)), columns),
               c(x = 0.1, ga = -0.2, gc = 0.3))
  expect_error(match_log_rr(beta[1:2], columns), "No log relative risk for: ga")
  expect_error(match_log_rr(c(beta, z = 1), columns), "z")
  expect_error(match_log_rr(unname(beta), columns), "named")
  expect_error(match_log_rr(c(x = NA, ga = 0, gc = 0), columns), "finite")
})
