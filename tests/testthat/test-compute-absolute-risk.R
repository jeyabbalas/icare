bc <- load_example_data()
inputs <- scenario_inputs(bc)
official <- readRDS(test_path("fixtures", "official-outputs.rds"))
fork <- readRDS(test_path("fixtures", "fork-outputs.rds"))
profiles <- bc$ref_cov_dat[inputs$profile_rows, ]

complete_risks <- function(..., profiles_used = profiles, ages = inputs$ages,
                           lengths = inputs$lengths) {
  do.call(computeAbsoluteRisk, c(covariate_model_args(bc, ...), list(
    apply.age.start = ages, apply.age.interval.length = lengths,
    apply.cov.profile = profiles_used, return.lp = TRUE, return.refs.risk = TRUE)))
}

test_that("risks of complete profiles equal those of earlier versions", {
  result <- complete_risks()
  expected <- official$complete_profiles
  expect_equal(as.numeric(result$risk), expected$risk, tolerance = 1e-12)
  expect_equal(as.numeric(result$lps), expected$lps, tolerance = 1e-12)
  expect_equal(thin(result$refs.risk), expected$refs_risk, tolerance = 1e-12)
  expect_equal(as.vector(result$beta.used), as.vector(expected$beta))
  expect_equal(rownames(result$beta.used), names(bc$bc_model_log_or))
})

test_that("non-uniform reference weights calibrate as in the fork (fork patch 2)", {
  result <- complete_risks(model.ref.dataset.weights = inputs$weights)
  expect_equal(as.numeric(result$risk), fork$weighted_reference$risk, tolerance = 1e-12)
  # Bioconductor's version scrambles the weights in the calibration
  expect_gt(max(abs(as.numeric(result$risk) / official$weighted_reference$risk - 1)), 1e-3)
  # weights are normalized: their scale does not matter
  rescaled <- complete_risks(model.ref.dataset.weights = 1000 * inputs$weights)
  expect_equal(rescaled$risk, result$risk, tolerance = 1e-14)
})

test_that("the SNP model reproduces earlier versions for genotyped profiles", {
  set.seed(1)
  result <- do.call(computeAbsoluteRisk, c(snp_model_args(bc), list(
    apply.age.start = 50, apply.age.interval.length = 30,
    apply.cov.profile = bc$new_cov_prof, apply.snp.profile = bc$new_snp_prof,
    return.lp = TRUE, return.refs.risk = TRUE)))
  # profile 3 has every SNP; the value of the earlier versions' unit test
  expect_equal(result$risk[3], official$snp_example$risk[3], tolerance = 1e-12)
  expect_equal(result$risk[3], 0.16893543, tolerance = 1e-7)
  expect_equal(result$lps[3], official$snp_example$lps[3], tolerance = 1e-12)
  expect_equal(thin(result$refs.risk), official$snp_example$refs_risk, tolerance = 1e-12)

  set.seed(1)
  windows <- do.call(computeAbsoluteRisk, c(snp_model_args(bc), list(
    apply.age.start = c(40, 50, 60), apply.age.interval.length = c(5, 10, 20),
    apply.cov.profile = bc$new_cov_prof[c(3, 3, 3), ],
    apply.snp.profile = bc$new_snp_prof[c(3, 3, 3), ], return.lp = TRUE,
    return.refs.risk = TRUE)))
  expect_equal(as.numeric(windows$risk), official$snp_windows$risk, tolerance = 1e-12)
  expect_equal(thin(windows$refs.risk), official$snp_windows$refs_risk, tolerance = 1e-12)
})

test_that("a SNP-only model reproduces earlier versions for genotyped profiles", {
  set.seed(1)
  result <- computeAbsoluteRisk(
    model.snp.info = bc$bc_72_snps, model.disease.incidence.rates = bc$bc_inc,
    model.competing.incidence.rates = bc$mort_inc, apply.age.start = 50,
    apply.age.interval.length = 30, apply.snp.profile = bc$new_snp_prof,
    return.refs.risk = TRUE)
  expect_equal(result$risk[3], official$snp_only$risk[3], tolerance = 1e-12)
  # one risk per simulated person of the first set (not averages over sets)
  expect_length(result$refs.risk, 10000L)
  expect_gt(stats::sd(result$refs.risk), 1.5 * stats::sd(official$snp_only$refs_risk))
  expect_message(
    population <- computeAbsoluteRisk(
      model.snp.info = bc$bc_72_snps, model.disease.incidence.rates = bc$bc_inc,
      apply.age.start = 50, apply.age.interval.length = 30),
    "10000 profiles")
  expect_length(population$risk, 10000L)
  expect_equal(stats::var(as.numeric(population$risk)), 0)
})

test_that("a single profile, extra columns and tibble-like inputs are accepted", {
  all_risks <- complete_risks(profiles_used = profiles[1:5, ], ages = inputs$ages[1:5],
                              lengths = inputs$lengths[1:5])
  one <- complete_risks(profiles_used = profiles[2, ], ages = inputs$ages[2],
                        lengths = inputs$lengths[2])
  expect_equal(one$risk[1], all_risks$risk[2])
  extra <- cbind(profiles[1:5, ], unused = letters[1:5])
  tibble_like <- structure(extra, class = c("tbl_df", "tbl", "data.frame"))
  expect_equal(complete_risks(profiles_used = tibble_like, ages = inputs$ages[1:5],
                              lengths = inputs$lengths[1:5])$risk,
               all_risks$risk)
  ages_as_matrix <- complete_risks(profiles_used = profiles[1:5, ],
                                   ages = matrix(inputs$ages[1:5]),
                                   lengths = inputs$lengths[1:5])
  expect_equal(ages_as_matrix$risk, all_risks$risk)
})

test_that("a single profile works with a model of continuous covariates only", {
  # (earlier versions stopped with "non-conformable arguments")
  info <- list(list(name = "famhist", type = "continuous"),
               list(name = "rd2_everhrt_e", type = "continuous"))
  args <- list(model.formula = diagnosis ~ famhist + rd2_everhrt_e,
               model.cov.info = info,
               model.log.RR = c(famhist = 0.6, rd2_everhrt_e = 0.2),
               model.ref.dataset = bc$ref_cov_dat[, c("famhist", "rd2_everhrt_e")],
               model.disease.incidence.rates = bc$bc_inc,
               model.competing.incidence.rates = bc$mort_inc,
               apply.age.start = 50, apply.age.interval.length = 10)
  profiles <- data.frame(famhist = c(0, 1, NA), rd2_everhrt_e = c(1, 0, 1))
  together <- do.call(computeAbsoluteRisk, c(args, list(apply.cov.profile = profiles)))
  for (i in 1:3) {
    alone <- do.call(computeAbsoluteRisk, c(args, list(apply.cov.profile = profiles[i, ])))
    expect_equal(alone$risk[1, 1], together$risk[i, 1])
  }
})

test_that("rates are needed only for the ages inside the intervals", {
  # The example mortality rates are NA from age 85. An interval ending at 85
  # does not include age 85, so its risk is defined (earlier versions gave NA,
  # or, when intervals differed, replaced every profile's risk by the average
  # of its percentile bin).
  expect_true(all(is.na(bc$mort_inc[bc$mort_inc[, 1] >= 85, 2])))
  profiles <- bc$ref_cov_dat[1:6, ]
  ending_at_85 <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = c(50, 60, 70, 50, 60, 70), apply.age.interval.length = c(35, 25, 15, 5, 5, 5),
    apply.cov.profile = profiles, return.lp = TRUE)))
  expect_false(anyNA(ending_at_85$risk))
  model <- do.call(build_risk_model, covariate_model_args(bc))
  prepared <- prepare_covariate_data(profiles, model$covariates$cov_info, "profile", TRUE)
  design <- build_design_matrix(model$covariates$spec, prepared)
  lp <- as.vector(design %*% model$beta)
  expect_equal(as.numeric(ending_at_85$lps), lp)
  expected <- vapply(1:6, function(i) {
    reference_absolute_risk(lp[i], ending_at_85$details$Int_Start[i],
                            ending_at_85$details$Int_End[i] - ending_at_85$details$Int_Start[i],
                            model$baseline, model$competing)[["risk"]]
  }, numeric(1))
  expect_equal(as.numeric(ending_at_85$risk), expected, tolerance = 1e-12)
  expect_error(do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 50, apply.age.interval.length = 36, apply.cov.profile = profiles))),
    "no rate for the age\\(s\\) 85")
})

test_that("the result holds the risks, the details and the calibration", {
  result <- complete_risks(profiles_used = profiles[1:4, ], ages = 40, lengths = 10)
  expect_s3_class(result, "icare")
  expect_named(result, c("risk", "details", "beta.used", "lps", "refs.risk"))
  expect_equal(dim(result$risk), c(4L, 1L))
  expect_equal(colnames(result$risk), "Risk_Estimate")
  expect_equal(names(result$details)[1:3], c("Int_Start", "Int_End", "Risk_Estimate"))
  expect_equal(result$details$Int_End, rep(50, 4))
  expect_equal(result$details$Risk_Estimate, as.numeric(result$risk))
  expect_equal(result$details[, -(1:3)], `rownames<-`(profiles[1:4, ], NULL))
  expect_length(result$refs.risk, nrow(bc$ref_cov_dat))
  calibration <- attr(result, "calibration")
  expect_named(calibration, c("age", "incidence.rate", "baseline.hazard", "expected.rr"))
  expect_equal(calibration$baseline.hazard * calibration$expected.rr,
               calibration$incidence.rate)
  plain <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 40, apply.age.interval.length = 10,
    apply.cov.profile = profiles[1:4, ])))
  expect_named(plain, c("risk", "details", "beta.used"))
})

test_that("the C-code switch is ignored with a warning, and unknown arguments are errors", {
  base <- complete_risks(profiles_used = profiles[1:3, ], ages = 40, lengths = 10)
  expect_warning(
    with_c <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
      apply.age.start = 40, apply.age.interval.length = 10,
      apply.cov.profile = profiles[1:3, ], use.c.code = 1))),
    "`use.c.code` is deprecated and ignored")
  expect_equal(with_c$risk, base$risk)
  expect_error(do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 40, apply.age.interval.length = 10,
    apply.cov.profile = profiles[1:3, ], use.C.code = 1))), "Unused argument")
})

test_that("invalid inputs give informative errors", {
  args <- c(covariate_model_args(bc), list(apply.cov.profile = profiles[1:3, ]))
  run <- function(...) do.call(computeAbsoluteRisk, override_args(args, list(...)))
  expect_error(run(apply.age.start = 80, apply.age.interval.length = 10),
               "has no rate for the age\\(s\\) 86-89")
  expect_error(run(apply.age.start = c(40, 50), apply.age.interval.length = 5),
               "apply.age.start")
  expect_error(run(apply.age.start = -1, apply.age.interval.length = 5), "negative")
  expect_error(run(apply.age.start = 40, apply.age.interval.length = 5, n.imp = 0), "n.imp")
  expect_error(do.call(computeAbsoluteRisk, override_args(args, list(
    apply.age.start = 40, apply.age.interval.length = 5,
    model.disease.incidence.rates = NULL))), "model.disease.incidence.rates")
  expect_error(computeAbsoluteRisk(model.disease.incidence.rates = bc$bc_inc,
                                   apply.age.start = 40, apply.age.interval.length = 5),
               "model.snp.info")
  wrong_level <- profiles[1:3, ]
  wrong_level$parity[2] <- 9
  expect_error(run(apply.age.start = 40, apply.age.interval.length = 5,
                   apply.cov.profile = wrong_level), "outside the levels")
})

test_that("SNP profiles without an SNP model are rejected", {
  bc <- load_example_data()
  expect_error(do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = 50, apply.age.interval.length = 30,
    apply.cov.profile = bc$new_cov_prof, apply.snp.profile = bc$new_snp_prof))),
    "the model has no SNPs")
})
