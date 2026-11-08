bc <- load_example_data()

split_args <- function(...) {
  override_args(list(
    model.formula = bc$bc_model_formula, model.cov.info = bc$bc_model_cov_info,
    model.snp.info = bc$bc_72_snps, model.log.RR = bc$bc_model_log_or,
    model.log.RR.2 = bc$bc_model_log_or_post_50, model.ref.dataset = bc$ref_cov_dat,
    model.ref.dataset.2 = bc$ref_cov_dat_post_50,
    model.disease.incidence.rates = bc$bc_inc,
    model.competing.incidence.rates = bc$mort_inc, model.bin.fh.name = "famhist",
    apply.age.start = 30, apply.age.interval.length = 40, cut.time = 50,
    apply.cov.profile = bc$new_cov_prof, apply.snp.profile = bc$new_snp_prof
  ), list(...))
}

test_that("the two parts chain through the survival to the cut age", {
  set.seed(1)
  result <- do.call(computeAbsoluteRiskSplitInterval,
                    split_args(return.lp = TRUE, return.refs.risk = TRUE))
  set.seed(1)
  model_1 <- do.call(build_risk_model, snp_model_args(bc))
  model_2 <- do.call(build_risk_model, snp_model_args(
    bc, model.log.RR = bc$bc_model_log_or_post_50,
    model.ref.dataset = bc$ref_cov_dat_post_50))
  part_1 <- predict_risks(model_1, prepare_profiles(model_1, bc$new_cov_prof, bc$new_snp_prof),
                          rep(30, 3), rep(20, 3), survival = TRUE)
  part_2 <- predict_risks(model_2, prepare_profiles(model_2, bc$new_cov_prof, bc$new_snp_prof),
                          rep(50, 3), rep(20, 3))
  expect_equal(as.numeric(result$risk), part_1$risk + part_1$survival * part_2$risk,
               tolerance = 1e-12)
  expect_equal(as.numeric(result$lps.1), part_1$lp, tolerance = 1e-12)
  expect_equal(as.numeric(result$lps.2), part_2$lp, tolerance = 1e-12)
  expect_equal(result$refs.risk.1, reference_risks(model_1, 30, 20), tolerance = 1e-12)
  expect_equal(result$refs.risk.2, reference_risks(model_2, 50, 20), tolerance = 1e-12)
  # Earlier versions' unit test combined the same parts as r1 + (1 - r1) r2,
  # which overstates the risk of the genotyped profile 3.
  r1 <- part_1$risk[3]
  r2 <- part_2$risk[3]
  expect_equal(r1 + (1 - r1) * r2, 0.14428262, tolerance = 1e-7)
  expect_lt(result$risk[3], 0.14428262)
  expect_s3_class(result, "icareSplit")
  expect_equal(names(result$details)[1:3], c("Int_Start", "Cut_Age", "Int_End"))
  expect_equal(result$details$Cut_Age, rep(50, 3))
  expect_equal(result$details$Risk_Estimate, as.numeric(result$risk))
})

test_that("the same model on both sides reproduces the risk over the whole interval", {
  profiles <- bc$ref_cov_dat[c(5, 500, 5000, 9000), ]
  profiles$bmi_dec[4] <- NA
  args <- c(covariate_model_args(bc), list(
    apply.cov.profile = profiles, apply.age.start = c(40, 45.5, 30, 42),
    apply.age.interval.length = c(20, 10, 40, 25)))
  whole <- as.numeric(do.call(computeAbsoluteRisk, args)$risk)
  split <- do.call(computeAbsoluteRiskSplitInterval,
                   c(args, list(cut.time = c(50, 47, 70, 50))))
  expect_equal(as.numeric(split$risk)[1:3], whole[1:3], tolerance = 1e-12)
  # with a missing covariate the parts are averaged separately
  expect_equal(split$risk[4], whole[4], tolerance = 1e-3)
  expect_equal(split$details$Cut_Age, c(50, 47, 70, 50))
  # a cut outside the interval leaves one model
  outside <- do.call(computeAbsoluteRiskSplitInterval,
                     c(args, list(cut.time = c(20, 90, 20, 90))))
  expect_equal(as.numeric(outside$risk)[1:3], whole[1:3], tolerance = 1e-12)
  expect_equal(outside$details$Cut_Age, c(40, 55.5, 30, 67))

  # the SNPs of a shared model are simulated once
  set.seed(1)
  whole_snp <- do.call(computeAbsoluteRisk, c(snp_model_args(bc), list(
    apply.age.start = 30, apply.age.interval.length = 40,
    apply.cov.profile = bc$new_cov_prof, apply.snp.profile = bc$new_snp_prof)))
  set.seed(1)
  split_snp <- do.call(computeAbsoluteRiskSplitInterval,
                       split_args(model.log.RR.2 = NULL, model.ref.dataset.2 = NULL))
  expect_equal(split_snp$risk[3], whole_snp$risk[3], tolerance = 1e-12)
})

test_that("the second part reuses the reference weights (fork patch 3)", {
  profiles <- bc$ref_cov_dat[c(5, 500, 5000), ]
  set.seed(2)
  counts <- sample(1:3, nrow(bc$ref_cov_dat), replace = TRUE)
  run <- function(...) {
    do.call(computeAbsoluteRiskSplitInterval, c(covariate_model_args(bc, ...), list(
      model.log.RR.2 = bc$bc_model_log_or_post_50, apply.cov.profile = profiles,
      apply.age.start = 40, apply.age.interval.length = 30, cut.time = 50)))
  }
  default <- run(model.ref.dataset.weights = counts)
  explicit <- run(model.ref.dataset.weights = counts,
                  model.ref.dataset.2 = bc$ref_cov_dat, model.ref.dataset.weights.2 = counts)
  expect_equal(default$risk, explicit$risk)
  replicated <- run(model.ref.dataset = bc$ref_cov_dat[rep(seq_along(counts), counts), ])
  expect_equal(default$risk, replicated$risk, tolerance = 1e-10)
  uniform <- run()
  expect_false(isTRUE(all.equal(default$risk, uniform$risk)))
})

test_that("without a cut age and a second model the risk is over the whole interval", {
  profiles <- bc$ref_cov_dat[1:3, ]
  args <- c(covariate_model_args(bc), list(apply.cov.profile = profiles,
                                           apply.age.start = 40,
                                           apply.age.interval.length = 20,
                                           return.refs.risk = TRUE))
  whole <- do.call(computeAbsoluteRisk, args)
  single <- do.call(computeAbsoluteRiskSplitInterval, args)
  expect_s3_class(single, "icareSplit")
  expect_equal(single$risk, whole$risk)
  expect_equal(single$refs.risk, whole$refs.risk)
})

test_that("a second covariate profile is reported with suffixed names", {
  profiles <- bc$ref_cov_dat[1:2, ]
  later <- profiles
  later$bmi_dec <- c(10, 10)
  result <- do.call(computeAbsoluteRiskSplitInterval, c(covariate_model_args(bc), list(
    apply.cov.profile = profiles, apply.cov.profile.2 = later, apply.age.start = 40,
    apply.age.interval.length = 30, cut.time = 50)))
  expect_true(all(paste0(names(profiles), ".2") %in% names(result$details)))
  expect_equal(result$details$bmi_dec.2, c(10, 10))
  expect_equal(names(result$details)[ncol(result$details)], "Risk_Estimate")
})

test_that("invalid cut ages and second models are rejected", {
  args <- c(covariate_model_args(bc), list(
    apply.cov.profile = bc$ref_cov_dat[1:2, ], apply.age.start = 40,
    apply.age.interval.length = 30))
  run <- function(...) do.call(computeAbsoluteRiskSplitInterval, c(args, list(...)))
  expect_error(run(model.log.RR.2 = bc$bc_model_log_or_post_50), "`cut.time` must be supplied")
  expect_error(run(model.ref.dataset.weights.2 = rep(1, nrow(bc$ref_cov_dat))),
               "`cut.time` must be supplied")
  expect_error(run(model.bin.fh.name.2 = "famhist"), "`cut.time` must be supplied")
  expect_equal(run(model.bin.fh.name.2 = NA)$risk, run()$risk)
  expect_error(run(cut.time = NA), "`cut.time` must be numeric")
  expect_error(run(cut.time = c(50, 55, 60)), "length 1 or one value per profile")
  expect_error(run(cut.time = 50, apply.cov.profile.2 = bc$ref_cov_dat[1:3, ]),
               "as many rows")
  expect_warning(run(cut.time = 50, use.c.code = 1), "deprecated")
})
