bc <- load_example_data()

test_that("percentile bins follow earlier versions", {
  lp <- c(0.123456, 0.2, 0.3, 0.4, 0.5)
  # sd 0.148, so the cut points are rounded down to 3 decimals
  cuts <- get_cutpoints(lp, seq(0, 1, 0.25))
  expect_equal(cuts, c(floor(c(0.123456, 0.2, 0.3, 0.4, 0.5) * 1000) / 1000, Inf))
  expect_equal(get_miss_perc(cuts, c(0.1, 0.123, 0.25, 0.5, 7)), c(1, 1, 2, 5, 5))
  expect_equal(get_cutpoints(c(1, 1, 1, 2), seq(0, 1, 0.25)), c(1, 1.25, 2, Inf))
  widen <- function(threshold) function(low, high) high - low > threshold
  expect_equal(get_new_endpoints(cuts, 2, widen(0.25)), cuts[c(1, 4)])
  expect_equal(get_new_endpoints(cuts, 2, widen(0.3)), cuts[c(1, 5)])
  expect_equal(get_new_endpoints(cuts, 2, widen(10)), cuts[c(1, 6)])
})

# Profiles of the example model with missing covariates: a main effect only,
# one side of the interaction, both sides of the interaction, a continuous
# covariate, several factors, everything, and nothing.
missing_profiles <- function(bc) {
  profiles <- bc$ref_cov_dat[c(11, 250, 3000, 4100, 5200, 6300, 7400, 8500), ]
  profiles$parity[1] <- NA
  profiles$bmi_dec[2] <- NA
  profiles[3, c("rd_menohrt", "bmi_dec")] <- NA
  profiles$famhist[4] <- NA
  profiles[5, c("menarche_dec", "height_dec", "alcoholdweek_dec")] <- NA
  profiles$rd_menohrt[6] <- NA
  profiles[7, ] <- NA
  rownames(profiles) <- NULL
  profiles
}
starts <- c(35, 42.5, 50, 61.3, 45, 55, 40, 30)
lengths <- c(10, 7.5, 20, 5, 15, 3, 30, 25)

test_that("missing covariates are completed with the donors' values (complete-data estimand)", {
  profiles <- missing_profiles(bc)
  for (weights in list(NULL, scenario_inputs(bc)$weights)) {
    args <- covariate_model_args(bc, model.ref.dataset.weights = weights)
    result <- do.call(computeAbsoluteRisk, c(args, list(
      apply.age.start = starts, apply.age.interval.length = lengths,
      apply.cov.profile = profiles, return.lp = TRUE)))
    model <- do.call(build_risk_model, args)
    prepared <- prepare_covariate_data(profiles, model$covariates$cov_info, "profile",
                                       allow_missing = TRUE)
    for (i in 1:7) {
      expected <- reference_missing_covariate_risk(model, bc$bc_model_formula,
                                                   prepared[i, ], starts[i], lengths[i])
      expect_equal(result$risk[i], expected[["risk"]], tolerance = 1e-12,
                   label = paste("risk of profile", i))
      expect_equal(result$lps[i], expected[["lp"]], tolerance = 1e-12,
                   label = paste("linear predictor of profile", i))
    }
  }
})

test_that("interactions with a continuous covariate and of three factors are completed exactly", {
  # famhist is continuous, so its interaction is rebuilt from the completed
  # data; the interactions of factors are looked up in tables by level.
  formula <- diagnosis ~ famhist * as.factor(bmi_dec) + as.factor(parity) + rd2_everhrt_e +
    as.factor(rd_menohrt) * as.factor(ever_smoke) * as.factor(alcoholdweek_dec)
  variables <- all.vars(formula)[-1]
  # (one element of the example spells the field "names", matched partially)
  info_names <- vapply(bc$bc_model_cov_info, function(x) x[["name", exact = FALSE]],
                       character(1))
  cov_info <- bc$bc_model_cov_info[match(variables, info_names)]
  reference <- bc$ref_cov_dat[, variables]
  prepared <- prepare_covariate_data(reference, normalize_cov_info(cov_info, formula),
                                     "reference", allow_missing = FALSE)
  columns <- build_design_spec(formula, prepared)$spec$columns
  set.seed(3)
  log_rr <- stats::setNames(stats::rnorm(length(columns), sd = 0.3), columns)
  args <- list(model.formula = formula, model.cov.info = cov_info, model.log.RR = log_rr,
               model.ref.dataset = reference, model.disease.incidence.rates = bc$bc_inc,
               model.competing.incidence.rates = bc$mort_inc)

  profiles <- reference[c(15, 260, 3010, 4120, 5230, 6340, 7450, 8560), ]
  profiles$famhist[1] <- NA                                   # the continuous side
  profiles$bmi_dec[2] <- NA                                   # the factor side
  profiles$ever_smoke[3] <- NA                                # three interactions at once
  profiles[4, c("rd_menohrt", "famhist")] <- NA               # a table and a rebuilt block
  profiles[5, c("rd_menohrt", "ever_smoke", "alcoholdweek_dec")] <- NA  # a whole interaction
  profiles[6, c("ever_smoke", "alcoholdweek_dec")] <- NA      # two of its three factors
  profiles$parity[7] <- NA                                    # a main effect only
  rownames(profiles) <- NULL                                  # (profile 8 is complete)
  start <- c(35, 42.5, 50, 61.3, 45, 55, 40, 30)
  length <- c(10, 7.5, 20, 5, 15, 3, 30, 25)

  result <- do.call(computeAbsoluteRisk, c(args, list(
    apply.age.start = start, apply.age.interval.length = length,
    apply.cov.profile = profiles, return.lp = TRUE)))
  model <- do.call(build_risk_model, args)
  expect_equal(lengths(model$population$interaction_tables) > 0,
               lengths(model$population$block_covariates) >= 2 &
                 !vapply(model$population$block_covariates,
                         function(v) "famhist" %in% v, logical(1)))
  profile_data <- prepare_covariate_data(profiles, model$covariates$cov_info, "profile",
                                         allow_missing = TRUE)
  for (i in 1:7) {
    expected <- reference_missing_covariate_risk(model, formula, profile_data[i, ],
                                                 start[i], length[i])
    expect_equal(result$risk[i], expected[["risk"]], tolerance = 1e-12,
                 label = paste("risk of profile", i))
    expect_equal(result$lps[i], expected[["lp"]], tolerance = 1e-12,
                 label = paste("linear predictor of profile", i))
  }
})

test_that("the interaction tables hold each block's contribution by level", {
  model <- do.call(build_risk_model, covariate_model_args(bc))
  population <- model$population
  tabulated <- which(lengths(population$interaction_tables) > 0)
  expect_identical(lapply(population$block_covariates[tabulated], sort),
                   list(c("bmi_dec", "rd_menohrt")))
  for (b in tabulated) {
    codes <- vapply(population$block_covariates[[b]], function(v) {
      as.integer(model$covariates$data[[v]])
    }, integer(nrow(model$covariates$data)))
    expect_equal(population$interaction_tables[[b]][codes], population$design_lp[[b]],
                 tolerance = 1e-15)
  }
})

test_that("a profile with no observed covariate gets the population risk over its window", {
  profiles <- missing_profiles(bc)[c(7, 7, 7), ]
  weights <- scenario_inputs(bc)$weights
  result <- do.call(computeAbsoluteRisk, c(
    covariate_model_args(bc, model.ref.dataset.weights = weights),
    list(apply.age.start = c(45, 60, 30), apply.age.interval.length = c(10, 5, 50),
         apply.cov.profile = profiles, return.refs.risk = TRUE)))
  expect_equal(result$risk[1], stats::weighted.mean(result$refs.risk, weights),
               tolerance = 1e-12)
  model <- do.call(build_risk_model, covariate_model_args(bc, model.ref.dataset.weights = weights))
  expect_equal(result$risk[2:3],
               c(stats::weighted.mean(reference_risks(model, 60, 5), weights),
                 stats::weighted.mean(reference_risks(model, 30, 50), weights)),
               tolerance = 1e-12)
})

test_that("a profile's risk does not depend on the other profiles (fork patch 4)", {
  profiles <- missing_profiles(bc)
  together <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = starts, apply.age.interval.length = lengths,
    apply.cov.profile = profiles)))
  alone <- vapply(seq_len(nrow(profiles)), function(i) {
    as.numeric(do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
      apply.age.start = starts[i], apply.age.interval.length = lengths[i],
      apply.cov.profile = profiles[i, ])))$risk)
  }, numeric(1))
  expect_equal(as.numeric(together$risk), alone, tolerance = 1e-14)
  shuffled <- c(5, 2, 8, 1, 7, 3, 6, 4)
  permuted <- do.call(computeAbsoluteRisk, c(covariate_model_args(bc), list(
    apply.age.start = starts[shuffled], apply.age.interval.length = lengths[shuffled],
    apply.cov.profile = profiles[shuffled, ])))
  expect_equal(as.numeric(permuted$risk), alone[shuffled], tolerance = 1e-14)
})

test_that("missing SNPs are completed with the donors' simulated genotypes", {
  set.seed(1)
  result <- do.call(computeAbsoluteRisk, c(snp_model_args(bc), list(
    apply.age.start = 50, apply.age.interval.length = 30,
    apply.cov.profile = bc$new_cov_prof, apply.snp.profile = bc$new_snp_prof,
    return.lp = TRUE)))
  set.seed(1)
  model <- suppressMessages(do.call(build_risk_model, snp_model_args(bc)))
  profiles <- prepare_profiles(model, bc$new_cov_prof, bc$new_snp_prof)
  for (i in 1:2) {
    expected <- reference_missing_column_risk(model, profiles$design[i, ], 50, 30)
    expect_equal(result$risk[i], expected[["risk"]], tolerance = 1e-12)
    expect_equal(result$lps[i], expected[["lp"]], tolerance = 1e-12)
  }
  # the values reported for this example in the changes from earlier versions
  expect_equal(as.numeric(result$risk), c(0.10227223, 0.09002969, 0.16893543),
               tolerance = 1e-7)
})
