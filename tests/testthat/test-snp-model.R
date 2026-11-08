bc <- load_example_data()

test_that("SNP information is validated", {
  snps <- normalize_snp_info(bc$bc_72_snps)
  expect_equal(unname(snps$log_or), log(bc$bc_72_snps$snp.odds.ratio))
  expect_error(normalize_snp_info(bc$bc_72_snps[, 1:2]), "snp.freq")
  bad <- bc$bc_72_snps
  bad$snp.freq[1] <- 1.5
  expect_error(normalize_snp_info(bad), "between 0 and 1")
  bad <- bc$bc_72_snps
  bad$snp.odds.ratio[1] <- 0
  expect_error(normalize_snp_info(bad), "positive")
  bad <- bc$bc_72_snps
  bad$snp.name[2] <- bad$snp.name[1]
  expect_error(normalize_snp_info(bad), "distinct")
})

test_that("SNP profiles are matched by name, or by position when unnamed", {
  names <- c("s1", "s2")
  by_name <- normalize_snp_profile(data.frame(s2 = c(0, 1), s1 = c(2, NA)), names)
  expect_equal(by_name, matrix(c(2, NA, 0, 1), 2, dimnames = list(NULL, names)))
  by_position <- normalize_snp_profile(matrix(c(2, NA, 0, 1), 2), names)
  expect_equal(by_position, by_name)
  expect_error(normalize_snp_profile(data.frame(a = 1), names), "one column per SNP")
  # A misspelled name must not silently switch to matching by position.
  expect_error(normalize_snp_profile(data.frame(s2 = c(0, 1), s_1 = c(2, NA)), names),
               "none for s1")
  expect_error(normalize_snp_profile(data.frame(s1 = 3, s2 = 0), names), "between 0 and 2")
  expect_error(normalize_snp_profile(data.frame(s1 = 1, s2 = 0), names, n = 2),
               "same number of rows")
})

test_that("simulated genotypes follow Hardy-Weinberg and the family-history tilt", {
  set.seed(11)
  log_or <- log(c(1.5, 0.8))
  freq <- c(0.3, 0.6)
  n <- 200000
  no_fh <- simulate_snp_genotypes(log_or, freq, rep(0, n))
  expect_equal(colMeans(no_fh), 2 * freq, tolerance = 0.01)
  expect_equal(mean(no_fh[, 1] == 1), 2 * 0.3 * 0.7, tolerance = 0.01)
  with_fh <- simulate_snp_genotypes(log_or, freq, rep(1, n))
  tilted <- exp(outer(log_or, (0:2) / 2)) *
    cbind((1 - freq)^2, 2 * freq * (1 - freq), freq^2)
  tilted <- tilted / rowSums(tilted)
  expect_equal(colMeans(with_fh), as.vector(tilted %*% 0:2), tolerance = 0.01)
})

test_that("simulation draws one uniform matrix, column by column", {
  set.seed(3)
  genotypes <- simulate_snp_genotypes(log(c(1.2, 1.1, 0.9)), c(0.2, 0.5, 0.7),
                                      c(0, 1, 0, 1))
  set.seed(3)
  u <- matrix(stats::runif(12), nrow = 4)
  probs <- c(0.2, 0.5, 0.7)
  expected_no_fh <- t(t(u) > (1 - probs)^2) + t(t(u) > 1 - probs^2)
  expect_equal(genotypes[c(1, 3), ], expected_no_fh[c(1, 3), ])
})

test_that("the family-history log relative risk is attenuated by half the PRS variance", {
  log_or <- log(c(1.2, 0.9))
  freq <- c(0.3, 0.4)
  expect_equal(family_history_attenuation(log_or, freq),
               sum(log_or^2 * 2 * freq * (1 - freq)) / 2)
  expect_equal(as_family_history(c(0, 1, NA), "x", allow_missing = TRUE), c(0, 1, NA))
  expect_error(as_family_history(c(0, 2), "x", allow_missing = FALSE), "binary")
})

test_that("the attenuation applies to family history whatever its coding", {
  bc <- load_example_data()
  # Family history as a factor with reference level "1": the design column
  # famhist0 indicates its absence, with the opposite log relative risk.
  cov_info <- bc$bc_model_cov_info
  fh <- which(vapply(cov_info, function(x) x[["name", exact = FALSE]], "") == "famhist")
  cov_info[[fh]] <- list(name = "famhist", type = "factor", levels = c("0", "1"), ref = "1")
  log_rr <- bc$bc_model_log_or
  names(log_rr)[names(log_rr) == "famhist"] <- "famhist0"
  log_rr["famhist0"] <- -log_rr["famhist0"]
  profiles <- bc$new_cov_prof
  profiles$famhist <- c(1, 0, 1)
  risks <- function(...) {
    set.seed(50)
    suppressMessages(do.call(computeAbsoluteRisk, c(snp_model_args(bc, ...), list(
      apply.age.start = 50, apply.age.interval.length = 30,
      apply.cov.profile = profiles, apply.snp.profile = bc$new_snp_prof[c(3, 3, 3), ]))))
  }
  usual <- risks()
  reversed <- risks(model.cov.info = cov_info, model.log.RR = log_rr)
  attenuation <- family_history_attenuation(log(bc$bc_72_snps$snp.odds.ratio),
                                            bc$bc_72_snps$snp.freq)
  expect_equal(usual$beta.used["famhist", 1],
               unname(bc$bc_model_log_or["famhist"]) - attenuation)
  expect_equal(reversed$beta.used["famhist0", 1], -usual$beta.used["famhist", 1])
  expect_equal(reversed$risk, usual$risk)
})
