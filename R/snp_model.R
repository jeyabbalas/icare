# SNP part of the relative-risk model
# ===================================
#
# Each SNP k adds log(OR_k) * g_k to the linear predictor, where g_k in
# {0, 1, 2} counts the risk alleles (dosages in [0, 2] are accepted). SNPs are
# assumed independent and in Hardy-Weinberg equilibrium with risk-allele
# frequency f_k, so that in the population without a family history
#
#   P(g_k = 0, 1, 2) = ((1 - f_k)^2, 2 f_k (1 - f_k), f_k^2).
#
# For women with a family history of the disease (in a first-degree relative),
# genotypes are drawn from the distribution tilted by half the odds ratio per
# allele (Mendelian transmission, rare disease):
#
#   P(g_k = g | FH = 1)  proportional to  OR_k^(g / 2) P(g_k = g).
#
# Because the SNPs explain part of the familial risk, the log relative risk of
# a binary family-history covariate is attenuated by half the variance of the
# polygenic score: beta_FH - sum_k log(OR_k)^2 2 f_k (1 - f_k) / 2.
#
# The reference population used to calibrate the baseline hazard and to
# impute missing SNPs consists of `n.imp` copies of the reference covariate
# data, each with independently simulated genotypes (10,000 simulated people
# per copy in a SNP-only model).

# Number of simulated people per copy in a SNP-only model (as in earlier versions).
snp_only_population_size <- 10000L

#' Validate `model.snp.info`
#' @return `list(names, log_or, freq)`.
#' @noRd
normalize_snp_info <- function(snp_info) {
  snp_info <- tryCatch(as_plain_data_frame(snp_info, "model.snp.info"), error = function(e) {
    abort("`model.snp.info` must be a data frame with columns 'snp.name', ",
          "'snp.odds.ratio' and 'snp.freq'.")
  })
  required <- c("snp.name", "snp.odds.ratio", "snp.freq")
  absent <- setdiff(required, names(snp_info))
  if (length(absent) > 0L) {
    abort("`model.snp.info` must have columns 'snp.name', 'snp.odds.ratio' and ",
          "'snp.freq'; missing ", format_names(absent), ".")
  }
  if (nrow(snp_info) == 0L) {
    abort("`model.snp.info` has no rows.")
  }
  snp_names <- as.character(snp_info$snp.name)
  if (anyNA(snp_names) || any(!nzchar(snp_names)) || anyDuplicated(snp_names)) {
    abort("The 'snp.name' column of `model.snp.info` must hold distinct, ",
          "non-missing names.")
  }
  odds_ratio <- snp_info$snp.odds.ratio
  freq <- snp_info$snp.freq
  if (!is.numeric(odds_ratio) || anyNA(odds_ratio) || any(!is.finite(odds_ratio)) ||
      any(odds_ratio <= 0)) {
    abort("The 'snp.odds.ratio' column of `model.snp.info` must hold positive, ",
          "finite odds ratios.")
  }
  if (!is.numeric(freq) || anyNA(freq) || any(freq < 0 | freq > 1)) {
    abort("The 'snp.freq' column of `model.snp.info` must hold allele ",
          "frequencies between 0 and 1.")
  }
  list(names = snp_names,
       log_or = stats::setNames(log(as.numeric(odds_ratio)), snp_names),
       freq = stats::setNames(as.numeric(freq), snp_names))
}

#' Validate `apply.snp.profile`
#'
#' Columns are matched to `model.snp.info` by name when all SNP names are
#' present, and by position otherwise. Genotypes are risk-allele counts (or
#' dosages) between 0 and 2; NA marks a missing genotype.
#' @return A numeric matrix with one column per SNP, in `model.snp.info` order.
#' @noRd
normalize_snp_profile <- function(profile, snp_names, n = NULL) {
  if (is.null(dim(profile))) {
    profile <- matrix(profile, ncol = 1L)
  }
  profile <- as_plain_data_frame(profile, "apply.snp.profile")
  if (all(snp_names %in% names(profile))) {
    profile <- profile[, snp_names, drop = FALSE]
  } else if (ncol(profile) != length(snp_names)) {
    abort("`apply.snp.profile` must have one column per SNP in ",
          "`model.snp.info` (", length(snp_names), "); it has ", ncol(profile), ".")
  }
  numeric_columns <- vapply(profile, function(x) is.numeric(x) || all(is.na(x)), logical(1))
  if (!all(numeric_columns)) {
    abort("`apply.snp.profile` must hold numeric genotypes (0, 1, 2 or NA).")
  }
  genotypes <- matrix(as.numeric(unlist(profile, use.names = FALSE)),
                      nrow = nrow(profile), dimnames = list(NULL, snp_names))
  if (any(genotypes < 0 | genotypes > 2, na.rm = TRUE)) {
    abort("Genotypes in `apply.snp.profile` must be risk-allele counts between ",
          "0 and 2 (or NA when missing).")
  }
  if (!is.null(n) && nrow(genotypes) != n) {
    abort("`apply.cov.profile` and `apply.snp.profile` must have the same ",
          "number of rows.")
  }
  genotypes
}

#' Simulate genotypes given family history
#'
#' Reproduces earlier versions of iCARE draw for draw: one uniform matrix of
#' size n x K is drawn column by column and thresholded at the cumulative
#' genotype probabilities of each person's family-history group.
#'
#' @param log_or,freq Per-SNP log odds ratios and risk-allele frequencies.
#' @param family_history 0/1 vector, one element per simulated person.
#' @return An n x K matrix of risk-allele counts.
#' @noRd
simulate_snp_genotypes <- function(log_or, freq, family_history) {
  n <- length(family_history)
  k <- length(log_or)
  probs_no_fh <- cbind((1 - freq)^2, 2 * freq * (1 - freq), freq^2)
  tilted <- exp(outer(log_or, c(0, 1, 2) / 2)) * probs_no_fh
  probs_fh <- tilted / rowSums(tilted)

  uniforms <- matrix(stats::runif(k * n), nrow = n, ncol = k)
  genotypes <- matrix(NA_real_, nrow = n, ncol = k)
  draw <- function(rows, probs) {
    u <- uniforms[rows, , drop = FALSE]
    t(t(u) > probs[, 1]) + t(t(u) > probs[, 1] + probs[, 2])
  }
  no_fh <- which(family_history == 0)
  with_fh <- which(family_history == 1)
  if (length(no_fh) > 0L) genotypes[no_fh, ] <- draw(no_fh, probs_no_fh)
  if (length(with_fh) > 0L) genotypes[with_fh, ] <- draw(with_fh, probs_fh)
  genotypes
}

#' Attenuation of the family-history log relative risk explained by the SNPs
#' @noRd
family_history_attenuation <- function(log_or, freq) {
  sum(log_or^2 * 2 * freq * (1 - freq)) / 2
}

#' Binary family-history values coded 0/1
#' @noRd
as_family_history <- function(x, name, allow_missing) {
  values <- if (is.factor(x)) as.character(x) else x
  ok <- values %in% c(0, 1, "0", "1") | (allow_missing & is.na(values))
  if (!all(ok)) {
    abort("The family-history covariate must be binary (0/1",
          if (allow_missing) " or NA" else "", ") when SNPs are in the model; ",
          "check `", name, "`.")
  }
  as.numeric(as.character(values))
}

#' Design column of the family-history main effect
#' @noRd
family_history_column <- function(spec, fh_name) {
  candidates <- which(vapply(spec$column_covariates, function(v) identical(v, fh_name),
                             logical(1)))
  if (length(candidates) != 1L) {
    abort("`model.bin.fh.name` ('", fh_name, "') must enter `model.formula` as ",
          "a binary main effect with a single design column, so that its log ",
          "relative risk can be adjusted for the SNPs.")
  }
  spec$columns[candidates]
}
