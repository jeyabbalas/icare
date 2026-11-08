# Absolute risk model: assembly and prediction
# ============================================
#
# build_risk_model() validates the model inputs once and assembles everything
# that does not depend on the profiles: the coefficients, the reference
# population (with simulated SNPs when needed), the rate tables and the
# calibrated baseline hazard. prepare_profiles() codes the profiles, and
# predict_risks() computes their risks over their windows. The public
# functions computeAbsoluteRisk(), computeAbsoluteRiskSplitInterval() and
# ModelValidation() are thin layers over these three steps.

#' Which kind of model the inputs describe
#'
#' The covariate part needs `model.formula`, `model.cov.info`, `model.log.RR`
#' and `model.ref.dataset` together; without all of them the model is
#' SNP-only and needs `model.snp.info`.
#' @return TRUE when the model has covariates.
#' @noRd
model_has_covariates <- function(formula, cov_info, log_rr, ref_dataset, snp_info) {
  supplied <- c(model.formula = !is.null(formula), model.cov.info = !is.null(cov_info),
                model.log.RR = !is.null(log_rr), model.ref.dataset = !is.null(ref_dataset))
  if (all(supplied)) {
    return(TRUE)
  }
  if (any(supplied)) {
    abort("A model with covariates needs all of `model.formula`, `model.cov.info`, ",
          "`model.log.RR` and `model.ref.dataset`; missing ",
          format_names(paste0("`", names(supplied)[!supplied], "`")),
          ". Leave all of them NULL for a SNP-only model.")
  }
  if (is.null(snp_info)) {
    abort("Specify a model: either `model.formula`, `model.cov.info`, ",
          "`model.log.RR` and `model.ref.dataset` (covariates), or ",
          "`model.snp.info` (SNPs), or both.")
  }
  FALSE
}

#' Reference weights normalized to sum to one
#' @noRd
normalize_reference_weights <- function(weights, n, name = "model.ref.dataset.weights") {
  if (is.null(weights)) {
    return(rep(1 / n, n))
  }
  if (is.data.frame(weights) || !is.numeric(as.vector(weights)) ||
      length(weights) != n) {
    abort("`", name, "` must be a numeric vector with one weight per row of ",
          "`model.ref.dataset` (", n, ").")
  }
  weights <- as.numeric(weights)
  if (anyNA(weights) || any(!is.finite(weights))) {
    abort("`", name, "` cannot contain missing or infinite values.")
  }
  if (any(weights < 0)) {
    abort("`", name, "` cannot contain negative values.")
  }
  if (sum(weights) <= 0) {
    abort("`", name, "` must have a positive sum.")
  }
  weights / sum(weights)
}

#' Validate the model inputs and calibrate the baseline hazard
#'
#' Simulating the SNP genotypes of the reference population is the only
#' random step; it draws exactly the random numbers earlier versions drew.
#' @noRd
build_risk_model <- function(model.formula = NULL, model.cov.info = NULL,
                             model.snp.info = NULL, model.log.RR = NULL,
                             model.ref.dataset = NULL, model.ref.dataset.weights = NULL,
                             model.disease.incidence.rates,
                             model.competing.incidence.rates = NULL,
                             model.bin.fh.name = NA, n.imp = 5) {
  has_covariates <- model_has_covariates(model.formula, model.cov.info, model.log.RR,
                                         model.ref.dataset, model.snp.info)
  has_snps <- !is.null(model.snp.info)
  n_imp <- as_count(n.imp, "n.imp", lower = 1, upper = 20)
  rates <- resolve_rate_tables(model.disease.incidence.rates,
                               model.competing.incidence.rates)

  model <- list(has_covariates = has_covariates, has_snps = has_snps, n_imp = n_imp,
                incidence = rates$incidence, competing = rates$competing,
                covariates = NULL, snps = NULL, family_history = NULL)

  # ----- Covariates
  beta_covariates <- numeric(0)
  if (has_covariates) {
    cov_info <- normalize_cov_info(model.cov.info, model.formula)
    reference_raw <- as_plain_data_frame(model.ref.dataset, "model.ref.dataset")
    reference <- prepare_covariate_data(reference_raw, cov_info, "model.ref.dataset",
                                        allow_missing = FALSE)
    if (nrow(reference) == 0L) {
      abort("`model.ref.dataset` has no rows.")
    }
    if (nrow(reference) < 200L) {
      warn("`model.ref.dataset` has only ", nrow(reference), " rows; the reference ",
           "sample should be large enough to represent the population.")
    }
    weights <- normalize_reference_weights(model.ref.dataset.weights, nrow(reference))
    setup <- build_design_spec(model.formula, reference)
    beta_covariates <- match_log_rr(model.log.RR, setup$spec$columns)
    model$covariates <- list(cov_info = cov_info, spec = setup$spec, data = reference,
                             design = setup$design, weights = weights)
  }

  # ----- SNPs and family history
  if (has_snps) {
    snps <- normalize_snp_info(model.snp.info)
    model$snps <- snps
    if (has_covariates && any(snps$names %in% names(model$covariates$cov_info))) {
      abort("SNP names in `model.snp.info` must differ from the covariate names.")
    }
    fh_name <- if (is_absent(model.bin.fh.name)) NULL else as.character(model.bin.fh.name)
    if (!is.null(fh_name) && has_covariates) {
      if (length(fh_name) != 1L || !fh_name %in% names(model$covariates$cov_info)) {
        abort("`model.bin.fh.name` must be the name of the family-history ",
              "covariate of the model (or NA when the model has none).")
      }
      fh_population <- as_family_history(reference_raw[[fh_name]], "model.ref.dataset",
                                         allow_missing = FALSE)
      fh_column <- family_history_column(model$covariates$spec, fh_name)
      attenuation <- family_history_attenuation(snps$log_or, snps$freq)
      beta_covariates[fh_column] <- beta_covariates[fh_column] - attenuation
      model$family_history <- list(name = fh_name, column = fh_column,
                                   attenuation = attenuation)
    } else if (has_covariates) {
      message("Note: SNP genotypes are imputed without regard to family history, ",
              "as `model.bin.fh.name` is NA.")
      fh_population <- rep(0, nrow(model$covariates$data))
    } else {
      fh_population <- rep(0, snp_only_population_size)
    }
  }

  # ----- Reference population and coefficients
  if (has_covariates && has_snps) {
    n_reference <- nrow(model$covariates$data)
    genotypes <- simulate_snp_genotypes(snps$log_or, snps$freq,
                                        rep(fh_population, n_imp))
    replicate_rows <- rep(seq_len(n_reference), n_imp)
    population <- list(
      design = cbind(genotypes, model$covariates$design[replicate_rows, , drop = FALSE]),
      weights = rep(model$covariates$weights, n_imp) / n_imp,
      reference_row = replicate_rows,
      block = rep(seq_len(n_imp), each = n_reference)
    )
    beta <- c(snps$log_or, beta_covariates)
  } else if (has_covariates) {
    n_reference <- nrow(model$covariates$data)
    population <- list(design = model$covariates$design,
                       weights = model$covariates$weights,
                       reference_row = seq_len(n_reference),
                       block = rep(1L, n_reference))
    beta <- beta_covariates
  } else {
    genotypes <- simulate_snp_genotypes(snps$log_or, snps$freq,
                                        rep(fh_population, n_imp))
    n_population <- nrow(genotypes)
    population <- list(design = genotypes,
                       weights = rep(1 / n_population, n_population),
                       reference_row = NULL,
                       block = rep(seq_len(n_imp), each = snp_only_population_size))
    beta <- snps$log_or
  }
  colnames(population$design) <- names(beta)
  population$lp <- as.vector(population$design %*% beta)
  model$population <- population
  model$beta <- beta
  model$column_variables <- c(if (has_snps) as.list(snps$names),
                              if (has_covariates) model$covariates$spec$column_covariates)

  # ----- Baseline hazard
  model$baseline <- calibrate_baseline_hazard(model$incidence, population$lp,
                                              population$weights)
  model
}

#' Code the profiles whose risks are requested
#'
#' @param n_default Number of profiles of a SNP-only model without
#'   `apply.snp.profile`.
#' @return A list with the profile design matrix (NA where a risk factor is
#'   missing), the missingness of each risk factor, the prepared covariate
#'   data, the genotypes and the original inputs (for the results).
#' @noRd
prepare_profiles <- function(model, apply.cov.profile = NULL, apply.snp.profile = NULL,
                             n_default = NULL) {
  out <- list(cov_profile = NULL, snp_profile = NULL, covariate_data = NULL)
  if (model$has_covariates) {
    if (is.null(apply.cov.profile)) {
      abort("`apply.cov.profile` must be supplied for a model with covariates.")
    }
    out$cov_profile <- as_plain_data_frame(apply.cov.profile, "apply.cov.profile")
    covariate_data <- prepare_covariate_data(out$cov_profile, model$covariates$cov_info,
                                             "apply.cov.profile", allow_missing = TRUE)
    if (!is.null(model$family_history)) {
      as_family_history(out$cov_profile[[model$family_history$name]],
                        "apply.cov.profile", allow_missing = TRUE)
    }
    out$covariate_data <- covariate_data
    n <- nrow(covariate_data)
    if (n == 0L) {
      abort("`apply.cov.profile` has no rows.")
    }
    cov_design <- build_design_matrix(model$covariates$spec, covariate_data)
    cov_missing <- is.na(as.matrix(covariate_data))
    colnames(cov_missing) <- names(covariate_data)
  } else {
    if (!is.null(apply.cov.profile)) {
      abort("`apply.cov.profile` was supplied, but the model has no covariates: ",
            "supply `model.formula`, `model.cov.info`, `model.log.RR` and ",
            "`model.ref.dataset` too.")
    }
    n <- NULL
    cov_design <- NULL
    cov_missing <- NULL
  }

  if (model$has_snps) {
    snp_names <- model$snps$names
    if (is.null(apply.snp.profile)) {
      if (is.null(n)) {
        n <- n_default
        message("Note: `apply.snp.profile` was not supplied; computing risks for ",
                n, " profiles with all SNPs missing (each is the population-average ",
                "risk). Supply `apply.snp.profile` for genotyped individuals.")
      } else {
        message("Note: `apply.snp.profile` was not supplied; all SNPs are imputed.")
      }
      genotypes <- matrix(NA_real_, nrow = n, ncol = length(snp_names),
                          dimnames = list(NULL, snp_names))
    } else {
      genotypes <- normalize_snp_profile(apply.snp.profile, snp_names, n)
      n <- nrow(genotypes)
      if (n == 0L) {
        abort("`apply.snp.profile` has no rows.")
      }
    }
    out$snp_profile <- genotypes
    snp_missing <- is.na(genotypes)
  } else {
    genotypes <- NULL
    snp_missing <- NULL
  }

  out$n <- n
  out$design <- cbind(genotypes, cov_design)
  colnames(out$design) <- names(model$beta)
  out$variable_missing <- cbind(snp_missing, cov_missing)
  out
}

#' Absolute risks of prepared profiles over their windows
#'
#' Profiles with every risk factor observed are evaluated directly; the others
#' are model-averaged over the reference population (missing_covariates.R).
#' @param survival Also return the event-free survival at the end of the window.
#' @return `list(risk, lp, survival, imputed)`; `imputed` flags the profiles
#'   with missing risk factors.
#' @noRd
predict_risks <- function(model, profiles, start, length, survival = FALSE) {
  needed <- exposed_ages(start, length)
  check_rate_coverage(model$incidence, needed, "model.disease.incidence.rates")
  check_rate_coverage(model$competing, needed, "model.competing.incidence.rates")

  n <- profiles$n
  risk <- lp <- surv <- rep(NA_real_, n)
  imputed <- rowSums(is.na(profiles$design)) > 0
  complete <- which(!imputed)
  if (length(complete) > 0L) {
    lp[complete] <- as.vector(profiles$design[complete, , drop = FALSE] %*% model$beta)
    res <- compute_absolute_risk_engine(lp[complete], start[complete], length[complete],
                                        model$baseline, model$competing, survival = TRUE)
    risk[complete] <- res$risk
    surv[complete] <- res$survival
  }
  incomplete <- which(imputed)
  if (length(incomplete) > 0L) {
    ctx <- list(design = profiles$design, variable_missing = profiles$variable_missing,
                column_variables = model$column_variables, population = model$population,
                beta = model$beta, start = start, length = length,
                baseline = model$baseline, competing = model$competing,
                survival = survival, covariate_data = profiles$covariate_data,
                reference_data = model$covariates$data,
                covariate_spec = model$covariates$spec)
    res <- impute_missing_risks(incomplete, ctx)
    risk[incomplete] <- res$risk
    lp[incomplete] <- res$lp
    surv[incomplete] <- res$survival
  }
  list(risk = risk, lp = lp, survival = surv, imputed = imputed)
}

#' Risks of the reference population over given windows
#'
#' For a model with covariates, one risk per row of `model.ref.dataset`; with
#' SNPs, each row's risk is averaged over its `n.imp` simulated genotype sets.
#' For a SNP-only model, the risks of the first set of simulated people.
#' @param start,length Windows: length 1, or one per reference row.
#' @noRd
reference_risks <- function(model, start, length) {
  population <- model$population
  first_block <- population$block == 1L
  n_reference <- sum(first_block)
  if (base::length(start) == 1L) start <- rep(start, n_reference)
  if (base::length(length) == 1L) length <- rep(length, n_reference)
  if (base::length(start) != n_reference || base::length(length) != n_reference) {
    abort("Internal error: reference windows must have length 1 or one per ",
          "reference row.")
  }
  needed <- exposed_ages(start, length)
  check_rate_coverage(model$incidence, needed, "model.disease.incidence.rates")
  check_rate_coverage(model$competing, needed, "model.competing.incidence.rates")

  if (!model$has_covariates) {
    return(compute_absolute_risk_engine(population$lp[first_block], start, length,
                                        model$baseline, model$competing))
  }
  row <- population$reference_row
  risk <- compute_absolute_risk_engine(population$lp, start[row], length[row],
                                       model$baseline, model$competing)
  as.vector(rowsum(risk, row, reorder = TRUE)) / model$n_imp^model$has_snps
}

#' Linear predictors of the reference rows (averaged over simulated SNPs)
#' @noRd
reference_linear_predictors <- function(model) {
  population <- model$population
  if (!model$has_covariates) {
    return(population$lp[population$block == 1L])
  }
  row <- population$reference_row
  as.vector(rowsum(population$lp, row, reorder = TRUE)) / model$n_imp^model$has_snps
}
