#' Absolute risk with a model that changes at a cut age
#'
#' @description
#' Estimates absolute risks over an age interval when the relative-risk model,
#' the covariate profiles or the reference population change at a cut age
#' (for example, different log relative risks before and after menopause).
#' The interval \eqn{[a, a + \tau)} is split at the cut age \eqn{c}; the first
#' model gives the risk \eqn{AR_1} over \eqn{[a, c)} and the probability
#' \eqn{S_1} of reaching age \eqn{c} alive and free of the disease, the second
#' model gives the risk \eqn{AR_2} over \eqn{[c, a + \tau)}, and the risk over
#' the whole interval is
#' \deqn{AR = AR_1 + S_1 AR_2.}
#' With the same model on both sides and a whole-number cut age, this equals
#' the risk from [computeAbsoluteRisk()] over the whole interval for profiles
#' without missing risk factors. (A fractional cut age splits a year of age,
#' within which the yearly risk model is not exactly additive, so the two can
#' then differ very slightly.) (Earlier versions combined the parts as
#' \eqn{AR_1 + (1 - AR_1) AR_2}, which ignores deaths from competing causes
#' before the cut age and overstates the risk.)
#'
#' @details
#' Arguments ending in `.2` describe the model after the cut age; each one
#' left `NULL` is taken from its counterpart before the cut age (the
#' reference weights are reused with the reference dataset). The disease and
#' competing incidence rates, the SNP information and `n.imp` are shared. If
#' `cut.time` is outside a profile's interval, that profile's risk comes from
#' one model only. If no `.2` argument and no `cut.time` is given, this is
#' [computeAbsoluteRisk()] over the whole interval.
#'
#' For profiles with missing risk factors, \eqn{AR_1}, \eqn{S_1} and
#' \eqn{AR_2} are each averaged over the reference population of their model
#' (the two models may have different reference datasets, so the missing
#' values cannot be completed jointly). The result can therefore differ
#' slightly from the risk over the whole interval even when both models are
#' the same; in the bundled example the relative difference is below 0.1%.
#'
#' @inheritParams computeAbsoluteRisk
#' @param model.bin.fh.name Name of the binary family-history covariate (see
#'   [computeAbsoluteRisk()]); `NULL` or `NA` when not applicable.
#' @param cut.time Cut age(s) at which the model changes: a single value or
#'   one per profile.
#' @param apply.cov.profile.2,model.formula.2,model.log.RR.2,model.ref.dataset.2,model.ref.dataset.weights.2,model.cov.info.2,model.bin.fh.name.2
#'   The model after the cut age (see Details).
#' @param return.refs.risk If `TRUE`, also return the risks of each model's
#'   reference population over the first profile's part of the interval.
#'
#' @return An object of class `"icareSplit"`, a list with
#' * `risk`: one-column matrix (`Risk_Estimate`) of absolute risks;
#' * `details`: data frame with the start of the interval (`Int_Start`), the
#'   age at which the second model starts (`Cut_Age`), the end of the
#'   interval (`Int_End`), the covariate profile and the risk;
#' * `lps.1`, `lps.2`: linear predictors of both models (if `return.lp`);
#' * `refs.risk.1`, `refs.risk.2`: reference risks of both models (if
#'   `return.refs.risk`).
#'
#' Without a cut age and second model, the result of [computeAbsoluteRisk()]
#' with class `"icareSplit"`.
#'
#' @seealso [computeAbsoluteRisk()]
#'
#' @examples
#' data(bc_data, package = "iCARE")
#'
#' # Risk from age 40 to 80 with different log relative risks after age 50
#' results <- computeAbsoluteRiskSplitInterval(
#'   model.formula = bc_model_formula,
#'   model.cov.info = bc_model_cov_info,
#'   model.log.RR = bc_model_log_or,
#'   model.log.RR.2 = bc_model_log_or_post_50,
#'   model.ref.dataset = ref_cov_dat,
#'   model.ref.dataset.2 = ref_cov_dat_post_50,
#'   model.disease.incidence.rates = bc_inc,
#'   model.competing.incidence.rates = mort_inc,
#'   apply.age.start = 40,
#'   apply.age.interval.length = 40,
#'   cut.time = 50,
#'   apply.cov.profile = new_cov_prof
#' )
#' results
#' @export
computeAbsoluteRiskSplitInterval <- function(apply.age.start, apply.age.interval.length,
                                             apply.cov.profile, model.formula,
                                             model.disease.incidence.rates, model.log.RR,
                                             model.ref.dataset,
                                             model.ref.dataset.weights = NULL,
                                             model.cov.info,
                                             model.competing.incidence.rates = NULL,
                                             return.lp = FALSE, apply.snp.profile = NULL,
                                             model.snp.info = NULL, model.bin.fh.name = NULL,
                                             cut.time = NULL, apply.cov.profile.2 = NULL,
                                             model.formula.2 = NULL, model.log.RR.2 = NULL,
                                             model.ref.dataset.2 = NULL,
                                             model.ref.dataset.weights.2 = NULL,
                                             model.cov.info.2 = NULL,
                                             model.bin.fh.name.2 = NULL, n.imp = 5,
                                             return.refs.risk = FALSE, ...) {
  check_removed_arguments(list(...), "computeAbsoluteRiskSplitInterval")
  # Covariate arguments may be omitted together for a SNP-only model.
  if (missing(apply.cov.profile)) apply.cov.profile <- NULL
  if (missing(model.formula)) model.formula <- NULL
  if (missing(model.log.RR)) model.log.RR <- NULL
  if (missing(model.ref.dataset)) model.ref.dataset <- NULL
  if (missing(model.cov.info)) model.cov.info <- NULL
  return.lp <- as_flag(return.lp, "return.lp")
  return.refs.risk <- as_flag(return.refs.risk, "return.refs.risk")
  if (is.null(model.bin.fh.name)) model.bin.fh.name <- NA

  second_model <- list(apply.cov.profile.2, model.formula.2, model.log.RR.2,
                       model.ref.dataset.2, model.ref.dataset.weights.2, model.cov.info.2,
                       if (!is_absent(model.bin.fh.name.2)) model.bin.fh.name.2)
  if (is.null(cut.time) && all(vapply(second_model, is.null, logical(1)))) {
    result <- computeAbsoluteRisk(
      model.formula = model.formula, model.cov.info = model.cov.info,
      model.snp.info = model.snp.info, model.log.RR = model.log.RR,
      model.ref.dataset = model.ref.dataset,
      model.ref.dataset.weights = model.ref.dataset.weights,
      model.disease.incidence.rates = model.disease.incidence.rates,
      model.competing.incidence.rates = model.competing.incidence.rates,
      model.bin.fh.name = model.bin.fh.name, n.imp = n.imp,
      apply.age.start = apply.age.start,
      apply.age.interval.length = apply.age.interval.length,
      apply.cov.profile = apply.cov.profile, apply.snp.profile = apply.snp.profile,
      return.lp = return.lp, return.refs.risk = return.refs.risk
    )
    class(result) <- "icareSplit"
    return(result)
  }
  if (is.null(cut.time)) {
    abort("`cut.time` must be supplied to use different models before and after ",
          "a cut age.")
  }

  # ----- The model after the cut age defaults to the model before it
  profile_2_supplied <- !is.null(apply.cov.profile.2)
  if (is.null(model.ref.dataset.2)) {
    model.ref.dataset.2 <- model.ref.dataset
    if (is.null(model.ref.dataset.weights.2)) {
      model.ref.dataset.weights.2 <- model.ref.dataset.weights
    }
  }
  apply.cov.profile.2 <- apply.cov.profile.2 %||% apply.cov.profile
  model.formula.2 <- model.formula.2 %||% model.formula
  model.log.RR.2 <- model.log.RR.2 %||% model.log.RR
  model.cov.info.2 <- model.cov.info.2 %||% model.cov.info
  if (is.null(model.bin.fh.name.2)) model.bin.fh.name.2 <- model.bin.fh.name

  build <- function(formula, cov_info, log_rr, ref, weights, fh) {
    build_risk_model(
      model.formula = formula, model.cov.info = cov_info,
      model.snp.info = model.snp.info, model.log.RR = log_rr,
      model.ref.dataset = ref, model.ref.dataset.weights = weights,
      model.disease.incidence.rates = model.disease.incidence.rates,
      model.competing.incidence.rates = model.competing.incidence.rates,
      model.bin.fh.name = fh, n.imp = n.imp
    )
  }
  model_1 <- build(model.formula, model.cov.info, model.log.RR, model.ref.dataset,
                   model.ref.dataset.weights, model.bin.fh.name)
  same_model <- identical(model.formula.2, model.formula) &&
    identical(model.cov.info.2, model.cov.info) &&
    identical(model.log.RR.2, model.log.RR) &&
    identical(model.ref.dataset.2, model.ref.dataset) &&
    identical(model.ref.dataset.weights.2, model.ref.dataset.weights) &&
    identical(model.bin.fh.name.2, model.bin.fh.name)
  # The same model on both sides is built (and its SNPs simulated) once, so
  # that splitting an interval reproduces the risk over the whole interval.
  model_2 <- if (same_model) {
    model_1
  } else {
    build(model.formula.2, model.cov.info.2, model.log.RR.2, model.ref.dataset.2,
          model.ref.dataset.weights.2, model.bin.fh.name.2)
  }

  # ----- Profiles and the two parts of each interval
  n_default <- default_profile_count(apply.age.start, apply.age.interval.length)
  profiles_1 <- prepare_profiles(model_1, apply.cov.profile, apply.snp.profile, n_default)
  profiles_2 <- prepare_profiles(model_2, apply.cov.profile.2, apply.snp.profile, n_default)
  if (profiles_2$n != profiles_1$n) {
    abort("`apply.cov.profile.2` must have as many rows as `apply.cov.profile`.")
  }
  profile_name <- if (model_1$has_covariates) "apply.cov.profile" else "apply.snp.profile"
  windows <- normalize_risk_windows(apply.age.start, apply.age.interval.length,
                                    profiles_1$n, profile_name)
  cut <- normalize_cut_ages(cut.time, profiles_1$n)
  length_1 <- pmin(pmax(cut - windows$start, 0), windows$length)
  start_2 <- windows$start + length_1
  length_2 <- windows$start + windows$length - start_2

  part_1 <- predict_risks(model_1, profiles_1, windows$start, length_1, survival = TRUE)
  part_2 <- predict_risks(model_2, profiles_2, start_2, length_2)
  risk <- part_1$risk + part_1$survival * part_2$risk

  # ----- Results
  details <- data.frame(Int_Start = windows$start, Cut_Age = start_2,
                        Int_End = windows$start + windows$length, check.names = FALSE)
  if (!is.null(profiles_1$snp_profile)) {
    details <- cbind(details, as.data.frame(profiles_1$snp_profile, optional = TRUE))
  }
  if (!is.null(profiles_1$cov_profile)) {
    profile <- profiles_1$cov_profile
    rownames(profile) <- NULL
    details <- cbind(details, profile)
  }
  if (profile_2_supplied) {
    profile <- profiles_2$cov_profile
    rownames(profile) <- NULL
    names(profile) <- paste0(names(profile), ".2")
    details <- cbind(details, profile)
  }
  details$Risk_Estimate <- risk

  result <- list(risk = matrix(risk, ncol = 1L, dimnames = list(NULL, "Risk_Estimate")),
                 details = details)
  if (return.lp) {
    result$lps.1 <- matrix(part_1$lp, ncol = 1L, dimnames = list(NULL, "Linear_Predictor"))
    result$lps.2 <- matrix(part_2$lp, ncol = 1L, dimnames = list(NULL, "Linear_Predictor"))
  }
  if (return.refs.risk) {
    result$refs.risk.1 <- reference_risks(model_1, windows$start[1], length_1[1])
    result$refs.risk.2 <- reference_risks(model_2, start_2[1], length_2[1])
  }
  class(result) <- "icareSplit"
  result
}

normalize_cut_ages <- function(cut.time, n) {
  cut.time <- as.vector(cut.time)
  if (!is.numeric(cut.time) || anyNA(cut.time) || any(!is.finite(cut.time))) {
    abort("`cut.time` must be numeric with no missing values.")
  }
  if (length(cut.time) == 1L) cut.time <- rep(cut.time, n)
  if (length(cut.time) != n) {
    abort("`cut.time` must have length 1 or one value per profile (", n, ").")
  }
  snap_to_whole(as.numeric(cut.time))
}
