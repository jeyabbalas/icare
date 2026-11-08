#' Validate an absolute risk model in a cohort or a sample from a cohort
#'
#' @description
#' Assesses the calibration and discrimination of an absolute risk model in
#' independent data from a cohort study (full cohort) or from a sample drawn
#' from a cohort with known inclusion probabilities, such as a nested
#' case-control study. Risks are validated over a fixed interval after study
#' entry (e.g. 5-year risks) or over each subject's observed follow-up.
#'
#' The function reports
#' * **calibration in the large**: the ratio of the expected to the observed
#'   number of cases (E/O) with its confidence interval;
#' * **calibration by risk group**: subjects are grouped by their risk score
#'   (deciles by default), and the observed and expected absolute risks,
#'   relative risks and E/O ratios of the groups are compared, with a
#'   Hosmer-Lemeshow test of the absolute risks and a test of the relative
#'   risks;
#' * **discrimination**: the area under the ROC curve (AUC) of the risk score;
#' * the study's **incidence rates by attained age**, to compare with the
#'   population incidence rates of the model.
#'
#' @details
#' **Follow-up.** Each subject is followed from `study.entry.age` for the
#' shorter of the risk prediction interval and the observed follow-up, and
#' counts as a case only if the disease develops within that time
#' (`Adjusted_Followup`, `Subject_Specific_Observed_Outcome`). When the risks
#' are computed here, each subject's risk covers that same time from
#' `study.entry.age`, so that a subject followed for 3 years of a 5-year
#' interval contributes a 3-year risk. To validate fixed-horizon risks,
#' compute them yourself and supply `predicted.risk` and `linear.predictor`.
#' Cases whose onset is recorded after exit are treated as non-cases censored
#' at exit, with a warning.
#'
#' **Sampled data.** When `study.data` has a `sampling.weights` column (each
#' subject's probability of inclusion in the sample), all estimates are
#' weighted by the inverse probabilities, and the variances add the
#' variability due to sampling.
#'
#' **Risk-score groups.** Subjects are grouped at quantiles of
#' `linear.predictor` (weighted quantiles in sampled data), or at
#' `linear.predictor.cutoffs`. Values tied at a quantile form their own group,
#' and groups left empty are dropped; the tests use the number of groups
#' actually formed, \eqn{G}.
#'
#' **Tests.** The Hosmer-Lemeshow statistic is
#' \eqn{\sum_g (O_g - E_g)^2 / V_g} on \eqn{G} degrees of freedom, where
#' \eqn{O_g} and \eqn{E_g} are the observed and expected absolute risks of
#' group \eqn{g}. With `hosmer.lemeshow.variance = "expected"` (the default),
#' \eqn{V_g} is the variance under the null hypothesis of calibration,
#' \eqn{E_g(1 - E_g)/N_g} (plus a sampling term in sampled data), as in the
#' classical Hosmer-Lemeshow test. With `"observed"`, \eqn{V_g} is the
#' variance of \eqn{O_g}, as in earlier versions of iCARE; that test rejects
#' too often when groups have few cases and is infinite when a group has
#' none. The relative-risk test compares the observed and expected relative
#' risks of the groups (each group's risk divided by the mean over groups) on
#' \eqn{G - 1} degrees of freedom; it needs a case in every group.
#'
#' **Confidence intervals** are 95% Wald intervals: on the log scale for E/O
#' ratios and relative risks, and on the risk scale for absolute risks and the
#' AUC. The E/O interval in a full cohort treats the expected risks as known;
#' in sampled data it also accounts for the sampling variability of the
#' expected risk and its covariance with the observed risk.
#'
#' **AUC.** The probability that a case has a higher risk score than a
#' non-case, counting ties as 1/2, weighted by the inverse inclusion
#' probabilities in sampled data, with the variance of DeLong et al. (1988).
#' The risk score is the linear predictor, which excludes age, so the AUC
#' measures discrimination by the risk factors within the age distribution of
#' the study.
#'
#' **Study incidence.** Onsets divided by person-years at each year of
#' attained age over the observed follow-up, with exact Poisson intervals in
#' a full cohort and weighted estimates with Wald intervals in sampled data.
#'
#' @param study.data Data frame with one row per subject and the columns
#'   * `observed.outcome`: 1 if the disease developed during follow-up, 0
#'     otherwise;
#'   * `study.entry.age`: age at entry into the study;
#'   * `study.exit.age`: age at the last follow-up;
#'   * `time.of.onset`: years from entry to disease onset (`Inf` or `NA` for
#'     subjects without the disease);
#'   * `sampling.weights` (optional): for a sample from the cohort, such as a
#'     nested case-control study, each subject's probability of inclusion in
#'     the sample.
#' @param total.followup.validation If `TRUE`, validate risks over each
#'   subject's observed follow-up; if `FALSE` (default), over
#'   `predicted.risk.interval`.
#' @param predicted.risk Optional vector of predicted absolute risks, one per
#'   subject. If `NULL`, the risks are computed from `iCARE.model.object`.
#' @param predicted.risk.interval Years after entry over which risks are
#'   validated (e.g. 5): one value, or one per subject. Not used when
#'   `total.followup.validation = TRUE`.
#' @param linear.predictor Optional vector of risk scores (the linear
#'   predictor \eqn{x'\beta}), one per subject. If `NULL`, computed from
#'   `iCARE.model.object`. A supplied `predicted.risk` or `linear.predictor`
#'   is always used as given; only a missing one is computed.
#' @param iCARE.model.object Named list of arguments of
#'   [computeAbsoluteRisk()] describing the model and the subjects:
#'   `model.formula`, `model.cov.info`, `model.snp.info`, `model.log.RR`,
#'   `model.ref.dataset`, `model.ref.dataset.weights`,
#'   `model.disease.incidence.rates`, `model.competing.incidence.rates`,
#'   `model.bin.fh.name`, `apply.cov.profile` and `apply.snp.profile` (the
#'   subjects' risk factors, one row per row of `study.data`) and `n.imp`.
#'   Settings that earlier versions needed here (`use.c.code`, `return.lp`,
#'   `return.refs.risk`) are accepted and ignored.
#' @param number.of.percentiles Number of risk-score groups formed at
#'   quantiles of the risk score (default 10, deciles).
#' @param reference.entry.age,reference.exit.age Optional ages (one value, or
#'   one per row of `model.ref.dataset`) over which to compute the absolute
#'   risks of the reference population, returned as `Reference_Absolute_Risk`.
#' @param predicted.risk.ref,linear.predictor.ref Optional absolute risks and
#'   risk scores of the reference population, returned as given.
#' @param linear.predictor.cutoffs Optional cut points of the risk score
#'   defining the risk groups, instead of quantiles.
#' @param dataset,model.name Names of the dataset and the model, shown by
#'   `print()` and [plotModelValidation()].
#' @param hosmer.lemeshow.variance Variance in the Hosmer-Lemeshow statistic:
#'   `"expected"` (default) or `"observed"` (earlier versions); see Details.
#'
#' @return An object of class `"icareValid"`: a list with
#' * `Subject_Specific_Observed_Outcome`, `Adjusted_Followup`: outcome and
#'   follow-up over the risk prediction interval;
#' * `Risk_Prediction_Interval`: the interval, e.g. `"5 years"` or
#'   `"Observed Followup"`;
#' * `Subject_Specific_Predicted_Absolute_Risk`, `Subject_Specific_Risk_Score`:
#'   each subject's predicted risk and risk score;
#' * `Reference_Absolute_Risk`, `Reference_Risk_Score`: those of the
#'   reference population (`NULL` unless supplied or requested);
#' * `Population_Incidence_Rate`: the model's disease incidence rates;
#' * `Study_Incidence_Rate`: data frame of the study's incidence rates by
#'   attained age (`ages`, `study_incidence`), their confidence limits
#'   (`lower`, `upper`), and the (weighted) `cases` and `person_years`;
#' * `Category_Results`: data frame with one row per risk-score group: the
#'   group (`Categories`); its observed and predicted absolute risks and
#'   relative risks with confidence intervals for the observed ones; its E/O
#'   ratio with confidence interval (`Expected_to_Observed_Ratio`,
#'   `CI_Expected_to_Observed_Ratio_Lower`, `CI_Expected_to_Observed_Ratio_Upper`);
#'   and its numbers of subjects and cases in the data;
#' * `Category_Specific_Observed_Absolute_Risk`,
#'   `Category_Specific_Predicted_Absolute_Risk`,
#'   `Category_Specific_Observed_Relative_Risk`,
#'   `Category_Specific_Predicted_Relative_Risk`,
#'   `Category_Specific_Expected_to_Observed_Ratio`: the same, as vectors;
#' * `Variance_Matrix_Absolute_Risk`, `Variance_Matrix_LogRelative_Risk`:
#'   covariance matrices of the groups' observed absolute risks and log
#'   relative risks;
#' * `Hosmer_Lemeshow_Results`, `HL_pvalue`: the Hosmer-Lemeshow test (an
#'   `htest` object) and its p-value;
#' * `RR_test_result`, `RR_test_pvalue`: the relative-risk test and its
#'   p-value;
#' * `AUC`, `Variance_AUC`, `CI_AUC`: the AUC, its variance and 95%
#'   confidence interval;
#' * `Overall_Expected_to_Observed_Ratio`,
#'   `CI_Overall_Expected_to_Observed_Ratio`: E/O and its 95% confidence
#'   interval;
#' * `input.args`: the model formula, the study data, the names and the
#'   settings, used by `print()` and `plot()`.
#'
#' @references
#' Pal Choudhury P, Maas P, Wilcox A, et al. (2020). iCARE: An R package to
#' build, validate and apply absolute risk models. *PLOS ONE* 15(2):
#' e0228198. \doi{10.1371/journal.pone.0228198}
#'
#' DeLong ER, DeLong DM, Clarke-Pearson DL (1988). Comparing the areas under
#' two or more correlated receiver operating characteristic curves: a
#' nonparametric approach. *Biometrics* 44(3): 837-845.
#'
#' Hosmer DW, Lemeshow S (1980). Goodness of fit tests for the multiple
#' logistic regression model. *Communications in Statistics - Theory and
#' Methods* 9(10): 1043-1069.
#'
#' @seealso [plotModelValidation()], [computeAbsoluteRisk()]
#'
#' @examples
#' data(bc_data, package = "iCARE")
#'
#' # 5-year risks in a random sample of the full validation cohort
#' set.seed(50)
#' cohort <- validation.cohort.data[sample(nrow(validation.cohort.data), 5000), ]
#' covariates <- all.vars(bc_model_formula)[-1]
#' validation <- ModelValidation(
#'   study.data = cohort,
#'   predicted.risk.interval = 5,
#'   iCARE.model.object = list(
#'     model.formula = bc_model_formula,
#'     model.cov.info = bc_model_cov_info,
#'     model.log.RR = bc_model_log_or,
#'     model.ref.dataset = ref_cov_dat,
#'     model.disease.incidence.rates = bc_inc,
#'     model.competing.incidence.rates = mort_inc,
#'     apply.cov.profile = cohort[, covariates]
#'   ),
#'   dataset = "Simulated cohort",
#'   model.name = "Breast cancer risk model"
#' )
#' validation
#' validation$Category_Results
#' @export
ModelValidation <- function(study.data,
                            total.followup.validation = FALSE,
                            predicted.risk = NULL,
                            predicted.risk.interval = NULL,
                            linear.predictor = NULL,
                            iCARE.model.object = list(
                              model.formula = NULL, model.cov.info = NULL,
                              model.snp.info = NULL, model.log.RR = NULL,
                              model.ref.dataset = NULL,
                              model.ref.dataset.weights = NULL,
                              model.disease.incidence.rates = NULL,
                              model.competing.incidence.rates = NULL,
                              model.bin.fh.name = NA,
                              apply.cov.profile = NULL,
                              apply.snp.profile = NULL, n.imp = 5),
                            number.of.percentiles = 10,
                            reference.entry.age = NULL,
                            reference.exit.age = NULL,
                            predicted.risk.ref = NULL,
                            linear.predictor.ref = NULL,
                            linear.predictor.cutoffs = NULL,
                            dataset = "Example Dataset",
                            model.name = "Example Risk Prediction Model",
                            hosmer.lemeshow.variance = c("expected", "observed")) {
  hosmer.lemeshow.variance <- match.arg(hosmer.lemeshow.variance)
  total.followup.validation <- as_flag(total.followup.validation,
                                       "total.followup.validation")
  model_inputs <- read_model_object(iCARE.model.object)

  # ----- Study data and follow-up over the risk prediction interval
  study <- validate_study_data(study.data)
  n <- length(study$outcome)
  interval <- normalize_prediction_interval(predicted.risk.interval,
                                            total.followup.validation, n)
  followup <- adjust_validation_followup(study, interval)
  outcome <- followup$outcome
  if (sum(outcome) == 0 || sum(outcome) == n) {
    abort("The study needs both cases and non-cases over the risk prediction ",
          "interval; it has ", sum(outcome), " cases among ", n, " subjects.")
  }

  # ----- Predicted risks and risk scores (only the missing ones are computed)
  model <- NULL
  get_model <- function(pieces, need_profile) {
    if (is.null(model)) {
      check_validation_model_object(model_inputs, pieces, need_profile)
      model <<- do.call(build_risk_model,
                         model_inputs[intersect(model_argument_names, names(model_inputs))])
    }
    model
  }
  missing_pieces <- c("predicted.risk", "linear.predictor")[
    c(is.null(predicted.risk), is.null(linear.predictor))]
  if (length(missing_pieces) > 0L) {
    model <- get_model(missing_pieces, need_profile = TRUE)
    profiles <- prepare_profiles(model, model_inputs$apply.cov.profile,
                                 model_inputs$apply.snp.profile, n_default = n)
    if (profiles$n != n) {
      abort("The profiles in `iCARE.model.object` (`apply.cov.profile`, ",
            "`apply.snp.profile`) must have one row per row of `study.data` (", n,
            "); they have ", profiles$n, ".")
    }
    prediction <- predict_risks(model, profiles, study$entry, followup$followup)
    if (is.null(predicted.risk)) predicted.risk <- prediction$risk
    if (is.null(linear.predictor)) linear.predictor <- prediction$lp
  }
  predicted.risk <- check_subject_values(predicted.risk, "predicted.risk", n,
                                         lower = 0, upper = 1)
  linear.predictor <- check_subject_values(linear.predictor, "linear.predictor", n)

  # ----- Reference population (optional)
  reference_missing <- c("predicted.risk.ref", "linear.predictor.ref")[
    c(is.null(predicted.risk.ref), is.null(linear.predictor.ref))]
  if (length(reference_missing) > 0L && !is.null(reference.entry.age) &&
      !is.null(reference.exit.age)) {
    model <- get_model(reference_missing, need_profile = FALSE)
    windows <- reference_windows(model, reference.entry.age, reference.exit.age)
    if (is.null(predicted.risk.ref)) {
      predicted.risk.ref <- reference_risks(model, windows$start, windows$length)
    }
    if (is.null(linear.predictor.ref)) {
      linear.predictor.ref <- reference_linear_predictors(model)
    }
  }

  # ----- Discrimination, calibration and incidence
  weights <- if (is.null(study$sampling_prob)) NULL else 1 / study$sampling_prob
  auc <- compute_auc(linear.predictor, outcome, study$sampling_prob)
  overall <- overall_expected_observed(predicted.risk, outcome, weights)
  number.of.percentiles <- as_count(number.of.percentiles, "number.of.percentiles",
                                    lower = 2)
  groups <- risk_score_groups(linear.predictor, weights, number.of.percentiles,
                              linear.predictor.cutoffs)
  calibration <- risk_group_calibration(predicted.risk, outcome, groups, weights,
                                        hosmer.lemeshow.variance)
  incidence <- study_incidence_by_age(study$entry, study$exit, followup$onset, weights)

  # ----- Results
  n_groups <- nlevels(groups)
  group_description <- paste0("observed and expected risks in ", n_groups,
                              " risk-score groups")
  hl_test <- structure(
    list(statistic = c(Chisquare = calibration$hl$statistic),
         parameter = c(df = calibration$hl$df),
         p.value = calibration$hl$p_value,
         method = "Hosmer and Lemeshow goodness of fit (GOF) test for Absolute Risk",
         data.name = paste0(group_description, "; variance from the ",
                            hosmer.lemeshow.variance, " risks")),
    class = "htest")
  rr_test <- structure(
    list(statistic = c(Chisquare = calibration$rr$statistic),
         parameter = c(df = calibration$rr$df),
         p.value = calibration$rr$p_value,
         method = "Goodness of fit (GOF) test for Relative Risk",
         data.name = group_description),
    class = "htest")

  category_results <- data.frame(
    Categories = levels(groups),
    Observed_Absolute_Risk = calibration$observed,
    Predicted_Absolute_Risk = calibration$expected,
    CI_Absolute_Risk_Lower = calibration$observed_lower,
    CI_Absolute_Risk_Upper = calibration$observed_upper,
    Observed_Relative_Risk = calibration$observed_rr,
    Predicted_Relative_Risk = calibration$expected_rr,
    CI_Relative_Risk_Lower = calibration$observed_rr_lower,
    CI_Relative_Risk_Upper = calibration$observed_rr_upper,
    Expected_to_Observed_Ratio = calibration$eo_ratio,
    CI_Expected_to_Observed_Ratio_Lower = calibration$eo_lower,
    CI_Expected_to_Observed_Ratio_Upper = calibration$eo_upper,
    Number_of_Subjects = calibration$subjects,
    Number_of_Cases = calibration$cases,
    stringsAsFactors = FALSE
  )

  result <- list(
    Subject_Specific_Observed_Outcome = outcome,
    Risk_Prediction_Interval = describe_prediction_interval(interval),
    Adjusted_Followup = followup$followup,
    Subject_Specific_Predicted_Absolute_Risk = predicted.risk,
    Reference_Absolute_Risk = predicted.risk.ref,
    Subject_Specific_Risk_Score = linear.predictor,
    Reference_Risk_Score = linear.predictor.ref,
    Population_Incidence_Rate = model_inputs$model.disease.incidence.rates,
    Study_Incidence_Rate = incidence,
    Category_Results = category_results,
    Category_Specific_Observed_Absolute_Risk = calibration$observed,
    Category_Specific_Predicted_Absolute_Risk = calibration$expected,
    Category_Specific_Observed_Relative_Risk = calibration$observed_rr,
    Category_Specific_Predicted_Relative_Risk = calibration$expected_rr,
    Category_Specific_Expected_to_Observed_Ratio = calibration$eo_ratio,
    Variance_Matrix_Absolute_Risk = calibration$variance_matrix,
    Variance_Matrix_LogRelative_Risk = calibration$variance_log_rr,
    Hosmer_Lemeshow_Results = hl_test,
    HL_pvalue = calibration$hl$p_value,
    RR_test_result = rr_test,
    RR_test_pvalue = calibration$rr$p_value,
    AUC = auc$auc,
    Variance_AUC = auc$variance,
    CI_AUC = auc$ci,
    Overall_Expected_to_Observed_Ratio = overall$ratio,
    CI_Overall_Expected_to_Observed_Ratio = overall$ci,
    input.args = list(model.formula = model_inputs$model.formula,
                      study.data = study.data, dataset = dataset,
                      model.name = model.name,
                      hosmer.lemeshow.variance = hosmer.lemeshow.variance)
  )
  class(result) <- "icareValid"
  result
}

# Arguments of build_risk_model() that iCARE.model.object can hold
model_argument_names <- c("model.formula", "model.cov.info", "model.snp.info",
                          "model.log.RR", "model.ref.dataset",
                          "model.ref.dataset.weights", "model.disease.incidence.rates",
                          "model.competing.incidence.rates", "model.bin.fh.name",
                          "n.imp")

#' Read `iCARE.model.object`, with defaults for the settings left out
#' @noRd
read_model_object <- function(object) {
  if (is.null(object)) object <- list()
  if (!is.list(object) || is.data.frame(object)) {
    abort("`iCARE.model.object` must be a named list of arguments of ",
          "computeAbsoluteRisk().")
  }
  known <- c(model_argument_names, "apply.cov.profile", "apply.snp.profile")
  ignored <- c("use.c.code", "return.lp", "return.refs.risk")
  supplied <- names(object)
  if (length(object) > 0L && (is.null(supplied) || any(!nzchar(supplied)))) {
    abort("Every element of `iCARE.model.object` must be named.")
  }
  unknown <- setdiff(supplied, c(known, ignored))
  if (length(unknown) > 0L) {
    warn("Ignoring unknown element(s) of `iCARE.model.object`: ",
         format_names(unknown), ". Names must match the arguments of ",
         "computeAbsoluteRisk() exactly.")
  }
  inputs <- object[intersect(supplied, known)]
  if (is.null(inputs$n.imp)) inputs$n.imp <- 5
  if (is.null(inputs$model.bin.fh.name)) inputs$model.bin.fh.name <- NA
  inputs
}

#' Stop with a clear message when the model object cannot give what is needed
#' @param pieces Names of the quantities to compute.
#' @param need_profile Whether the subjects' profiles are needed.
#' @noRd
check_validation_model_object <- function(inputs, pieces, need_profile) {
  has <- function(name) !is.null(inputs[[name]])
  covariate_args <- c("model.formula", "model.cov.info", "model.log.RR",
                      "model.ref.dataset", if (need_profile) "apply.cov.profile")
  complete <- has("model.disease.incidence.rates") &&
    (all(vapply(covariate_args, has, logical(1))) || has("model.snp.info"))
  if (!complete) {
    pieces_text <- paste0("`", pieces, "`", collapse = " and ")
    abort(pieces_text, if (length(pieces) > 1L) " were" else " was", " not supplied, ",
          "so ", if (length(pieces) > 1L) "they" else "it", " must be computed from ",
          "`iCARE.model.object`, which needs `model.disease.incidence.rates` and ",
          "either all of ", format_names(paste0("`", covariate_args, "`")), " or ",
          "`model.snp.info`. Supply ", pieces_text, " directly (to group subjects by ",
          "absolute risk, use `linear.predictor = log(predicted.risk)`), or complete ",
          "`iCARE.model.object`.")
  }
  invisible(TRUE)
}

#' Check one value per subject: numeric, finite, within bounds
#' @noRd
check_subject_values <- function(x, name, n, lower = -Inf, upper = Inf) {
  x <- as.vector(x)
  if (!is.numeric(x) || length(x) != n) {
    abort("`", name, "` must be a numeric vector with one value per row of ",
          "`study.data` (", n, "); it has ", length(x), ".")
  }
  if (anyNA(x) || any(!is.finite(x))) {
    abort("`", name, "` cannot contain missing or infinite values.")
  }
  if (any(x < lower | x > upper)) {
    abort("`", name, "` must lie between ", lower, " and ", upper, ".")
  }
  as.numeric(x)
}

#' Windows of the reference population from entry and exit ages
#' @noRd
reference_windows <- function(model, entry, exit) {
  n_reference <- sum(model$population$block == 1L)
  check <- function(x, name) {
    x <- as.vector(x)
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
        !length(x) %in% c(1L, n_reference)) {
      abort("`", name, "` must be numeric with one value, or one per row of the ",
            "reference population (", n_reference, ").")
    }
    as.numeric(x)
  }
  entry <- check(entry, "reference.entry.age")
  exit <- check(exit, "reference.exit.age")
  if (any(exit <= entry)) {
    abort("`reference.exit.age` must be greater than `reference.entry.age`.")
  }
  length <- exit - entry
  if (base::length(entry) != base::length(length)) {
    entry <- rep_len(entry, base::length(length))
  }
  normalize_risk_windows(entry, length, base::length(length), "model.ref.dataset")
}
