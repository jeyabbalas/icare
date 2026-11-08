#' Build and apply an absolute risk model
#'
#' @description
#' Estimates the absolute risk of a disease over an age interval for each
#' covariate (and/or SNP) profile. The model combines
#'
#' * a model for the **relative risks** of the risk factors (log relative
#'   risks of the design-matrix columns of `model.formula`, and/or SNP odds
#'   ratios),
#' * a **reference dataset** of risk factors representing the population in
#'   which the model is applied,
#' * the population's age-specific **disease incidence rates**, and
#' * optionally, age-specific **competing mortality rates**.
#'
#' The baseline hazard is calibrated so that the model reproduces the
#' population incidence rates, and risks of profiles with missing risk
#' factors are averaged coherently over the reference population.
#'
#' @details
#' **Model.** The hazard of a person with risk factors \eqn{z} at age
#' \eqn{t} is \eqn{\lambda_0(t)\exp(\beta' x(z))}. The absolute risk over
#' \eqn{[a, a + \tau)} is computed over whole years of age \eqn{t}, as in
#' earlier versions of iCARE:
#' \deqn{AR = \sum_t e_t \lambda_0(t) RR \exp\{-\sum_{u<t} e_u
#'   (\lambda_0(u) RR + m(u))\},}
#' where \eqn{RR = \exp(\beta' x(z))}, \eqn{m} is the competing mortality
#' rate and \eqn{e_t} is the time the interval spends in the year of age
#' \eqn{[t, t + 1)}: 1 inside the interval for whole-number ages, and the
#' corresponding fraction of a year at fractional ends.
#'
#' **Calibration.** The baseline hazard solves
#' \eqn{\lambda(t) = \lambda_0(t) E[RR \mid T \ge t]}, the expectation being
#' over the (weighted) reference population still free of the disease at age
#' \eqn{t}, starting from the youngest age of the incidence table.
#'
#' **Missing risk factors.** For a profile with missing risk factors, the risk
#' is \eqn{E[AR(z_o, Z_m) \mid z_o]}: reference rows are grouped into 100
#' percentile bins of the linear predictor of the observed risk factors, and
#' the risk is the weighted mean, over the rows in the profile's bin, of the
#' risk of the profile's own observed risk factors completed with each row's
#' values of the missing ones. A profile with no observed risk factor gets the
#' population-average risk.
#'
#' **SNPs.** SNP genotypes of the reference population are simulated under
#' Hardy-Weinberg equilibrium from the allele frequencies (`n.imp` sets),
#' conditionally on family history when `model.bin.fh.name` is given, in which
#' case the family-history log relative risk is attenuated by half the
#' variance of the polygenic score. Missing genotypes in `apply.snp.profile`
#' are handled like other missing risk factors. Results depend on the random
#' number generator; call [set.seed()] first for reproducible results.
#'
#' @param model.formula A formula such as `Y ~ x1 + as.factor(x2) + x1:x3`
#'   describing the design matrix of the relative-risk model. The response is
#'   ignored. Leave `NULL` (with the other covariate arguments) for a SNP-only
#'   model.
#' @param model.cov.info A list with one element per covariate of
#'   `model.formula`, each a list with fields `name` (the covariate name),
#'   `type` (`"continuous"` or `"factor"`) and, for factors, `levels` (all
#'   levels) and optionally `ref` (the referent level).
#' @param model.snp.info A data frame with columns `snp.name`,
#'   `snp.odds.ratio` (per-allele odds ratio) and `snp.freq` (risk-allele
#'   frequency).
#' @param model.log.RR Log relative risks of the design-matrix columns (no
#'   intercept), as a named vector or a one-column matrix with row names; the
#'   names must be the column names that `model.formula` generates.
#' @param model.ref.dataset A data frame of risk factors for a sample
#'   representing the population, with no missing values. Columns are matched
#'   to the covariates by name.
#' @param model.ref.dataset.weights Optional non-negative sampling weights of
#'   the rows of `model.ref.dataset` (normalized internally).
#' @param model.disease.incidence.rates Age-specific annual incidence rates:
#'   a two-column table (integer age, rate) or a three-column table (start
#'   age, end age, annual rate in that age band). Must cover consecutive ages
#'   including every year of age used by the risk intervals.
#' @param model.competing.incidence.rates Age-specific annual rates of
#'   competing mortality, in the same layouts. `NULL` means no competing
#'   mortality.
#' @param model.bin.fh.name Name of a binary (0/1) family-history covariate,
#'   used with `model.snp.info` to simulate genotypes given family history
#'   and to attenuate its log relative risk. `NA` when not applicable.
#' @param n.imp Number of simulated SNP genotype sets of the reference
#'   population (1 to 20).
#' @param apply.age.start Age(s) at the start of the risk interval: a single
#'   value or one per profile. Fractional ages are allowed.
#' @param apply.age.interval.length Length(s) of the risk interval in years:
#'   a single value or one per profile. Fractional lengths are allowed.
#' @param apply.cov.profile A data frame of covariate profiles (columns
#'   matched by name; extra columns are ignored); missing values are allowed.
#' @param apply.snp.profile A matrix or data frame of genotypes (risk-allele
#'   counts 0, 1, 2 or `NA`), one column per SNP of `model.snp.info` (matched
#'   by name when the names are present, otherwise by position). If `NULL`,
#'   all SNPs are treated as missing.
#' @param return.lp If `TRUE`, also return the linear predictors.
#' @param return.refs.risk If `TRUE`, also return the risks of the reference
#'   population over the first profile's interval.
#' @param ... Not used. The removed argument `use.c.code` is accepted with a
#'   deprecation warning and ignored.
#'
#' @return An object of class `"icare"`, a list with
#' * `risk`: one-column matrix (`Risk_Estimate`) of absolute risks;
#' * `details`: data frame with the interval (`Int_Start`, `Int_End`), the
#'   risk (`Risk_Estimate`) and the profile (genotypes, then covariates);
#' * `beta.used`: one-column matrix of the log relative risks used (SNPs
#'   first; the family-history coefficient after attenuation);
#' * `lps`: one-column matrix of linear predictors (if `return.lp`); for
#'   profiles with missing risk factors, the model-averaged linear predictor;
#' * `refs.risk`: risks of the reference rows (if `return.refs.risk`).
#'
#' The attribute `"calibration"` holds the calibrated baseline hazard: a data
#' frame with the `age`, the population `incidence.rate`, the
#' `baseline.hazard` and the survivor-weighted mean relative risk
#' `expected.rr` at each age of the incidence table.
#'
#' @references Pal Choudhury P, Maas P, Wilcox A, et al. (2020). iCARE: An R
#'   package to build, validate and apply absolute risk models. *PLoS ONE*
#'   15(2): e0228198. \doi{10.1371/journal.pone.0228198}
#'
#' @seealso [computeAbsoluteRiskSplitInterval()] for models that change at a
#'   cut age, [ModelValidation()] to validate a model.
#'
#' @examples
#' data(bc_data, package = "iCARE")
#'
#' # Covariates and 72 SNPs, with family-history adjustment
#' set.seed(50)
#' results <- computeAbsoluteRisk(
#'   model.formula = bc_model_formula,
#'   model.cov.info = bc_model_cov_info,
#'   model.snp.info = bc_72_snps,
#'   model.log.RR = bc_model_log_or,
#'   model.ref.dataset = ref_cov_dat,
#'   model.disease.incidence.rates = bc_inc,
#'   model.competing.incidence.rates = mort_inc,
#'   model.bin.fh.name = "famhist",
#'   apply.age.start = 50,
#'   apply.age.interval.length = 30,
#'   apply.cov.profile = new_cov_prof,
#'   apply.snp.profile = new_snp_prof,
#'   return.refs.risk = TRUE
#' )
#' results
#' summary(results)
#'
#' # A SNP-only model
#' snp_only <- computeAbsoluteRisk(
#'   model.snp.info = bc_72_snps,
#'   model.disease.incidence.rates = bc_inc,
#'   model.competing.incidence.rates = mort_inc,
#'   apply.age.start = 50,
#'   apply.age.interval.length = 30,
#'   apply.snp.profile = new_snp_prof
#' )
#' snp_only$risk
#' @export
computeAbsoluteRisk <- function(model.formula = NULL, model.cov.info = NULL,
                                model.snp.info = NULL, model.log.RR = NULL,
                                model.ref.dataset = NULL,
                                model.ref.dataset.weights = NULL,
                                model.disease.incidence.rates,
                                model.competing.incidence.rates = NULL,
                                model.bin.fh.name = NA, n.imp = 5,
                                apply.age.start, apply.age.interval.length,
                                apply.cov.profile = NULL,
                                apply.snp.profile = NULL,
                                return.lp = FALSE, return.refs.risk = FALSE, ...) {
  check_removed_arguments(list(...), "computeAbsoluteRisk")
  return.lp <- as_flag(return.lp, "return.lp")
  return.refs.risk <- as_flag(return.refs.risk, "return.refs.risk")
  if (missing(model.disease.incidence.rates)) {
    abort("`model.disease.incidence.rates` must be supplied.")
  }
  if (missing(apply.age.start) || missing(apply.age.interval.length)) {
    abort("`apply.age.start` and `apply.age.interval.length` must be supplied.")
  }

  model <- build_risk_model(
    model.formula = model.formula, model.cov.info = model.cov.info,
    model.snp.info = model.snp.info, model.log.RR = model.log.RR,
    model.ref.dataset = model.ref.dataset,
    model.ref.dataset.weights = model.ref.dataset.weights,
    model.disease.incidence.rates = model.disease.incidence.rates,
    model.competing.incidence.rates = model.competing.incidence.rates,
    model.bin.fh.name = model.bin.fh.name, n.imp = n.imp
  )
  profiles <- prepare_profiles(
    model, apply.cov.profile, apply.snp.profile,
    n_default = default_profile_count(apply.age.start, apply.age.interval.length)
  )
  windows <- normalize_risk_windows(
    apply.age.start, apply.age.interval.length, profiles$n,
    if (model$has_covariates) "apply.cov.profile" else "apply.snp.profile"
  )
  prediction <- predict_risks(model, profiles, windows$start, windows$length)

  result <- package_absolute_risks(model, profiles, windows, prediction, return.lp)
  if (return.refs.risk) {
    result$refs.risk <- reference_risks(model, windows$start[1], windows$length[1])
  }
  attr(result, "calibration") <- calibration_table(model)
  class(result) <- "icare"
  result
}

# Number of profiles of a SNP-only model without genotypes (earlier versions:
# 10,000 for a single interval, otherwise one per interval).
default_profile_count <- function(start, length) {
  n <- max(base::length(start), base::length(length))
  if (n == 1L) snp_only_population_size else n
}

#' Assemble the results list of computeAbsoluteRisk()
#' @noRd
package_absolute_risks <- function(model, profiles, windows, prediction, return.lp) {
  n <- profiles$n
  details <- data.frame(Int_Start = windows$start,
                        Int_End = windows$start + windows$length,
                        Risk_Estimate = prediction$risk,
                        check.names = FALSE)
  if (!is.null(profiles$snp_profile)) {
    details <- cbind(details, as.data.frame(profiles$snp_profile, optional = TRUE))
  }
  if (!is.null(profiles$cov_profile)) {
    cov_profile <- profiles$cov_profile
    rownames(cov_profile) <- NULL
    details <- cbind(details, cov_profile)
  }
  result <- list(
    risk = matrix(prediction$risk, ncol = 1L, dimnames = list(NULL, "Risk_Estimate")),
    details = details,
    beta.used = matrix(model$beta, ncol = 1L,
                       dimnames = list(names(model$beta), "log OR used"))
  )
  if (return.lp) {
    result$lps <- matrix(prediction$lp, ncol = 1L, dimnames = list(NULL, "Linear_Predictor"))
  }
  result
}

calibration_table <- function(model) {
  data.frame(age = model$baseline$age,
             incidence.rate = model$incidence$rate,
             baseline.hazard = model$baseline$hazard,
             expected.rr = model$baseline$expected_rr)
}
