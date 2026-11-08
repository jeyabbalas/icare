# Record the outputs of an earlier version of iCARE for the regression tests
# ==========================================================================
#
# Usage (from the package root):
#
#   Rscript tools/fixtures/reference-outputs.R <library> <kind>
#
# <library> is a library holding the earlier iCARE, and <kind> is
# * "official": Bioconductor iCARE 1.40.0, installed with
#     BiocManager::install("iCARE", lib = "<library>");
# * "fork": iCARE 1.33.1.9002, the patched development version in which the
#     maintainer first fixed the defects listed in vignette("differences").
#
# The outputs go to tests/testthat/fixtures/<kind>-outputs.rds. Each version
# is run in its own R process, because all versions are named iCARE.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[2] %in% c("official", "fork"))
library_path <- args[1]
kind <- args[2]
suppressPackageStartupMessages(library(iCARE, lib.loc = library_path))
message("Recording outputs of iCARE ", utils::packageVersion("iCARE", lib.loc = library_path))

source(file.path("tests", "testthat", "helper-scenarios.R"))
bc <- load_example_data()
inputs <- scenario_inputs(bc)

quietly <- function(expr) {
  result <- NULL
  utils::capture.output(result <- suppressMessages(suppressWarnings(expr)))
  result
}
risk_args <- function(args) {
  c(args, list(use.c.code = 0))
}
keep <- function(result) {
  list(risk = as.numeric(result$risk),
       lps = if (!is.null(result$lps)) as.numeric(result$lps),
       refs_risk = thin(if (!is.null(result$refs.risk)) as.numeric(result$refs.risk)),
       beta = result$beta.used)
}
outputs <- list(version = as.character(utils::packageVersion("iCARE", lib.loc = library_path)))

# ----- Absolute risks
profiles <- bc$ref_cov_dat[inputs$profile_rows, ]
outputs$complete_profiles <- keep(quietly(do.call(computeAbsoluteRisk, risk_args(c(
  covariate_model_args(bc),
  list(apply.age.start = inputs$ages, apply.age.interval.length = inputs$lengths,
       apply.cov.profile = profiles, return.lp = TRUE, return.refs.risk = TRUE))))))

outputs$weighted_reference <- keep(quietly(do.call(computeAbsoluteRisk, risk_args(c(
  covariate_model_args(bc, model.ref.dataset.weights = inputs$weights),
  list(apply.age.start = inputs$ages, apply.age.interval.length = inputs$lengths,
       apply.cov.profile = profiles, return.lp = TRUE))))))

set.seed(1)
outputs$snp_example <- keep(quietly(do.call(computeAbsoluteRisk, risk_args(c(
  snp_model_args(bc),
  list(apply.age.start = 50, apply.age.interval.length = 30,
       apply.cov.profile = bc$new_cov_prof, apply.snp.profile = bc$new_snp_prof,
       return.lp = TRUE, return.refs.risk = TRUE))))))

set.seed(1)
outputs$snp_windows <- keep(quietly(do.call(computeAbsoluteRisk, risk_args(c(
  snp_model_args(bc),
  list(apply.age.start = c(40, 50, 60), apply.age.interval.length = c(5, 10, 20),
       apply.cov.profile = bc$new_cov_prof[c(3, 3, 3), ],
       apply.snp.profile = bc$new_snp_prof[c(3, 3, 3), ],
       return.lp = TRUE, return.refs.risk = TRUE))))))

set.seed(1)
outputs$snp_only <- keep(quietly(computeAbsoluteRisk(
  model.snp.info = bc$bc_72_snps, model.disease.incidence.rates = bc$bc_inc,
  model.competing.incidence.rates = bc$mort_inc, apply.age.start = 50,
  apply.age.interval.length = 30, apply.snp.profile = bc$new_snp_prof,
  return.refs.risk = TRUE, use.c.code = 0)))

# ----- Model validation, with the risks of this version
cohort <- inputs$cohort
ncc <- inputs$ncc
risks_of <- function(data, interval) {
  observed <- data$study.exit.age - data$study.entry.age
  followup <- if (is.null(interval)) observed else pmin(interval, observed)
  result <- quietly(do.call(computeAbsoluteRisk, risk_args(c(
    covariate_model_args(bc),
    list(apply.age.start = data$study.entry.age, apply.age.interval.length = followup,
         apply.cov.profile = data[, inputs$covariates], return.lp = TRUE)))))
  list(risk = as.numeric(result$risk), lp = as.numeric(result$lps))
}
cohort_risk <- risks_of(cohort, 5)
ncc_risk <- risks_of(ncc, NULL)
validate <- function(...) {
  result <- tryCatch(quietly(ModelValidation(...)), error = function(e) conditionMessage(e))
  if (is.character(result)) {
    return(result)
  }
  numeric_column <- function(name) as.numeric(as.character(result$Category_Results[[name]]))
  list(
    auc = result$AUC, variance_auc = result$Variance_AUC,
    eo = result$Overall_Expected_to_Observed_Ratio,
    eo_ci = result$CI_Overall_Expected_to_Observed_Ratio,
    hl = unname(result$Hosmer_Lemeshow_Results$statistic),
    hl_df = unname(result$Hosmer_Lemeshow_Results$parameter),
    hl_p = result$HL_pvalue,
    rr = unname(result$RR_test_result$statistic),
    rr_df = unname(result$RR_test_result$parameter),
    rr_p = result$RR_test_pvalue,
    categories = as.character(result$Category_Results$Categories),
    observed = numeric_column("Observed_Absolute_Risk"),
    expected = numeric_column("Predicted_Absolute_Risk"),
    observed_lower = numeric_column("CI_Absolute_Risk_Lower"),
    observed_upper = numeric_column("CI_Absolute_Risk_Upper"),
    observed_rr = numeric_column("Observed_Relative_Risk"),
    expected_rr = numeric_column("Predicted_Relative_Risk"),
    observed_rr_lower = numeric_column("CI_Relative_Risk_Lower"),
    observed_rr_upper = numeric_column("CI_Relative_Risk_Upper"),
    variance_log_rr = result$Variance_Matrix_LogRelative_Risk,
    predicted_risk = thin(result$Subject_Specific_Predicted_Absolute_Risk),
    cases = sum(result$Subject_Specific_Observed_Outcome),
    followup_sum = sum(result$Adjusted_Followup)
  )
}
outputs$validation <- list(
  cohort = validate(study.data = cohort, predicted.risk.interval = 5,
                    predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp),
  ncc = validate(study.data = ncc, total.followup.validation = TRUE,
                 predicted.risk = ncc_risk$risk, linear.predictor = ncc_risk$lp),
  cutoffs = validate(study.data = cohort, predicted.risk.interval = 5,
                     predicted.risk = cohort_risk$risk, linear.predictor = cohort_risk$lp,
                     linear.predictor.cutoffs = c(-0.4, -0.2, 0, 0.2, 0.4)),
  tied_cohort = validate(study.data = cohort, predicted.risk.interval = 5,
                         predicted.risk = cohort_risk$risk,
                         linear.predictor = round(cohort_risk$lp, 1)),
  tied_ncc = validate(study.data = ncc, total.followup.validation = TRUE,
                      predicted.risk = ncc_risk$risk, linear.predictor = round(ncc_risk$lp, 1)),
  internal = validate(study.data = cohort, predicted.risk.interval = 5,
                      iCARE.model.object = c(covariate_model_args(bc), list(
                        apply.cov.profile = cohort[, inputs$covariates], n.imp = 5,
                        use.c.code = 1, return.lp = TRUE, return.refs.risk = TRUE)))
)

# ----- Risk-score groups (the fork's weighted.quantcut())
if (kind == "fork") {
  outputs$groups <- lapply(group_cases(), function(case) {
    groups <- droplevels(iCARE:::weighted.quantcut(case$x, weights = case$weights,
                                                   q = case$q))
    list(code = as.integer(groups), labels = levels(groups))
  })
}

path <- file.path("tests", "testthat", "fixtures", paste0(kind, "-outputs.rds"))
saveRDS(outputs, path, version = 2, compress = "xz")
message("Saved ", path, " (", format(file.size(path), big.mark = ","), " bytes)")
