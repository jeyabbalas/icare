#' iCARE: Individualized Coherent Absolute Risk Estimation
#'
#' @description
#' iCARE builds, applies and validates models for the absolute risk of a
#' disease. An absolute risk model combines
#'
#' 1. a model for the **relative risks** of the risk factors, typically
#'    estimated in a cohort or case-control study,
#' 2. a **reference dataset** of risk factors representing the population in
#'    which the model is applied,
#' 3. the population's age-specific **disease incidence rates**, and
#' 4. optionally, age-specific rates of **competing mortality**.
#'
#' The main functions are
#'
#' * [computeAbsoluteRisk()]: estimate absolute risks for covariate and/or
#'   SNP profiles, including profiles with missing risk factors;
#' * [computeAbsoluteRiskSplitInterval()]: the same when the relative-risk
#'   model or the risk-factor distribution changes at a cut age;
#' * [ModelValidation()]: validate a model in a cohort or nested
#'   case-control study (calibration and discrimination);
#' * [plotModelValidation()]: publication-ready figures of a validation.
#'
#' @section Unofficial distribution:
#' This is an **unofficial**, refactored distribution of iCARE. The official
#' package is distributed by Bioconductor
#' (<https://bioconductor.org/packages/iCARE/>) and on GitHub
#' (<https://github.com/parichoy/iCARE>). This version keeps the user-facing
#' interface, fixes known defects and has no compiled code; see
#' `vignette("differences", package = "iCARE")` for every difference.
#'
#' @references Pal Choudhury P, Maas P, Wilcox A, et al. (2020). iCARE: An R
#'   package to build, validate and apply absolute risk models. *PLoS ONE*
#'   15(2): e0228198. \doi{10.1371/journal.pone.0228198}
#'
#' @keywords internal
"_PACKAGE"
