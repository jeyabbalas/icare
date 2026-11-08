# Validation results shared by the plot and print tests
example_validations <- function() {
  bc <- load_example_data()
  inputs <- scenario_inputs(bc)
  cohort <- inputs$cohort
  ncc <- inputs$ncc
  covariates <- inputs$covariates
  object <- covariate_model_args(bc)
  list(
    cohort = ModelValidation(
      study.data = cohort, predicted.risk.interval = 5,
      iCARE.model.object = c(object, list(apply.cov.profile = cohort[, covariates])),
      dataset = "Example cohort", model.name = "Breast cancer model"),
    ncc = ModelValidation(
      study.data = ncc, total.followup.validation = TRUE,
      iCARE.model.object = c(object, list(apply.cov.profile = ncc[, covariates])),
      dataset = "Example nested case-control study", model.name = "Breast cancer model")
  )
}
