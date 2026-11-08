# Example figure of the README and the website's home page
# ========================================================
#
# Draws man/figures/README-validation.png: the validation figure of the
# README's quick start, for the 5-year risks of the example cohort.
#
# Usage, from the root of the repository with this distribution installed:
#
#   Rscript tools/figures/readme-figure.R

suppressPackageStartupMessages(library(iCARE))
data("bc_data", package = "iCARE")

cohort <- validation.cohort.data
validation <- ModelValidation(
  study.data = cohort,
  predicted.risk.interval = 5,
  iCARE.model.object = list(
    model.formula = bc_model_formula,
    model.cov.info = bc_model_cov_info,
    model.log.RR = bc_model_log_or,
    model.ref.dataset = ref_cov_dat,
    model.disease.incidence.rates = bc_inc,
    model.competing.incidence.rates = mort_inc,
    apply.cov.profile = cohort[, all.vars(bc_model_formula)[-1]]
  ),
  dataset = "Simulated cohort",
  model.name = "Breast cancer risk model"
)
figure <- plotModelValidation(validation.results = validation, draw = FALSE)

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
ggplot2::ggsave("man/figures/README-validation.png", figure, device = device,
                width = 8, height = 8.6, dpi = 150, bg = "white")
