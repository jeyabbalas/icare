# Time iCARE on standard workloads
# ================================
#
# Usage (from the package root):
#
#   Rscript tools/benchmarks/run-benchmarks.R <library|dev> <label> [use.c.code]
#
# <library> is a library holding a version of iCARE ("dev" loads this source
# tree with pkgload; to time this version byte-compiled, as users run it,
# install it into a library with R CMD INSTALL --library=<library> . and pass
# that library), <label> names the results file
# tools/benchmarks/results-<label>.csv, and use.c.code (0 or 1) is passed to
# versions that have compiled code. Each version must run in its own R
# process, because all versions are named iCARE. Timings are the median of
# three runs (one run when the first takes longer than 20 seconds); memory is
# what R allocated during one run (allocations inside compiled code are not
# counted).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
library_path <- args[1]
label <- args[2]
use_c_code <- if (length(args) >= 3L) as.numeric(args[3]) else NULL

if (library_path == "dev") {
  pkgload::load_all(".", quiet = TRUE)
  version <- paste(read.dcf("DESCRIPTION", fields = "Version")[1, 1], "(dev)")
} else {
  suppressPackageStartupMessages(library(iCARE, lib.loc = library_path))
  version <- as.character(utils::packageVersion("iCARE", lib.loc = library_path))
}
message("Benchmarking iCARE ", version, " as '", label, "'",
        if (!is.null(use_c_code)) paste0(" with use.c.code = ", use_c_code))

env <- new.env()
utils::data("bc_data", package = "iCARE", envir = env)
bc <- as.list(env)

# ----- Inputs (the same in every version) -----------------------------------
set.seed(2025)
n_complete <- 10000L
complete_profiles <- bc$ref_cov_dat[sample(nrow(bc$ref_cov_dat), n_complete, TRUE), ]
complete_ages <- sample(30:70, n_complete, TRUE)
complete_lengths <- sample(1:15, n_complete, TRUE)

n_missing <- 2000L
missing_profiles <- bc$ref_cov_dat[sample(nrow(bc$ref_cov_dat), n_missing, TRUE), ]
for (v in names(missing_profiles)) {
  missing_profiles[[v]][stats::runif(n_missing) < 0.15] <- NA
}
missing_ages <- sample(30:70, n_missing, TRUE)
missing_lengths <- sample(1:15, n_missing, TRUE)

n_snp <- 1000L
snp_cov_profiles <- bc$ref_cov_dat[sample(nrow(bc$ref_cov_dat), n_snp, TRUE), ]
snp_profiles <- vapply(bc$bc_72_snps$snp.freq, function(f) {
  as.numeric(stats::rbinom(n_snp, 2, f))
}, numeric(n_snp))
snp_profiles[stats::runif(length(snp_profiles)) < 0.3] <- NA
snp_profiles <- as.data.frame(snp_profiles)
names(snp_profiles) <- bc$bc_72_snps$snp.name

n_split <- 2000L
split_profiles <- bc$ref_cov_dat[sample(nrow(bc$ref_cov_dat), n_split, TRUE), ]

cohort <- bc$validation.cohort.data
cohort <- cohort[cohort$study.exit.age < 85, ]
ncc <- bc$validation.nested.case.control.data
cohort$inclusion <- as.numeric(cohort$id %in% ncc$id)
cohort$observed.followup <- cohort$study.exit.age - cohort$study.entry.age
selection <- stats::glm(inclusion ~ observed.outcome * (study.entry.age + observed.followup),
                        data = cohort, family = stats::binomial())
ncc$sampling.weights <- unname(stats::fitted(selection)[match(ncc$id, cohort$id)])
ncc <- ncc[!is.na(ncc$sampling.weights) & ncc$study.exit.age < 85, ]
covariates <- all.vars(bc$bc_model_formula)[-1]
cohort_10k <- cohort[seq_len(10000L), ]

# ----- Helpers ---------------------------------------------------------------
quietly <- function(expr) {
  result <- NULL
  utils::capture.output(result <- suppressMessages(suppressWarnings(expr)))
  result
}
with_c <- function(args) c(args, if (!is.null(use_c_code)) list(use.c.code = use_c_code))
model_args <- function(...) {
  utils::modifyList(list(
    model.formula = bc$bc_model_formula, model.cov.info = bc$bc_model_cov_info,
    model.log.RR = bc$bc_model_log_or, model.ref.dataset = bc$ref_cov_dat,
    model.disease.incidence.rates = bc$bc_inc,
    model.competing.incidence.rates = bc$mort_inc
  ), list(...))
}
absolute_risk <- function(...) quietly(do.call(computeAbsoluteRisk, with_c(list(...))))

# Risks used as inputs of the validation workloads (not timed)
five_year <- pmin(5, cohort$study.exit.age - cohort$study.entry.age)
cohort_risk <- do.call(absolute_risk, c(model_args(), list(
  apply.age.start = cohort$study.entry.age, apply.age.interval.length = five_year,
  apply.cov.profile = cohort[, covariates], return.lp = TRUE)))
ncc_followup <- ncc$study.exit.age - ncc$study.entry.age
ncc_risk <- do.call(absolute_risk, c(model_args(), list(
  apply.age.start = ncc$study.entry.age, apply.age.interval.length = ncc_followup,
  apply.cov.profile = ncc[, covariates], return.lp = TRUE)))
model_object <- c(model_args(), list(
  model.snp.info = NULL, model.ref.dataset.weights = NULL, model.bin.fh.name = NA,
  apply.cov.profile = cohort_10k[, covariates], apply.snp.profile = NULL, n.imp = 5,
  use.c.code = if (is.null(use_c_code)) 1 else use_c_code, return.lp = TRUE,
  return.refs.risk = FALSE))

# ----- Workloads ---------------------------------------------------------------
workloads <- list(
  risk_complete = list(
    description = "computeAbsoluteRisk: 10,000 complete covariate profiles",
    run = function() do.call(absolute_risk, c(model_args(), list(
      apply.age.start = complete_ages, apply.age.interval.length = complete_lengths,
      apply.cov.profile = complete_profiles)))),
  risk_missing = list(
    description = "computeAbsoluteRisk: 2,000 profiles with missing covariates",
    run = function() do.call(absolute_risk, c(model_args(), list(
      apply.age.start = missing_ages, apply.age.interval.length = missing_lengths,
      apply.cov.profile = missing_profiles)))),
  risk_snps = list(
    description = "computeAbsoluteRisk: covariates and 72 SNPs, 1,000 profiles (30% of genotypes missing)",
    run = function() {
      set.seed(1)
      do.call(absolute_risk, c(model_args(model.snp.info = bc$bc_72_snps,
                                          model.bin.fh.name = "famhist"), list(
        apply.age.start = 50, apply.age.interval.length = 20,
        apply.cov.profile = snp_cov_profiles, apply.snp.profile = snp_profiles)))
    }),
  split_interval = list(
    description = "computeAbsoluteRiskSplitInterval: 2,000 profiles, model changing at 50",
    run = function() quietly(do.call(computeAbsoluteRiskSplitInterval, with_c(c(model_args(), list(
      model.log.RR.2 = bc$bc_model_log_or_post_50, model.ref.dataset.2 = bc$ref_cov_dat_post_50,
      apply.age.start = 40, apply.age.interval.length = 30, cut.time = 50,
      apply.cov.profile = split_profiles)))))),
  validation_cohort = list(
    description = sprintf("ModelValidation: full cohort (%s subjects), risks supplied",
                          format(nrow(cohort), big.mark = ",")),
    run = function() quietly(ModelValidation(
      study.data = cohort, predicted.risk.interval = 5,
      predicted.risk = as.numeric(cohort_risk$risk),
      linear.predictor = as.numeric(cohort_risk$lps)))),
  validation_ncc = list(
    description = sprintf("ModelValidation: nested case-control study (%s subjects), risks supplied",
                          format(nrow(ncc), big.mark = ",")),
    run = function() quietly(ModelValidation(
      study.data = ncc, total.followup.validation = TRUE,
      predicted.risk = as.numeric(ncc_risk$risk),
      linear.predictor = as.numeric(ncc_risk$lps)))),
  validation_internal = list(
    description = "ModelValidation: 10,000 cohort subjects, risks computed from the model",
    run = function() quietly(ModelValidation(
      study.data = cohort_10k, predicted.risk.interval = 5,
      iCARE.model.object = model_object)))
)

# ----- Timing -------------------------------------------------------------------
time_workload <- function(run) {
  gc()
  memory <- tryCatch(
    as.numeric(bench::mark(run(), iterations = 1, check = FALSE, filter_gc = FALSE,
                           memory = TRUE)$mem_alloc),
    error = function(e) NA_real_)
  first <- system.time(run())[["elapsed"]]
  times <- if (first > 20) first else c(first, replicate(2, system.time(run())[["elapsed"]]))
  c(seconds = stats::median(times), memory_mb = memory / 2^20)
}

results <- do.call(rbind, lapply(names(workloads), function(name) {
  message("  ", name, " ...")
  outcome <- tryCatch(time_workload(workloads[[name]]$run),
                      error = function(e) {
                        message("    failed: ", conditionMessage(e))
                        c(seconds = NA_real_, memory_mb = NA_real_)
                      })
  data.frame(label = label, version = version,
             use_c_code = if (is.null(use_c_code)) NA else use_c_code,
             workload = name, description = workloads[[name]]$description,
             seconds = outcome[["seconds"]], memory_mb = outcome[["memory_mb"]])
}))
print(results[, c("workload", "seconds", "memory_mb")], row.names = FALSE)
utils::write.csv(results, file.path("tools", "benchmarks", paste0("results-", label, ".csv")),
                 row.names = FALSE)
