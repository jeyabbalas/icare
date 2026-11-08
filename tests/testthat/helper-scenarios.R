# Inputs of the regression scenarios
#
# The same inputs are used by the tests and by tools/fixtures/reference-outputs.R,
# which records the outputs of earlier versions of iCARE in
# tests/testthat/fixtures/. This file uses base R and the example data only.

load_example_data <- function() {
  env <- new.env()
  utils::data("bc_data", package = "iCARE", envir = env)
  env
}

scenario_inputs <- function(bc) {
  # Complete covariate profiles with various ages and interval lengths
  set.seed(1)
  rows <- sample(nrow(bc$ref_cov_dat), 200)
  ages <- sample(30:65, 200, replace = TRUE)
  lengths <- sample(1:15, 200, replace = TRUE)
  # Non-uniform reference weights
  set.seed(2)
  weights <- stats::runif(nrow(bc$ref_cov_dat), 0.5, 3)

  # A random sample of the validation cohort (exit before age 85, the last age
  # of the example incidence rates)
  set.seed(50)
  cohort_rows <- sample(nrow(bc$validation.cohort.data), 5000)
  cohort_rows <- cohort_rows[bc$validation.cohort.data$study.exit.age[cohort_rows] < 85]

  # Inclusion probabilities of the nested case-control study (as in the
  # vignette of earlier versions)
  cohort <- bc$validation.cohort.data
  cohort$inclusion <- as.numeric(cohort$id %in% bc$validation.nested.case.control.data$id)
  cohort$observed.followup <- cohort$study.exit.age - cohort$study.entry.age
  selection <- stats::glm(inclusion ~ observed.outcome * (study.entry.age + observed.followup),
                          data = cohort, family = stats::binomial())
  ncc <- bc$validation.nested.case.control.data
  ncc$sampling.weights <- unname(stats::fitted(selection)[match(ncc$id, cohort$id)])
  ncc <- ncc[ncc$study.exit.age < 85, ]

  list(profile_rows = rows, ages = ages, lengths = lengths, weights = weights,
       cohort = bc$validation.cohort.data[cohort_rows, ], ncc = ncc,
       covariates = all.vars(bc$bc_model_formula)[-1])
}

# Replace (or, with NULL, remove) named arguments; unlike modifyList(), data
# frames and other lists are replaced whole
override_args <- function(args, changes) {
  for (name in names(changes)) args[[name]] <- changes[[name]]
  args
}

# Arguments of the covariate model of the examples
covariate_model_args <- function(bc, ...) {
  override_args(list(
    model.formula = bc$bc_model_formula,
    model.cov.info = bc$bc_model_cov_info,
    model.log.RR = bc$bc_model_log_or,
    model.ref.dataset = bc$ref_cov_dat,
    model.disease.incidence.rates = bc$bc_inc,
    model.competing.incidence.rates = bc$mort_inc
  ), list(...))
}

# Arguments of the covariate and SNP model of the examples
snp_model_args <- function(bc, ...) {
  override_args(covariate_model_args(
    bc, model.snp.info = bc$bc_72_snps, model.bin.fh.name = "famhist"
  ), list(...))
}

# Risk scores with and without ties, weighted and unweighted, to compare the
# risk-score groups with those of earlier versions
group_cases <- function() {
  set.seed(2024)
  lapply(seq_len(120), function(i) {
    n <- sample(c(20, 57, 200, 1000), 1)
    x <- switch(i %% 6 + 1,
                stats::rnorm(n),
                round(stats::rnorm(n), 1),
                sample(c(0, 0.3, 0.7, 1.2), n, TRUE, prob = c(0.6, 0.2, 0.15, 0.05)),
                c(rep(0, n %/% 2), stats::rnorm(n - n %/% 2)),
                round(stats::rexp(n), 2),
                log(stats::runif(n)) / 3)
    weights <- if (i %% 2 == 0) NULL else 1 / stats::runif(n, 0.03, 1)
    list(x = x, weights = weights, q = sample(c(4, 5, 10, 20), 1))
  })
}

# Every 10th element, to keep the recorded reference risks small
thin <- function(x) if (is.null(x)) NULL else x[seq(1, length(x), by = 10)]
