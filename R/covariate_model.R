# Covariate part of the relative-risk model
# =========================================
#
# The relative risk of a covariate profile z is exp(x(z)' beta), where x(z) is
# the row of the design matrix that `model.formula` generates for z (without
# the intercept, which is absorbed into the baseline hazard) and beta is
# `model.log.RR`.
#
# `model.cov.info` declares each covariate as "continuous" or "factor"; factor
# covariates get the declared levels, with the optional referent level `ref`
# first, and are coded with treatment contrasts. The design is set up once on
# the reference dataset; its terms (including the data-dependent bases of
# functions such as poly() or ns(), through their "predvars") and factor levels
# are then reused for every other data set, so that profiles and the reference
# population are always coded identically. Missing covariate values in a
# profile give NA in exactly the design columns that depend on them.

#' Validate `model.cov.info` against the formula
#'
#' Each element is a list with `name` and `type` ("continuous" or "factor");
#' factors need `levels` and may have `ref`. As in earlier versions, `name` is
#' matched partially, so an element spelled `names` is accepted.
#'
#' @return A list of normalized entries named by covariate, in the order of
#'   `all.vars()` of the formula's right-hand side: `name`, `type`, `levels`
#'   (character, referent level first; NULL for continuous covariates).
#' @noRd
normalize_cov_info <- function(cov_info, formula) {
  if (!is.list(cov_info) || is.data.frame(cov_info) || length(cov_info) == 0L) {
    abort("`model.cov.info` must be a list with one list per covariate.")
  }
  entries <- lapply(seq_along(cov_info), function(i) {
    entry <- cov_info[[i]]
    if (!is.list(entry)) {
      abort("Element ", i, " of `model.cov.info` must be a list with fields ",
            "'name' and 'type' (and 'levels' for factors).")
    }
    name <- entry[["name", exact = FALSE]]
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
      abort("Element ", i, " of `model.cov.info` needs a 'name' field holding ",
            "the covariate name.")
    }
    type <- entry[["type", exact = FALSE]]
    if (!is.character(type) || length(type) != 1L || !type %in% c("continuous", "factor")) {
      abort("The 'type' of covariate '", name, "' in `model.cov.info` must be ",
            "\"continuous\" or \"factor\".")
    }
    levels <- entry[["levels", exact = FALSE]]
    if (type == "continuous") {
      if (!is.null(levels)) {
        warn("Covariate '", name, "' is declared continuous, so the 'levels' ",
             "given for it in `model.cov.info` are not used.")
      }
      return(list(name = name, type = type, levels = NULL))
    }
    if (is.null(levels) || !is.atomic(levels) || length(levels) == 0L || anyNA(levels)) {
      abort("Factor covariate '", name, "' needs a 'levels' field in ",
            "`model.cov.info`: a vector of its non-missing levels.")
    }
    levels <- as.character(levels)
    if (anyDuplicated(levels)) {
      abort("The levels of factor covariate '", name, "' in `model.cov.info` ",
            "must be distinct.")
    }
    ref <- entry[["ref", exact = FALSE]]
    if (!is.null(ref)) {
      ref <- as.character(ref)
      if (length(ref) != 1L || !ref %in% levels) {
        abort("The 'ref' of covariate '", name, "' in `model.cov.info` must be ",
              "one of its levels.")
      }
      levels <- c(ref, setdiff(levels, ref))
    }
    list(name = name, type = type, levels = levels)
  })
  names(entries) <- vapply(entries, `[[`, character(1), "name")
  if (anyDuplicated(names(entries))) {
    abort("`model.cov.info` describes the covariate(s) ",
          format_names(unique(names(entries)[duplicated(names(entries))])),
          " more than once.")
  }

  variables <- formula_covariates(formula)
  missing_info <- setdiff(variables, names(entries))
  unused_info <- setdiff(names(entries), variables)
  if (length(missing_info) > 0L || length(unused_info) > 0L) {
    abort("`model.cov.info` must describe exactly the covariates of ",
          "`model.formula`.",
          if (length(missing_info)) paste0("\n  Not described: ", format_names(missing_info)),
          if (length(unused_info)) paste0("\n  Not in the formula: ", format_names(unused_info)))
  }
  entries[variables]
}

# Covariates named on the right-hand side of the formula, in order.
formula_covariates <- function(formula) {
  check_model_formula(formula)
  all.vars(stats::delete.response(stats::terms(formula)))
}

check_model_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    abort("`model.formula` must be a formula such as `Y ~ x1 + x2`.")
  }
  rhs_terms <- stats::terms(formula)
  if (attr(rhs_terms, "intercept") != 1L) {
    abort("`model.formula` must keep the intercept (no `- 1` or `+ 0`): the ",
          "intercept is absorbed into the baseline hazard, and dropping it ",
          "would change how factors are coded.")
  }
  if (length(attr(rhs_terms, "term.labels")) == 0L) {
    abort("`model.formula` has no covariates; use `model.snp.info` alone for a ",
          "SNP-only model.")
  }
  invisible(TRUE)
}

#' Select and code the covariate columns of a data set
#'
#' Columns are selected by name, so extra columns are ignored. Factor
#' covariates become factors with the declared levels; continuous covariates
#' must be numeric (or logical).
#' @noRd
prepare_covariate_data <- function(data, cov_info, name, allow_missing) {
  data <- as_plain_data_frame(data, name)
  variables <- names(cov_info)
  absent <- setdiff(variables, names(data))
  if (length(absent) > 0L) {
    abort("`", name, "` lacks the covariate column(s) ", format_names(absent),
          ". Its columns must be named as the covariates of `model.formula`.")
  }
  out <- data[, variables, drop = FALSE]
  for (variable in variables) {
    entry <- cov_info[[variable]]
    x <- out[[variable]]
    if (entry$type == "factor") {
      values <- as.character(x)
      outside <- setdiff(unique(values[!is.na(values)]), entry$levels)
      if (length(outside) > 0L) {
        abort("Covariate '", variable, "' in `", name, "` has value(s) ",
              format_names(outside), " outside the levels declared in ",
              "`model.cov.info` (", format_names(entry$levels), ").")
      }
      out[[variable]] <- factor(values, levels = entry$levels)
    } else {
      if (!is.numeric(x) && !is.logical(x)) {
        abort("Covariate '", variable, "' is declared continuous in ",
              "`model.cov.info`, so it must be numeric in `", name, "`, not ",
              class(x)[1], ".")
      }
      x <- as.numeric(x)
      if (any(is.infinite(x))) {
        abort("Covariate '", variable, "' in `", name, "` has infinite values.")
      }
      out[[variable]] <- x
    }
  }
  if (!allow_missing && anyNA(out)) {
    bad <- variables[vapply(out, anyNA, logical(1))]
    abort("`", name, "` cannot contain missing values; missing in ",
          format_names(bad), ".")
  }
  rownames(out) <- NULL
  out
}

#' Set up the covariate design on the reference data
#'
#' @return A list with the design `spec` (terms with predvars, factor levels,
#'   contrasts, column names and the covariates each column depends on) and
#'   the reference `design` matrix.
#' @noRd
build_design_spec <- function(formula, reference) {
  rhs <- stats::delete.response(stats::terms(formula))
  # (The covariates are complete; missing values here would come from the
  # formula's transformations and are reported below.)
  frame <- stats::model.frame(rhs, data = reference, na.action = stats::na.pass)
  model_terms <- attr(frame, "terms")
  factor_columns <- names(frame)[vapply(frame, is.factor, logical(1))]
  contrasts <- if (length(factor_columns) > 0L) {
    stats::setNames(lapply(factor_columns, function(v) {
      if (is.ordered(frame[[v]])) "contr.poly" else "contr.treatment"
    }), factor_columns)
  } else {
    NULL
  }
  spec <- list(
    terms = model_terms,
    xlevels = stats::.getXlevels(model_terms, frame),
    contrasts = contrasts
  )
  full <- stats::model.matrix(model_terms, frame, contrasts.arg = contrasts)
  keep <- colnames(full) != "(Intercept)"
  design <- full[, keep, drop = FALSE]
  undefined <- rowSums(!is.finite(design)) > 0
  if (any(undefined)) {
    abort("`model.formula` gives infinite or undefined values (for example, ",
          "log(0)) in row(s) ", format_names(which(undefined)),
          " of `model.ref.dataset`.")
  }
  spec$columns <- colnames(design)
  spec$column_covariates <- column_covariates(model_terms, attr(full, "assign")[keep])
  dimnames(design) <- list(NULL, spec$columns)
  attr(design, "assign") <- NULL
  attr(design, "contrasts") <- NULL
  list(spec = spec, design = design)
}

# For every design column, the data covariates it is computed from.
column_covariates <- function(model_terms, assign) {
  factors <- attr(model_terms, "factors")
  frame_variables <- rownames(factors)
  covariates_of_frame_variable <- lapply(frame_variables, function(v) {
    all.vars(str2lang(v))
  })
  lapply(assign, function(term) {
    used <- factors[, term] > 0
    unique(unlist(covariates_of_frame_variable[used], use.names = FALSE))
  })
}

#' Design matrix of a prepared covariate data set
#'
#' Rows with missing covariates are kept, with NA in the columns that depend on
#' them.
#' @noRd
build_design_matrix <- function(spec, data) {
  if (nrow(data) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = length(spec$columns),
                  dimnames = list(NULL, spec$columns)))
  }
  frame <- stats::model.frame(spec$terms, data = data, na.action = stats::na.pass,
                              xlev = spec$xlevels)
  full <- stats::model.matrix(spec$terms, frame, contrasts.arg = spec$contrasts)
  design <- full[, colnames(full) != "(Intercept)", drop = FALSE]
  if (!identical(colnames(design), spec$columns)) {
    abort("Internal error: the design columns of a data set differ from those ",
          "of the reference dataset.")
  }
  dimnames(design) <- list(NULL, spec$columns)
  attr(design, "assign") <- NULL
  attr(design, "contrasts") <- NULL
  design
}

#' Log relative risks matched to the design columns
#'
#' Accepts a named numeric vector or a one-column matrix with row names. The
#' names must be exactly the design columns; they may come in any order.
#' @return A numeric vector named and ordered as `columns`.
#' @noRd
match_log_rr <- function(log_rr, columns, name = "model.log.RR") {
  if (is.data.frame(log_rr)) {
    abort("`", name, "` must be a named numeric vector or a one-column matrix ",
          "with row names, not a data frame.")
  }
  if (is.matrix(log_rr)) {
    if (ncol(log_rr) != 1L) {
      abort("`", name, "` must be a named vector or a one-column matrix.")
    }
    log_rr <- stats::setNames(as.vector(log_rr), rownames(log_rr))
  }
  if (!is.numeric(log_rr) || length(log_rr) == 0L) {
    abort("`", name, "` must be a named numeric vector of log relative risks.")
  }
  rr_names <- names(log_rr)
  design_text <- paste0("\n  Design columns: ", format_names(columns, 12L))
  if (is.null(rr_names) || anyNA(rr_names) || any(!nzchar(rr_names))) {
    abort("`", name, "` must be named after the design-matrix columns of ",
          "`model.formula`.", design_text)
  }
  if (anyDuplicated(rr_names)) {
    abort("`", name, "` has duplicated names: ",
          format_names(unique(rr_names[duplicated(rr_names)])), ".")
  }
  missing_names <- setdiff(columns, rr_names)
  extra_names <- setdiff(rr_names, columns)
  if (length(missing_names) > 0L || length(extra_names) > 0L) {
    abort("The names of `", name, "` must match the design-matrix columns ",
          "generated by `model.formula` and `model.cov.info`.",
          if (length(missing_names)) paste0("\n  No log relative risk for: ", format_names(missing_names, 12L)),
          if (length(extra_names)) paste0("\n  Not a design column: ", format_names(extra_names, 12L)),
          design_text)
  }
  log_rr <- log_rr[columns]
  if (any(!is.finite(log_rr))) {
    abort("`", name, "` must contain finite values.")
  }
  log_rr
}
