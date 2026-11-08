# Age-specific rate tables
# ========================
#
# iCARE works in whole years of age. The disease incidence rates lambda(t) and
# the competing mortality rates m(t) are annual hazards indexed by the integer
# age t, i.e. they apply to the year of age [t, t + 1).
#
# Users may give a table in either of two layouts:
#
#   * two columns  (age, rate): one row per integer age;
#   * three columns (start, end, rate): one row per age band. The rate is the
#     annual rate that applies to every year of age in the band. Bands are
#     read as inclusive integer ranges when consecutive rows satisfy
#     start[i + 1] == end[i] + 1 (e.g. 50-54, 55-59), and as half-open ranges
#     [start, end) when they share their end points, start[i + 1] == end[i]
#     (e.g. 50-55, 55-60), which is the layout older versions of iCARE
#     required. A single band is read as inclusive.
#
# Older versions divided a band's rate by the number of ages in it and let
# shared end points overwrite each other, which shrank the rates several-fold
# (the 5-year layout above gave each age a sixth of the rate).

#' Normalize one rate table to yearly rates (ages with a missing rate are dropped)
#'
#' @param rates A data frame, tibble or matrix with 2 or 3 numeric columns.
#' @param name Argument name used in error messages.
#' @param consecutive Require the yearly ages to be consecutive.
#' @return A data frame with integer-valued numeric columns `age` and `rate`,
#'   sorted by age.
#' @noRd
normalize_rate_table <- function(rates, name, consecutive = FALSE) {
  rates <- tryCatch(as_plain_data_frame(rates, name), error = function(e) {
    abort("`", name, "` must be a data frame or matrix with two columns ",
          "(age, rate) or three columns (start age, end age, rate).")
  })
  if (!ncol(rates) %in% c(2L, 3L)) {
    abort("`", name, "` must have two columns (age, rate) or three columns ",
          "(start age, end age, rate); it has ", ncol(rates), ".")
  }
  if (nrow(rates) == 0L) {
    abort("`", name, "` has no rows.")
  }
  numeric_columns <- vapply(rates, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    abort("All columns of `", name, "` must be numeric.")
  }
  age_columns <- as.matrix(rates[, -ncol(rates), drop = FALSE])
  if (anyNA(age_columns) || any(!is.finite(age_columns))) {
    abort("The ages in `", name, "` cannot be missing or infinite.")
  }
  rate <- rates[[ncol(rates)]]
  if (any(!is.finite(rate) & !is.na(rate))) {
    abort("The rates in `", name, "` cannot be infinite.")
  }
  if (any(rate < 0 | rate > 1, na.rm = TRUE)) {
    abort("The rates in `", name, "` must be annual rates between 0 and 1.")
  }

  yearly <- if (ncol(rates) == 2L) {
    yearly_from_ages(rates[[1]], rate, name)
  } else {
    yearly_from_bands(rates[[1]], rates[[2]], rate, name)
  }
  # Ages whose rate is missing are treated as not covered by the table.
  yearly <- yearly[!is.na(yearly$rate), , drop = FALSE]
  rownames(yearly) <- NULL
  if (nrow(yearly) == 0L) {
    abort("`", name, "` has no non-missing rates.")
  }

  if (consecutive && any(diff(yearly$age) != 1)) {
    gaps <- setdiff(seq(min(yearly$age), max(yearly$age)), yearly$age)
    abort("`", name, "` must cover consecutive ages, because the baseline ",
          "hazard is calibrated by accumulating it from the youngest age; ",
          "rates are missing for ages ", format_age_ranges(gaps), ".")
  }
  yearly
}

yearly_from_ages <- function(age, rate, name) {
  if (any(age %% 1 != 0)) {
    abort("The ages in the first column of `", name, "` must be whole numbers.")
  }
  if (anyDuplicated(age)) {
    abort("`", name, "` lists the age(s) ", format_age_ranges(age[duplicated(age)]),
          " more than once.")
  }
  ord <- order(age)
  data.frame(age = as.numeric(age[ord]), rate = as.numeric(rate[ord]))
}

yearly_from_bands <- function(start, end, rate, name) {
  if (any(start %% 1 != 0) || any(end %% 1 != 0)) {
    abort("The start and end ages of the bands in `", name, "` must be whole ",
          "numbers.")
  }
  if (any(end < start)) {
    abort("Every band in `", name, "` must end at or after its start age.")
  }
  ord <- order(start)
  start <- start[ord]
  end <- end[ord]
  rate <- rate[ord]
  n <- length(start)

  last_age <- end
  if (n > 1L) {
    if (all(start[-1] == end[-n] + 1)) {
      # Inclusive bands: start-end covers ages start, ..., end.
    } else if (all(start[-1] == end[-n]) && all(end > start)) {
      # Bands sharing end points are half-open: [start, end).
      last_age <- end - 1
    } else {
      abort("The age bands in `", name, "` must be contiguous and must not ",
            "overlap: either inclusive bands such as 50-54, 55-59 (each band ",
            "starts one year after the previous one ends) or bands sharing ",
            "their end points such as 50-55, 55-60 (read as [50, 55), [55, 60)).")
    }
  }
  widths <- last_age - start + 1
  ages <- unlist(Map(seq, start, last_age), use.names = FALSE)
  data.frame(age = as.numeric(ages), rate = as.numeric(rep(rate, widths)))
}

#' Resolve the disease and competing rate tables of a model
#'
#' The incidence table must cover consecutive ages (it defines the ages over
#' which the baseline hazard is calibrated). Without competing rates, the
#' competing hazard is zero at every incidence age.
#' @noRd
resolve_rate_tables <- function(incidence, competing) {
  if (is.null(incidence)) {
    abort("`model.disease.incidence.rates` must be supplied.")
  }
  incidence <- normalize_rate_table(incidence, "model.disease.incidence.rates",
                                    consecutive = TRUE)
  competing <- if (is_absent(competing)) {
    data.frame(age = incidence$age, rate = 0)
  } else {
    normalize_rate_table(competing, "model.competing.incidence.rates")
  }
  list(incidence = incidence, competing = competing)
}

#' Check that a rate table covers the ages needed by the risk windows
#' @param needed Integer ages that must be present.
#' @noRd
check_rate_coverage <- function(table, needed, name) {
  absent <- setdiff(needed, table$age)
  if (length(absent) > 0L) {
    abort("`", name, "` has no rate for the age(s) ", format_age_ranges(absent),
          ", which the prediction intervals need. Rates are required for every ",
          "year of age from floor(apply.age.start) to ",
          "ceiling(apply.age.start + apply.age.interval.length) - 1.")
  }
  invisible(TRUE)
}

# Rates of `table` at integer `ages`; ages absent from the table get `fill`.
rates_at <- function(table, ages, fill = NA_real_) {
  out <- table$rate[match(ages, table$age)]
  out[is.na(out)] <- fill
  out
}
