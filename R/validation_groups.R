# Risk-score groups for calibration
# =================================
#
# Subjects are grouped by their risk score (linear predictor), either at
# user-supplied cut-offs or at quantiles of the score. The quantile groups
# follow gtools::quantcut() (type-7 quantiles) for a full cohort and use
# weighted quantiles (inverse-probability weights, "i/n" interpolation of the
# weighted empirical distribution, as Hmisc::wtd.quantile()) for a sampled
# cohort, with the repairs of the earlier fork of iCARE: values tied at
# repeated quantiles form their own groups instead of stopping with "breaks
# are not unique", and no subject is ever left without a group. Weighted
# quantiles are computed on the exact score values (Hmisc rounds tied values
# to 15 significant digits, which could leave the extreme subjects out).

#' Weighted quantiles by linear interpolation of the weighted ECDF
#'
#' The ECDF takes the value F_k = (cumulative weight up to the k-th smallest
#' distinct value) / (total weight) at the k-th smallest distinct value, and
#' F = 0 at the smallest value; quantiles interpolate linearly between these
#' points.
#' @noRd
weighted_quantile_interpolated <- function(x, weights, probs) {
  keep <- weights > 0
  x <- x[keep]
  weights <- weights[keep]
  values <- sort(unique(x))
  totals <- as.vector(rowsum(weights, match(x, values), reorder = TRUE))
  cumulative <- cumsum(totals)
  # Dividing by the last cumulative sum makes the ECDF end at exactly 1, so
  # that the top quantile is exactly max(x).
  ecdf <- cumulative / cumulative[length(cumulative)]
  stats::approx(c(0, ecdf), c(values[1], values), xout = probs, rule = 2,
                ties = "ordered")$y
}

#' Group subjects by their risk score
#'
#' @param score Risk scores.
#' @param weights Inverse-probability weights, or NULL for a full cohort.
#' @param n_groups Number of quantile groups (`number.of.percentiles`).
#' @param cutoffs Optional cut-offs (`linear.predictor.cutoffs`).
#' @return A factor with one level per non-empty group, ordered by score.
#' @noRd
risk_score_groups <- function(score, weights = NULL, n_groups = 10L, cutoffs = NULL) {
  groups <- if (!is.null(cutoffs)) {
    if (!is.numeric(cutoffs) || anyNA(cutoffs)) {
      abort("`linear.predictor.cutoffs` must be numeric with no missing values.")
    }
    # Cut-offs outside the range of the scores bound empty groups. (Kept as
    # breaks, a cut-off below the minimum would put the subjects at the
    # minimum in a group of their own, as earlier versions did.)
    lowest <- min(score)
    highest <- max(score)
    inside <- sort(unique(cutoffs[cutoffs > lowest & cutoffs < highest]))
    if (lowest == highest) {
      factor(rep(format(lowest), length(score)))
    } else {
      cut(score, breaks = c(lowest, inside, highest), include.lowest = TRUE)
    }
  } else {
    quantile_groups(score, weights, n_groups)
  }
  groups <- droplevels(groups)
  if (anyNA(groups)) {
    abort("Internal error: some subjects were not assigned a risk-score group.")
  }
  if (nlevels(groups) < 2L) {
    abort("The risk scores formed ", nlevels(groups), " non-empty group; the ",
          "calibration tests need at least 2. Check `linear.predictor.cutoffs` and ",
          "`number.of.percentiles`, or whether `linear.predictor` takes a single value.")
  }
  groups
}

quantile_groups <- function(x, weights, n_groups) {
  probs <- seq(0, 1, length.out = n_groups + 1L)
  quant <- if (is.null(weights)) {
    stats::quantile(x, probs, names = FALSE)
  } else {
    weighted_quantile_interpolated(x, weights, probs)
  }
  if (!anyDuplicated(quant)) {
    return(cut(x, quant, include.lowest = TRUE))
  }

  # Repeated quantiles: values equal to a repeated quantile form their own
  # group "[x]"; the other values are cut at the distinct quantiles moved up
  # to the next observed value.
  tied <- x %in% unique(quant[duplicated(quant)])
  labels <- ifelse(tied, paste0("[", as.character(x), "]"), NA_character_)
  reposition <- function(q) {
    above <- x >= q
    if (any(above)) min(x[above]) else q
  }
  breaks <- unique(vapply(unique(quant), reposition, numeric(1)))
  if (any(!tied)) {
    labels[!tied] <- if (length(breaks) >= 2L) {
      as.character(cut(x[!tied], breaks = breaks, include.lowest = TRUE))
    } else {
      paste0("[", as.character(x[!tied]), "]")
    }
  }
  levels_in_order <- unique(labels[order(x)])
  groups <- factor(labels, levels = levels_in_order)
  levels(groups) <- tidy_group_labels(levels_in_order)
  groups
}

# Labels as gtools::quantcut() makes them: a tied-value group is labelled by
# its value, the first interval is closed, the others are open on the left
# (as cut() makes them), and an interval followed by a tied-value group at
# its upper end is open on the right, e.g. "(1,2)" followed by "2".
tidy_group_labels <- function(labels) {
  bounds <- lapply(strsplit(labels, "[^0-9eE+.-]+"), function(parts) {
    parts <- parts[nzchar(parts)]
    if (length(parts) == 1L) rep(parts, 2L) else parts[1:2]
  })
  lower <- vapply(bounds, `[`, character(1), 1L)
  upper <- vapply(bounds, `[`, character(1), 2L)
  point <- lower == upper
  n <- length(labels)
  closed_lower <- c(TRUE, rep(FALSE, n - 1L))
  closed_upper <- c(!(upper[-n] == lower[-1] & point[-1]), TRUE)
  ifelse(point, lower,
         paste0(ifelse(closed_lower, "[", "("), lower, ",", upper,
                ifelse(closed_upper, "]", ")")))
}
