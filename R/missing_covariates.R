# Risks for profiles with missing risk factors
# ============================================
#
# iCARE's estimand for a profile with observed risk factors z_o and missing
# risk factors Z_m is the model-averaged risk (Pal Choudhury et al. 2020,
# equation 3)
#
#   AR(z_o) = E[ AR(z_o, Z_m) | z_o ],
#
# where the conditional distribution of Z_m given z_o is approximated
# non-parametrically from the reference population: reference rows are
# grouped into percentile bins (100 bins) of the linear predictor of the
# observed risk factors, and the rows in the profile's bin are the donors.
# The risk is the weighted mean, over donors, of the risk of the profile's
# OWN observed risk factors completed with each donor's values of the missing
# ones. (Earlier versions averaged the donors' own complete-data risks, which
# replaces the profile's observed risk factors by the donors' and misses the
# estimand whenever the bin is wide or a missing factor has no effect.)
#
# Missing SNPs are handled the same way: the reference population carries
# simulated genotypes, so donors supply genotypes drawn from the population
# distribution given the profile's observed risk factors. A profile with no
# observed risk factor gets the weighted mean over the whole population.
#
# The percentile bins, the rounding of their cut points and the widening of
# empty bins are those of earlier versions (get_cutpoints(), get_miss_perc(),
# get_newEndpoints()); bins are formed from all reference rows, unweighted, and
# a bin counts as empty when its donors have no positive weight.

missing_data_bins <- 100L
# Upper bound on (profile, donor) pairs evaluated in one vectorized pass.
pair_batch_size <- 2000000L

#' Cut points of the percentile bins of a reference linear predictor
#'
#' Quantiles (type 7) at `probs`, rounded down to a precision of 0.1% of the
#' standard deviation, de-duplicated, with Inf appended. Bins are the
#' half-open intervals from cut_k to cut_k+1.
#' @noRd
get_cutpoints <- function(ref_lp, probs) {
  se <- stats::sd(ref_lp, na.rm = TRUE)
  if (is.na(se) || se < 1e-12) se <- 1e-12
  precision <- 1 / (0.001 * se)
  digits <- sum(precision >= 10^(1:16))
  cuts <- stats::quantile(ref_lp, probs = probs, names = FALSE)
  cuts <- floor(cuts * 10^digits) / 10^digits
  cuts <- c(unique(cuts), Inf)
  if (length(cuts) < 2L) {
    abort("Internal error: could not compute the percentile cut points of the ",
          "reference linear predictor.")
  }
  cuts
}

#' Bin of each linear predictor (vectorized get_miss_perc())
#' @noRd
get_miss_perc <- function(cutpoints, lp) {
  ncuts <- length(cutpoints)
  bin <- findInterval(as.numeric(lp), cutpoints)
  bin[bin == ncuts] <- ncuts - 1L
  bin[bin == 0L] <- 1L
  bin
}

#' Widen an empty bin one cut point at a time on both sides
#'
#' @param has_donors function(low, high) telling whether the half-open range from low to high holds
#'   usable donors.
#' @return c(low, high).
#' @noRd
get_new_endpoints <- function(cutpoints, bin, has_donors) {
  ncuts <- length(cutpoints)
  left <- bin - 1L
  right <- bin + 2L
  for (i in seq_len(ncuts)) {
    left <- max(left, 1L)
    right <- min(right, ncuts)
    low <- cutpoints[left]
    high <- cutpoints[right]
    if (has_donors(low, high)) break
    left <- left - 1L
    right <- right + 1L
  }
  c(low, high)
}

#' Model-averaged risks for profiles with missing risk factors
#'
#' @param rows Indices (into the profile set) of the profiles to impute.
#' @param ctx Prediction context (see `predict_risks()`): profile design and
#'   missingness, population, coefficients, windows, rates.
#' @return `list(risk, lp, survival)` for `rows`.
#' @noRd
impute_missing_risks <- function(rows, ctx) {
  n_rows <- length(rows)
  out <- list(risk = numeric(n_rows), lp = numeric(n_rows), survival = numeric(n_rows))
  if (n_rows == 0L) {
    return(out)
  }
  present <- !is.na(ctx$design[rows, , drop = FALSE])
  missing_variables <- ctx$variable_missing[rows, , drop = FALSE]
  pattern <- paste(pattern_strings(present), pattern_strings(missing_variables))
  for (group in split(seq_len(n_rows), pattern)) {
    res <- impute_pattern(rows[group], present[group[1], ],
                          missing_variables[group[1], ], ctx)
    out$risk[group] <- res$risk
    out$lp[group] <- res$lp
    out$survival[group] <- res$survival
  }
  out
}

# One string per row of a logical matrix, e.g. "0110".
pattern_strings <- function(x) {
  if (ncol(x) == 0L) {
    return(rep("", nrow(x)))
  }
  do.call(paste0, as.data.frame(x * 1L))
}

#' Profiles sharing one pattern of observed design columns and missing variables
#' @noRd
impute_pattern <- function(rows, present, missing_variables, ctx) {
  population <- ctx$population
  beta <- ctx$beta
  present_columns <- which(present)
  n_population <- length(population$weights)

  # ----- Donors: the profile's percentile bin of the observed linear predictor
  if (length(present_columns) > 0L) {
    ref_lp <- as.vector(population$design[, present_columns, drop = FALSE] %*%
                          beta[present_columns])
    observed_lp <- as.vector(ctx$design[rows, present_columns, drop = FALSE] %*%
                               beta[present_columns])
    donor_ranges <- match_donor_ranges(ref_lp, observed_lp, population$weights)
  } else {
    observed_lp <- numeric(length(rows))
    donor_ranges <- list(order = seq_len(n_population),
                         low = rep(1L, length(rows)),
                         high = rep(n_population, length(rows)))
  }

  # ----- How the completed profile is formed column by column
  missing_set <- names(missing_variables)[missing_variables]
  absent_columns <- which(!present)
  mixed <- vapply(absent_columns, function(c) {
    vars <- ctx$column_variables[[c]]
    n_missing <- sum(vars %in% missing_set)
    n_missing > 0L && n_missing < length(vars)
  }, logical(1))
  donor_columns <- absent_columns[!mixed]
  mixed_columns <- absent_columns[mixed]
  donor_lp <- if (length(donor_columns) > 0L) {
    as.vector(population$design[, donor_columns, drop = FALSE] %*% beta[donor_columns])
  } else {
    numeric(n_population)
  }

  if (length(mixed_columns) == 0L) {
    impute_additive(rows, observed_lp, donor_lp, donor_ranges, ctx)
  } else {
    impute_with_recomputation(rows, observed_lp, donor_lp, donor_ranges,
                              mixed_columns, missing_set, ctx)
  }
}

#' Donor ranges (positions in the order of the reference linear predictor)
#' @noRd
match_donor_ranges <- function(ref_lp, observed_lp, weights) {
  ord <- order(ref_lp)
  sorted <- ref_lp[ord]
  cumulative_weight <- c(0, cumsum(weights[ord]))
  count_below <- function(x) findInterval(x, sorted, left.open = TRUE)
  has_donors <- function(low, high) {
    lo <- count_below(low) + 1L
    hi <- count_below(high)
    hi >= lo & (cumulative_weight[hi + 1L] - cumulative_weight[lo]) > 0
  }

  cuts <- get_cutpoints(ref_lp, seq(0, 1, 1 / missing_data_bins))
  bin <- get_miss_perc(cuts, observed_lp)
  low <- cuts[bin]
  high <- cuts[bin + 1L]
  for (i in which(!has_donors(low, high))) {
    widened <- get_new_endpoints(cuts, bin[i], has_donors)
    low[i] <- widened[1]
    high[i] <- widened[2]
  }
  if (any(!has_donors(low, high))) {
    abort("Internal error: no reference rows with positive weight matched a ",
          "profile with missing risk factors.")
  }
  list(order = ord, low = count_below(low) + 1L, high = count_below(high))
}

#' Additive completion: the completed linear predictor is the profile's
#' observed part plus the donor's part for the missing risk factors
#' @noRd
impute_additive <- function(rows, observed_lp, donor_lp, donor_ranges, ctx) {
  weights <- ctx$population$weights
  n_rows <- length(rows)

  # Profiles in the same bin share donors; within a bin, identical profiles
  # (same observed linear predictor and window) need one evaluation, and
  # donors with the same missing-part linear predictor are pooled.
  bin_key <- paste(donor_ranges$low, donor_ranges$high)
  pending <- new_pair_batch()
  flush <- function(batch) {
    if (length(batch$lp) == 0L) {
      return(invisible(NULL))
    }
    res <- compute_absolute_risk_engine(batch$lp, batch$start, batch$length,
                                        ctx$baseline, ctx$competing,
                                        survival = ctx$survival)
    pair_risk <- if (ctx$survival) res$risk else res
    target <- batch$target
    total_weight <- rowsum(batch$weight, target, reorder = TRUE)
    index <- as.integer(rownames(total_weight))
    risk_sum <- rowsum(batch$weight * pair_risk, target, reorder = TRUE)
    lp_sum <- rowsum(batch$weight * batch$lp, target, reorder = TRUE)
    unique_risk[index] <<- risk_sum[, 1] / total_weight[, 1]
    unique_lp[index] <<- lp_sum[, 1] / total_weight[, 1]
    if (ctx$survival) {
      surv_sum <- rowsum(batch$weight * res$survival, target, reorder = TRUE)
      unique_survival[index] <<- surv_sum[, 1] / total_weight[, 1]
    }
    invisible(NULL)
  }

  # Unique (bin, observed lp, window) combinations.
  combo_frame <- data.frame(bin = bin_key, lp = observed_lp,
                            start = ctx$start[rows], length = ctx$length[rows])
  combo_id <- match_rows(combo_frame)
  combos <- which(!duplicated(combo_id))
  unique_risk <- unique_lp <- unique_survival <- numeric(length(combos))
  combo_index <- match(combo_id, combo_id[combos])

  for (bin in unique(bin_key[combos])) {
    members <- combos[bin_key[combos] == bin]
    positions <- seq(donor_ranges$low[members[1]], donor_ranges$high[members[1]])
    donors <- donor_ranges$order[positions]
    donors <- donors[weights[donors] > 0]
    donor_values <- donor_lp[donors]
    distinct <- unique(donor_values)
    pooled_weight <- as.vector(rowsum(weights[donors], match(donor_values, distinct),
                                      reorder = TRUE))
    k <- length(distinct)
    target <- match(members, combos)
    pending <- add_to_pair_batch(
      pending,
      lp = rep(observed_lp[members], each = k) + rep(distinct, times = length(members)),
      start = rep(ctx$start[rows[members]], each = k),
      length = rep(ctx$length[rows[members]], each = k),
      weight = rep(pooled_weight, times = length(members)),
      target = rep(target, each = k)
    )
    if (length(pending$lp) >= pair_batch_size) {
      flush(pending)
      pending <- new_pair_batch()
    }
  }
  flush(pending)

  list(risk = unique_risk[combo_index], lp = unique_lp[combo_index],
       survival = if (ctx$survival) unique_survival[combo_index] else rep(NA_real_, n_rows))
}

#' Completion when some design columns combine observed and missing risk
#' factors (interactions): those columns are recomputed from the completed
#' covariate data
#' @noRd
impute_with_recomputation <- function(rows, observed_lp, donor_lp, donor_ranges,
                                      mixed_columns, missing_set, ctx) {
  population <- ctx$population
  weights <- population$weights
  beta <- ctx$beta
  missing_covariates <- intersect(missing_set, names(ctx$covariate_data))
  n_rows <- length(rows)
  risk <- lp <- survival <- numeric(n_rows)

  # Batches of profiles whose donor sets together stay below the batch size.
  sizes <- donor_ranges$high - donor_ranges$low + 1L
  batch_id <- (cumsum(sizes) - 1) %/% pair_batch_size
  for (batch in split(seq_len(n_rows), batch_id)) {
    donor_list <- lapply(batch, function(i) {
      d <- donor_ranges$order[seq(donor_ranges$low[i], donor_ranges$high[i])]
      d[weights[d] > 0]
    })
    counts <- lengths(donor_list)
    donors <- unlist(donor_list, use.names = FALSE)
    profile_index <- rep(batch, counts)
    completed <- ctx$covariate_data[rows[profile_index], , drop = FALSE]
    reference_rows <- population$reference_row[donors]
    for (v in missing_covariates) {
      completed[[v]] <- ctx$reference_data[[v]][reference_rows]
    }
    design <- build_design_matrix(ctx$covariate_spec, completed)
    mixed_names <- colnames(ctx$design)[mixed_columns]
    pair_lp <- observed_lp[profile_index] + donor_lp[donors] +
      as.vector(design[, mixed_names, drop = FALSE] %*% beta[mixed_columns])
    res <- compute_absolute_risk_engine(pair_lp, ctx$start[rows[profile_index]],
                                        ctx$length[rows[profile_index]],
                                        ctx$baseline, ctx$competing,
                                        survival = ctx$survival)
    pair_risk <- if (ctx$survival) res$risk else res
    w <- weights[donors]
    total_weight <- rowsum(w, profile_index, reorder = TRUE)[, 1]
    risk[batch] <- rowsum(w * pair_risk, profile_index, reorder = TRUE)[, 1] / total_weight
    lp[batch] <- rowsum(w * pair_lp, profile_index, reorder = TRUE)[, 1] / total_weight
    if (ctx$survival) {
      survival[batch] <- rowsum(w * res$survival, profile_index, reorder = TRUE)[, 1] /
        total_weight
    }
  }
  list(risk = risk, lp = lp,
       survival = if (ctx$survival) survival else rep(NA_real_, n_rows))
}

new_pair_batch <- function() {
  list(lp = numeric(0), start = numeric(0), length = numeric(0),
       weight = numeric(0), target = integer(0))
}

add_to_pair_batch <- function(batch, lp, start, length, weight, target) {
  list(lp = c(batch$lp, lp), start = c(batch$start, start),
       length = c(batch$length, length), weight = c(batch$weight, weight),
       target = c(batch$target, target))
}

# Integer id per row of a data frame; equal rows get equal ids.
match_rows <- function(frame) {
  key <- do.call(paste, c(lapply(frame, function(x) {
    if (is.numeric(x)) sprintf("%.17g", x) else as.character(x)
  }), sep = "\r"))
  match(key, unique(key))
}
