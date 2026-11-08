# Study data and follow-up for model validation
# =============================================
#
# Subject i enters the study at age a_i, leaves it at age b_i and, if a case,
# develops the disease t_i years after entry (t_i = Inf for non-cases). The
# risk being validated covers tau_i years from entry: `predicted.risk.interval`,
# or the observed follow-up b_i - a_i when `total.followup.validation = TRUE`.
#
# Each subject is followed for min(tau_i, b_i - a_i) years and counts as a case
# only if the disease develops within that time. An onset recorded after exit
# is a data error; such a subject is censored at exit (a non-case).

#' Validate `study.data`
#'
#' @return A list with the numeric columns `outcome`, `entry`, `exit`, `onset`
#'   and `sampling_prob` (NULL for a full cohort).
#' @noRd
validate_study_data <- function(study.data) {
  data <- tryCatch(as_plain_data_frame(study.data, "study.data"), error = function(e) {
    abort("`study.data` must be a data frame.")
  })
  required <- c("observed.outcome", "study.entry.age", "study.exit.age", "time.of.onset")
  absent <- setdiff(required, names(data))
  if (length(absent) > 0L) {
    abort("`study.data` must have the columns ", format_names(required), "; missing ",
          format_names(absent), ".")
  }
  if (nrow(data) == 0L) {
    abort("`study.data` has no rows.")
  }
  numeric_column <- function(name) {
    x <- data[[name]]
    if (is.logical(x)) x <- as.numeric(x)
    if (!is.numeric(x) || anyNA(x)) {
      abort("`study.data$", name, "` must be numeric with no missing values.")
    }
    as.numeric(x)
  }
  outcome <- numeric_column("observed.outcome")
  if (!all(outcome %in% c(0, 1))) {
    abort("`study.data$observed.outcome` must be 1 for cases and 0 for non-cases.")
  }
  entry <- numeric_column("study.entry.age")
  exit <- numeric_column("study.exit.age")
  if (any(!is.finite(entry)) || any(!is.finite(exit))) {
    abort("`study.data$study.entry.age` and `study.data$study.exit.age` must be finite.")
  }
  if (any(entry >= exit)) {
    abort("`study.data$study.entry.age` must be lower than `study.data$study.exit.age` ",
          "for every subject; it is not in ", sum(entry >= exit), " row(s), the first ",
          "being row ", which(entry >= exit)[1], ".")
  }
  onset <- data[["time.of.onset"]]
  if (is.logical(onset)) onset <- as.numeric(onset)
  if (!is.numeric(onset)) {
    abort("`study.data$time.of.onset` must be numeric.")
  }
  onset <- as.numeric(onset)
  if (any(outcome == 1 & !is.finite(onset))) {
    abort("Every case (`observed.outcome` = 1) needs a finite `time.of.onset`.")
  }
  # The onset of a non-case plays no role (it may be Inf or NA).
  onset[outcome == 0] <- Inf
  if (any(onset < 0)) {
    abort("`study.data$time.of.onset` must be the time from study entry to disease ",
          "onset (non-negative) for cases.")
  }

  sampling_prob <- NULL
  if (!is.null(data[["sampling.weights"]])) {
    sampling_prob <- data[["sampling.weights"]]
    if (!is.numeric(sampling_prob) || anyNA(sampling_prob) ||
        any(sampling_prob <= 0 | sampling_prob > 1)) {
      abort("`study.data$sampling.weights` must hold each subject's probability of ",
            "inclusion in the sample, between 0 (excluded) and 1.")
    }
    sampling_prob <- as.numeric(sampling_prob)
  }
  list(outcome = outcome, entry = entry, exit = exit, onset = onset,
       sampling_prob = sampling_prob)
}

#' Follow-up and outcome over the risk prediction interval
#'
#' @param study Output of validate_study_data().
#' @param interval Length of the risk prediction interval per subject, or NULL
#'   for the observed follow-up.
#' @return A list with the adjusted `outcome` and `followup`, the `onset`
#'   (Inf for subjects censored at exit) and the `observed_followup`.
#' @noRd
adjust_validation_followup <- function(study, interval) {
  observed_followup <- study$exit - study$entry
  outcome <- study$outcome
  onset <- study$onset
  # Onsets are compared with the follow-up up to rounding error: with ages
  # recorded as decimals, exit - entry can fall just below the recorded onset
  # time of a case diagnosed at exit (70.94 - 61.71 < 9.23 in floating point),
  # and earlier versions turned such cases into non-cases.
  late <- outcome == 1 & onset > observed_followup + age_snap_tolerance
  if (any(late)) {
    warn(sum(late), " case(s) have `time.of.onset` after `study.exit.age`; they are ",
         "treated as non-cases censored at exit.")
    # Reset the onset as well as the outcome (earlier versions kept the onset).
    onset[late] <- Inf
    outcome[late] <- 0
  }
  at_exit <- outcome == 1 & onset > observed_followup
  onset[at_exit] <- observed_followup[at_exit]
  if (is.null(interval)) interval <- observed_followup
  followup <- pmin(interval, observed_followup)
  # Disease after the end of the risk prediction interval is not an event; an
  # onset at its end (up to rounding error) is.
  outcome[onset > interval + age_snap_tolerance] <- 0
  list(outcome = outcome, onset = onset, followup = followup,
       observed_followup = observed_followup)
}

#' Validate `predicted.risk.interval`
#' @return One interval length per subject, or NULL for the observed follow-up.
#' @noRd
normalize_prediction_interval <- function(interval, total_followup, n) {
  if (total_followup) {
    return(NULL)
  }
  if (is.null(interval)) {
    abort("Supply `predicted.risk.interval` (e.g. 5 to validate 5-year risks), or set ",
          "`total.followup.validation = TRUE` to validate risks over each subject's ",
          "observed follow-up.")
  }
  interval <- as.vector(interval)
  if (!is.numeric(interval) || anyNA(interval) || any(!is.finite(interval)) ||
      any(interval <= 0)) {
    abort("`predicted.risk.interval` must hold positive numbers of years.")
  }
  if (length(interval) == 1L) interval <- rep(interval, n)
  if (length(interval) != n) {
    abort("`predicted.risk.interval` must have length 1 or one value per row of ",
          "`study.data` (", n, ").")
  }
  as.numeric(interval)
}

#' Readable form of the interval label, for reports and figures
#' @noRd
display_prediction_interval <- function(label) {
  if (identical(label, "Observed Followup")) "observed follow-up" else label
}

#' Label of the risk prediction interval
#' @noRd
describe_prediction_interval <- function(interval) {
  if (is.null(interval)) {
    return("Observed Followup")
  }
  values <- unique(interval)
  if (length(values) == 1L) {
    return(paste(format(values), "years"))
  }
  paste0("Subject-specific (", format(min(values)), " to ", format(max(values)), " years)")
}
