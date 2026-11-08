# S3 methods for the results of computeAbsoluteRisk(),
# computeAbsoluteRiskSplitInterval() and ModelValidation()

#' Print, summarize and plot absolute risk results
#'
#' Methods for the results of [computeAbsoluteRisk()] (class `"icare"`) and
#' [computeAbsoluteRiskSplitInterval()] (class `"icareSplit"`). `print()`
#' shows the risk interval and the risk of the first profiles, `summary()`
#' summarizes the risks, and `plot()` draws their distribution.
#'
#' @param x,object Result of [computeAbsoluteRisk()] or
#'   [computeAbsoluteRiskSplitInterval()].
#' @param n Number of profiles to print.
#' @param base.size Base font size of the plot, in points.
#' @param ... Passed to [summary()] for `summary()`; ignored otherwise.
#' @return `print()` returns `x` invisibly; `summary()` a summary table of
#'   the risks; `plot()` a ggplot object, invisibly, after drawing it.
#' @name icare-methods
#' @examples
#' data(bc_data, package = "iCARE")
#' results <- computeAbsoluteRisk(
#'   model.formula = bc_model_formula,
#'   model.cov.info = bc_model_cov_info,
#'   model.log.RR = bc_model_log_or,
#'   model.ref.dataset = ref_cov_dat,
#'   model.disease.incidence.rates = bc_inc,
#'   model.competing.incidence.rates = mort_inc,
#'   apply.age.start = 50,
#'   apply.age.interval.length = 30,
#'   apply.cov.profile = ref_cov_dat[1:100, ]
#' )
#' results
#' summary(results)
#' plot(results)
NULL

#' @rdname icare-methods
#' @export
print.icare <- function(x, n = 10, ...) {
  print_risk_details(x, n, c("Int_Start", "Int_End", "Risk_Estimate"))
}

#' @rdname icare-methods
#' @export
print.icareSplit <- function(x, n = 10, ...) {
  columns <- c("Int_Start", if ("Cut_Age" %in% names(x$details)) "Cut_Age", "Int_End",
               "Risk_Estimate")
  print_risk_details(x, n, columns)
}

print_risk_details <- function(x, n, columns) {
  details <- as.data.frame(x$details)
  n_profiles <- nrow(details)
  cat("Absolute risks of ", format_count(n_profiles), " profile",
      if (n_profiles != 1L) "s", "\n\n", sep = "")
  shown <- details[seq_len(min(n, n_profiles)), intersect(columns, names(details)),
                   drop = FALSE]
  rownames(shown) <- NULL
  print(shown)
  if (n_profiles > n) {
    cat("... ", format_count(n_profiles - n), " more profile",
        if (n_profiles - n != 1L) "s", "\n", sep = "")
  }
  cat("Profiles and risks: `$details`; risks: `$risk`.\n")
  invisible(x)
}

#' @rdname icare-methods
#' @export
summary.icare <- function(object, ...) {
  summary(object$risk, ...)
}

#' @rdname icare-methods
#' @export
summary.icareSplit <- function(object, ...) {
  summary(object$risk, ...)
}

#' @rdname icare-methods
#' @export
plot.icare <- function(x, base.size = 11, ...) {
  plot_risk_distribution(x, base.size)
}

#' @rdname icare-methods
#' @export
plot.icareSplit <- function(x, base.size = 11, ...) {
  plot_risk_distribution(x, base.size)
}

plot_risk_distribution <- function(x, base.size) {
  risk <- 100 * as.numeric(x$risk)
  risk <- risk[is.finite(risk)]
  if (length(risk) == 0L) {
    abort("There are no risks to plot.")
  }
  bins <- max(10L, min(50L, ceiling(sqrt(length(risk)))))
  figure <- ggplot2::ggplot(data.frame(risk = risk), ggplot2::aes(x = .data$risk)) +
    ggplot2::geom_histogram(bins = bins, fill = plot_colours[["primary"]],
                            colour = "white", linewidth = 0.3) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.05))) +
    ggplot2::labs(title = "Distribution of absolute risks", x = "Absolute risk (%)",
                  y = "Number of profiles") +
    plot_theme(base.size) +
    ggplot2::theme(aspect.ratio = 0.62, panel.grid.major.x = ggplot2::element_blank())
  print(figure)
  invisible(figure)
}

#' Print and plot model validation results
#'
#' `print()` reports the study, calibration in the large (E/O), the
#' Hosmer-Lemeshow and relative-risk tests, the AUC and the results of each
#' risk-score group of a [ModelValidation()] result. `plot()` draws them with
#' [plotModelValidation()].
#'
#' @param x Result of [ModelValidation()].
#' @param digits Number of decimals of ratios and the AUC.
#' @param ... Passed to [plotModelValidation()] by `plot()`; ignored by
#'   `print()`.
#' @return `print()` returns `x` invisibly; `plot()` returns the figure as
#'   [plotModelValidation()] does.
#' @name icareValid-methods
NULL

#' @rdname icareValid-methods
#' @export
print.icareValid <- function(x, digits = 3, ...) {
  args <- x$input.args
  study <- args$study.data
  outcome <- x$Subject_Specific_Observed_Outcome
  weights <- if (is.null(study$sampling.weights)) NULL else 1 / study$sampling.weights
  formula_text <- if (is.null(args$model.formula)) {
    "not given (risks supplied, or a SNP-only model)"
  } else {
    paste(deparse(args$model.formula, width.cutoff = 500L), collapse = " ")
  }
  line <- function(label, value) {
    cat(formatC(label, width = -27), value, "\n", sep = "")
  }
  range_text <- function(v, digits) {
    sprintf("mean %s (range %s-%s)", format_number(mean(v), digits),
            format(min(v)), format(max(v)))
  }

  cat("iCARE model validation\n\n")
  line("Dataset:", args$dataset)
  line("Model:", args$model.name)
  line("Model formula:", paste(strwrap(formula_text, width = 70, exdent = 27),
                               collapse = "\n"))
  line("Risk prediction interval:", x$Risk_Prediction_Interval)
  line("Subjects (cases):", paste0(format_count(length(outcome)), " (",
                                   format_count(sum(outcome)), ")"))
  if (!is.null(weights)) {
    line("Weighted subjects (cases):", paste0(format_count(sum(weights)), " (",
                                              format_count(sum(weights * outcome)), ")"))
  }
  line("Follow-up (years):", range_text(x$Adjusted_Followup, 2))
  if (!is.null(study$study.entry.age)) {
    line("Age at entry (years):", range_text(study$study.entry.age, 1))
  }

  test_text <- function(test) {
    statistic <- unname(test$statistic)
    if (length(statistic) != 1L || is.na(statistic)) {
      return("not computed")
    }
    paste0("chi-square ", format_number(statistic, 2), ", df ", test$parameter,
           ", p ", format_p_value(test$p.value))
  }
  cat("\nCalibration\n")
  line("  E/O:", format_estimate(x$Overall_Expected_to_Observed_Ratio,
                                 x$CI_Overall_Expected_to_Observed_Ratio, digits))
  hl_variance <- args$hosmer.lemeshow.variance
  line("  Hosmer-Lemeshow test:", paste0(
    test_text(x$Hosmer_Lemeshow_Results),
    if (!is.null(hl_variance)) paste0(" (", hl_variance, "-risk variance)")))
  line("  Relative-risk test:", test_text(x$RR_test_result))
  cat("\nDiscrimination\n")
  line("  AUC:", format_estimate(x$AUC, x$CI_AUC, digits))

  categories <- x$Category_Results
  if (!is.null(categories$Expected_to_Observed_Ratio)) {
    percent <- function(v) format_number(100 * v, 2)
    table <- data.frame(
      Group = categories$Categories,
      Subjects = format_count(categories$Number_of_Subjects),
      Cases = format_count(categories$Number_of_Cases),
      `Expected %` = percent(categories$Predicted_Absolute_Risk),
      `Observed % (95% CI)` = paste0(
        percent(categories$Observed_Absolute_Risk), " (",
        percent(categories$CI_Absolute_Risk_Lower), " to ",
        percent(categories$CI_Absolute_Risk_Upper), ")"),
      `E/O (95% CI)` = vapply(seq_len(nrow(categories)), function(g) {
        format_estimate(categories$Expected_to_Observed_Ratio[g],
                        c(categories$CI_Expected_to_Observed_Ratio_Lower[g],
                          categories$CI_Expected_to_Observed_Ratio_Upper[g]), 2)
      }, character(1)),
      check.names = FALSE
    )
    table$`E/O (95% CI)` <- sub("95% CI ", "", table$`E/O (95% CI)`, fixed = TRUE)
    cat("\nRisk-score groups\n")
    print(table, right = TRUE, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname icareValid-methods
#' @export
plot.icareValid <- function(x, ...) {
  plotModelValidation(x$input.args$study.data, x, ...)
}
