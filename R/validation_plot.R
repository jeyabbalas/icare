#' Plot the results of a model validation
#'
#' @description
#' Draws a publication-ready figure of the results of [ModelValidation()]:
#' * **A. Absolute risk calibration**: observed against expected absolute
#'   risk (%) of each risk-score group, with 95% confidence intervals, E/O
#'   and the Hosmer-Lemeshow test;
#' * **B. Relative risk calibration**: observed against expected relative
#'   risks on log scales, with the relative-risk test;
#' * **C. Discrimination**: the distributions of the risk score in cases and
#'   non-cases, with the AUC;
#' * **D. Incidence by attained age**: the study's incidence rates, with 95%
#'   confidence intervals, against the population incidence rates of the
#'   model.
#'
#' An optional panel shows the E/O ratio of each risk-score group
#' (`panels = "eo"`). Every panel is square, and the calibration panels have
#' the same limits on both axes, so that the line of perfect calibration runs
#' at 45 degrees.
#'
#' @param study.data The `study.data` given to [ModelValidation()]; by default
#'   the copy stored in `validation.results`. Used for the sampling weights
#'   and the ages at entry.
#' @param validation.results Result of [ModelValidation()].
#' @param dataset,model.name Names shown in the subtitle and title; by
#'   default those given to [ModelValidation()].
#' @param x.lim.absrisk,y.lim.absrisk Axis limits (in %) of the absolute risk
#'   panel; by default the same for both axes, from 0.
#' @param x.lab.absrisk,y.lab.absrisk Axis labels of the absolute risk panel.
#' @param x.lim.RR,y.lim.RR Axis limits of the relative risk panel; by
#'   default the same for both axes.
#' @param x.lab.RR,y.lab.RR Axis labels of the relative risk panel.
#' @param risk.score.plot.kernel Kernel of the risk-score densities (see
#'   [stats::density()]).
#' @param risk.score.plot.bandwidth Bandwidth of the risk-score densities: a
#'   number, or a rule of [stats::density()] (`"nrd0"`, `"nrd"`, `"ucv"`,
#'   `"bcv"`, `"SJ"`). In sampled data, `"nrd0"` and `"nrd"` use the weighted
#'   standard deviation and interquartile range and the effective sample size.
#' @param risk.score.plot.percent.smooth Ignored; kept for compatibility with
#'   earlier versions.
#' @param panels Panels to draw, in order: any of `"absolute"`, `"relative"`,
#'   `"discrimination"`, `"incidence"` and `"eo"` (E/O by risk-score group).
#' @param base.size Base font size, in points.
#' @param draw If `TRUE` (default), draw the figure.
#'
#' @return The figure, a `patchwork` object (invisibly when `draw = TRUE`). Save
#'   it with [ggplot2::ggsave()], e.g. `ggsave("validation.pdf", figure,
#'   width = 8, height = 8.6)` for the four default panels.
#'
#' @seealso [ModelValidation()]
#'
#' @examples
#' data(bc_data, package = "iCARE")
#' set.seed(50)
#' cohort <- validation.cohort.data[sample(nrow(validation.cohort.data), 5000), ]
#' validation <- ModelValidation(
#'   study.data = cohort,
#'   predicted.risk.interval = 5,
#'   iCARE.model.object = list(
#'     model.formula = bc_model_formula,
#'     model.cov.info = bc_model_cov_info,
#'     model.log.RR = bc_model_log_or,
#'     model.ref.dataset = ref_cov_dat,
#'     model.disease.incidence.rates = bc_inc,
#'     model.competing.incidence.rates = mort_inc,
#'     apply.cov.profile = cohort[, all.vars(bc_model_formula)[-1]]
#'   ),
#'   dataset = "Simulated cohort",
#'   model.name = "Breast cancer risk model"
#' )
#' figure <- plotModelValidation(validation.results = validation)
#'
#' # E/O by risk-score group next to the absolute risk calibration
#' plotModelValidation(validation.results = validation, panels = c("absolute", "eo"))
#' @export
plotModelValidation <- function(study.data = NULL, validation.results,
                                dataset = NULL, model.name = NULL,
                                x.lim.absrisk = NULL, y.lim.absrisk = NULL,
                                x.lab.absrisk = "Expected absolute risk (%)",
                                y.lab.absrisk = "Observed absolute risk (%)",
                                x.lim.RR = NULL, y.lim.RR = NULL,
                                x.lab.RR = "Expected relative risk",
                                y.lab.RR = "Observed relative risk",
                                risk.score.plot.kernel = "gaussian",
                                risk.score.plot.bandwidth = "nrd0",
                                risk.score.plot.percent.smooth = 50,
                                panels = c("absolute", "relative", "discrimination",
                                           "incidence"),
                                base.size = 10, draw = TRUE) {
  if (!inherits(validation.results, "icareValid")) {
    abort("`validation.results` must be the result of ModelValidation().")
  }
  panels <- unique(match.arg(panels, c("absolute", "relative", "discrimination",
                                       "incidence", "eo"), several.ok = TRUE))
  if (!is.numeric(base.size) || length(base.size) != 1L || !is.finite(base.size) ||
      base.size <= 4) {
    abort("`base.size` must be a font size in points, greater than 4.")
  }
  draw <- as_flag(draw, "draw")
  check_axis_limits(x.lim.absrisk, "x.lim.absrisk")
  check_axis_limits(y.lim.absrisk, "y.lim.absrisk")
  check_axis_limits(x.lim.RR, "x.lim.RR", positive = TRUE)
  check_axis_limits(y.lim.RR, "y.lim.RR", positive = TRUE)

  data <- validation_plot_data(validation.results, study.data)
  dataset <- dataset %||% data$dataset
  model.name <- model.name %||% data$model.name

  plots <- lapply(panels, function(panel) {
    switch(panel,
           absolute = panel_absolute_risk(data, x.lim.absrisk, y.lim.absrisk,
                                          x.lab.absrisk, y.lab.absrisk, base.size),
           relative = panel_relative_risk(data, x.lim.RR, y.lim.RR, x.lab.RR, y.lab.RR,
                                          base.size),
           discrimination = panel_discrimination(data, risk.score.plot.kernel,
                                                 risk.score.plot.bandwidth, base.size),
           incidence = panel_incidence(data, base.size),
           eo = panel_expected_observed(data, base.size))
  })

  figure <- patchwork::wrap_plots(plots, ncol = min(2L, length(plots))) +
    patchwork::plot_annotation(
      title = model.name,
      subtitle = validation_subtitle(data, dataset),
      caption = validation_caption(data, panels),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = base.size + 3,
                                           colour = plot_colours[["text"]]),
        plot.subtitle = ggplot2::element_text(size = base.size,
                                              colour = plot_colours[["text_secondary"]],
                                              lineheight = 1.2,
                                              margin = ggplot2::margin(t = 3, b = 6)),
        plot.caption = ggplot2::element_text(size = base.size - 1.5, hjust = 0,
                                             colour = plot_colours[["text_secondary"]],
                                             lineheight = 1.2,
                                             margin = ggplot2::margin(t = 8)),
        plot.title.position = "plot",
        plot.caption.position = "plot"
      )
    ) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = base.size + 1,
                                                    colour = plot_colours[["text"]]))
  if (draw) {
    print(figure)
    return(invisible(figure))
  }
  figure
}

# ----- Data -------------------------------------------------------------------

#' Everything the panels need, from a ModelValidation() result
#' @noRd
validation_plot_data <- function(results, study.data) {
  study.data <- study.data %||% results$input.args$study.data
  outcome <- results$Subject_Specific_Observed_Outcome
  n <- length(outcome)
  sampling_prob <- NULL
  entry <- NULL
  if (!is.null(study.data)) {
    study.data <- as_plain_data_frame(study.data, "study.data")
    if (nrow(study.data) != n) {
      abort("`study.data` must be the data given to ModelValidation() (", n,
            " rows); it has ", nrow(study.data), ".")
    }
    sampling_prob <- study.data[["sampling.weights"]]
    entry <- study.data[["study.entry.age"]]
  }
  weights <- if (is.null(sampling_prob)) rep(1, n) else 1 / sampling_prob

  categories <- results$Category_Results
  column <- function(name) {
    if (is.null(categories[[name]])) return(NULL)
    as.numeric(as.character(categories[[name]]))
  }
  groups <- data.frame(
    label = as.character(categories$Categories),
    expected = column("Predicted_Absolute_Risk"),
    observed = column("Observed_Absolute_Risk"),
    lower = column("CI_Absolute_Risk_Lower"),
    upper = column("CI_Absolute_Risk_Upper"),
    expected_rr = column("Predicted_Relative_Risk"),
    observed_rr = column("Observed_Relative_Risk"),
    lower_rr = column("CI_Relative_Risk_Lower"),
    upper_rr = column("CI_Relative_Risk_Upper"),
    stringsAsFactors = FALSE
  )
  if (!is.null(column("Expected_to_Observed_Ratio"))) {
    groups$eo <- column("Expected_to_Observed_Ratio")
    groups$lower_eo <- column("CI_Expected_to_Observed_Ratio_Lower")
    groups$upper_eo <- column("CI_Expected_to_Observed_Ratio_Upper")
  }

  list(
    groups = groups,
    score = results$Subject_Specific_Risk_Score,
    outcome = outcome, weights = weights, sampled = !is.null(sampling_prob),
    entry = entry, followup = results$Adjusted_Followup,
    eo = results$Overall_Expected_to_Observed_Ratio,
    eo_ci = results$CI_Overall_Expected_to_Observed_Ratio,
    auc = results$AUC, auc_ci = results$CI_AUC,
    hl = results$Hosmer_Lemeshow_Results, rr = results$RR_test_result,
    hl_variance = results$input.args$hosmer.lemeshow.variance %||% "observed",
    interval = results$Risk_Prediction_Interval,
    population = results$Population_Incidence_Rate,
    study_incidence = results$Study_Incidence_Rate,
    dataset = results$input.args$dataset %||% "",
    model.name = results$input.args$model.name %||% ""
  )
}

check_axis_limits <- function(limits, name, positive = FALSE) {
  if (is.null(limits)) {
    return(invisible(TRUE))
  }
  if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits) ||
      any(!is.finite(limits)) || limits[1] >= limits[2] || (positive && limits[1] <= 0)) {
    abort("`", name, "` must be two increasing", if (positive) " positive" else "",
          " numbers.")
  }
  invisible(TRUE)
}

# ----- Text -------------------------------------------------------------------

format_number <- function(x, digits) formatC(x, format = "f", digits = digits)

format_count <- function(x) format(round(x), big.mark = ",", trim = TRUE)

# Estimate with its interval, e.g. "0.95 (95% CI 0.78-1.16)"
format_estimate <- function(estimate, ci, digits) {
  if (!is.finite(estimate)) {
    return("not estimable")
  }
  paste0(format_number(estimate, digits), " (95% CI ", format_number(ci[1], digits),
         "-", format_number(ci[2], digits), ")")
}

format_p_value <- function(p) {
  if (length(p) != 1L || is.na(p)) {
    return("not computed")
  }
  if (p < 0.001) {
    return("< 0.001")
  }
  paste("=", format_number(p, if (p < 0.01) 3 else 2))
}

# A chi-square test, e.g. "chi-square 16.6, df 10, p = 0.08"
format_test <- function(test, digits = 1) {
  statistic <- unname(test$statistic)
  if (length(statistic) != 1L || !is.finite(statistic)) {
    return("not computed")
  }
  paste0("chi-square ", format_number(statistic, digits), ", df ", test$parameter,
         ", p ", format_p_value(test$p.value))
}

validation_subtitle <- function(data, dataset) {
  cases <- sum(data$outcome)
  counts <- paste0(format_count(length(data$outcome)), " subjects, ",
                   format_count(cases), " cases")
  if (data$sampled) {
    counts <- paste0(counts, " (weighted: ", format_count(sum(data$weights)),
                     " subjects, ", format_count(sum(data$weights * data$outcome)),
                     " cases)")
  }
  details <- c(paste("Risk interval:", data$interval),
               sprintf("Mean follow-up %.1f years", weighted_average(data$followup,
                                                                     data$weights)))
  if (!is.null(data$entry)) {
    details <- c(details, sprintf("Age at entry %s-%s years", format(min(data$entry)),
                                  format(max(data$entry))))
  }
  lines <- c(if (nzchar(dataset)) dataset, counts, paste(details, collapse = "  |  "))
  paste(lines, collapse = "\n")
}

validation_caption <- function(data, panels) {
  n_groups <- nrow(data$groups)
  notes <- character(0)
  if (any(panels %in% c("absolute", "relative", "eo"))) {
    notes <- c(notes, paste0(
      "Points: ", n_groups, " groups of the risk score; bars: 95% confidence intervals",
      if (data$sampled) " (weighted by inverse inclusion probabilities)", "."))
  }
  if ("absolute" %in% panels) {
    notes <- c(notes, paste0("Hosmer-Lemeshow test with the variance of the ",
                             data$hl_variance, " risks."))
  }
  if ("incidence" %in% panels) {
    notes <- c(notes, "Incidence: cases per person-year at each year of attained age.")
  }
  paste(strwrap(paste(notes, collapse = " "), width = 110), collapse = "\n")
}

# ----- Look -------------------------------------------------------------------

# Validated categorical pair (blue, orange) for the two-series panels; the
# neutral grey marks reference lines only.
plot_colours <- c(primary = "#2a78d6", secondary = "#eb6834", reference = "#8c8b87",
                  text = "#1a1a19", text_secondary = "#52514e", grid = "#ebebe8")

plot_theme <- function(base.size) {
  ggplot2::theme_minimal(base_size = base.size) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(face = "bold", size = base.size,
                                         colour = plot_colours[["text"]],
                                         margin = ggplot2::margin(b = 6)),
      plot.title.position = "panel",
      plot.subtitle = ggplot2::element_text(size = base.size - 1, lineheight = 1.2,
                                            colour = plot_colours[["text_secondary"]],
                                            margin = ggplot2::margin(b = 6)),
      panel.grid.major = ggplot2::element_line(colour = plot_colours[["grid"]],
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = plot_colours[["text_secondary"]],
                                        linewidth = 0.3),
      axis.ticks = ggplot2::element_line(colour = plot_colours[["text_secondary"]],
                                         linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = plot_colours[["text_secondary"]],
                                        size = base.size - 1),
      axis.title = ggplot2::element_text(colour = plot_colours[["text"]],
                                         size = base.size - 0.5),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(colour = plot_colours[["text"]],
                                          size = base.size - 1),
      legend.key.size = ggplot2::unit(0.9, "lines"),
      legend.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 8, 4, 4)
    )
}

# Points with vertical interval bars
point_interval_layers <- function() {
  list(
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                            colour = plot_colours[["primary"]], linewidth = 0.55,
                            na.rm = TRUE),
    ggplot2::geom_point(shape = 21, size = 2.6, stroke = 0.6, colour = "white",
                        fill = plot_colours[["primary"]], na.rm = TRUE)
  )
}

identity_line <- function() {
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                       colour = plot_colours[["reference"]], linewidth = 0.45)
}

# Nice breaks for a log axis within the limits; thinning keeps 1 when shown
log_breaks <- function(limits, n = 6) {
  candidates <- c(outer(c(1, 2, 5), 10^(-6:6)))
  span <- diff(log10(limits))
  if (span < 1) candidates <- c(candidates, outer(c(1.5, 3, 4, 6, 8), 10^(-6:6)))
  if (span < 0.4) candidates <- c(candidates, outer(c(1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2.5),
                                                    10^(-6:6)))
  breaks <- sort(unique(signif(candidates, 6)))
  breaks <- breaks[breaks >= limits[1] & breaks <= limits[2]]
  while (length(breaks) > n) {
    anchor <- match(1, breaks, nomatch = 1L)
    breaks <- breaks[(seq_along(breaks) - anchor) %% 2L == 0L]
  }
  breaks
}

plain_numbers <- function(x) format(x, trim = TRUE, drop0trailing = TRUE, scientific = FALSE)

# ----- Panels -----------------------------------------------------------------

panel_absolute_risk <- function(data, xlim, ylim, xlab, ylab, base.size) {
  groups <- data$groups
  frame <- data.frame(x = 100 * groups$expected, y = 100 * groups$observed,
                      lower = 100 * pmax(groups$lower, 0), upper = 100 * groups$upper)
  default <- c(0, 1.06 * max(c(frame$x, frame$y, frame$upper), na.rm = TRUE))
  xlim <- xlim %||% default
  ylim <- ylim %||% default
  statistics <- paste0("E/O ", format_estimate(data$eo, data$eo_ci, 2),
                       "\nHosmer-Lemeshow ", format_test(data$hl))
  ggplot2::ggplot(frame, ggplot2::aes(x = .data$x, y = .data$y)) +
    identity_line() +
    point_interval_layers() +
    ggplot2::scale_x_continuous(breaks = pretty(xlim, 5), labels = plain_numbers) +
    ggplot2::scale_y_continuous(breaks = pretty(ylim, 5), labels = plain_numbers) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(title = "Absolute risk calibration", subtitle = statistics,
                  x = xlab, y = ylab) +
    plot_theme(base.size)
}

panel_relative_risk <- function(data, xlim, ylim, xlab, ylab, base.size) {
  groups <- data$groups
  frame <- data.frame(x = groups$expected_rr, y = groups$observed_rr,
                      lower = groups$lower_rr, upper = groups$upper_rr)
  frame <- frame[is.finite(frame$y) & frame$y > 0, , drop = FALSE]
  values <- c(frame$x, frame$y, frame$lower, frame$upper)
  values <- values[is.finite(values) & values > 0]
  span <- range(log(values))
  default <- exp(span + c(-1, 1) * 0.06 * max(diff(span), log(1.5)))
  xlim <- xlim %||% default
  ylim <- ylim %||% default
  ggplot2::ggplot(frame, ggplot2::aes(x = .data$x, y = .data$y)) +
    identity_line() +
    point_interval_layers() +
    ggplot2::scale_x_log10(breaks = log_breaks(xlim), labels = plain_numbers) +
    ggplot2::scale_y_log10(breaks = log_breaks(ylim), labels = plain_numbers) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(title = "Relative risk calibration",
                  subtitle = paste("Relative-risk test", format_test(data$rr)),
                  x = paste(xlab, "(log scale)"), y = paste(ylab, "(log scale)")) +
    plot_theme(base.size)
}

#' Bandwidth of a weighted kernel density
#'
#' Rules "nrd0" and "nrd" use the weighted standard deviation and
#' interquartile range with the Kish effective sample size; with equal
#' weights they equal stats::bw.nrd0() and stats::bw.nrd(). Other rules
#' ignore the weights.
#' @noRd
density_bandwidth <- function(x, w, rule) {
  if (is.numeric(rule)) {
    return(rule)
  }
  rule <- tolower(rule)
  if (!rule %in% c("nrd0", "nrd")) {
    selector <- switch(rule, ucv = stats::bw.ucv, bcv = stats::bw.bcv,
                       sj = , `sj-ste` = function(x) stats::bw.SJ(x, method = "ste"),
                       `sj-dpi` = function(x) stats::bw.SJ(x, method = "dpi"))
    return(selector(x))
  }
  if (length(unique(w)) == 1L) {
    return(if (rule == "nrd0") stats::bw.nrd0(x) else stats::bw.nrd(x))
  }
  share <- w / sum(w)
  sd_w <- sqrt(sum(share * (x - sum(share * x))^2))
  quartiles <- weighted_quantile_interpolated(x, w, c(0.25, 0.75))
  spread <- min(sd_w, diff(quartiles) / 1.34)
  if (!(spread > 0)) spread <- if (sd_w > 0) sd_w else max(abs(x[1]), 1) * 0.1
  effective_n <- sum(w)^2 / sum(w^2)
  (if (rule == "nrd0") 0.9 else 1.06) * spread * effective_n^(-0.2)
}

panel_discrimination <- function(data, kernel, bandwidth, base.size) {
  valid_rules <- c("nrd0", "nrd", "ucv", "bcv", "sj", "sj-ste", "sj-dpi")
  if (!((is.numeric(bandwidth) && length(bandwidth) == 1L && is.finite(bandwidth) &&
         bandwidth > 0) ||
        (is.character(bandwidth) && length(bandwidth) == 1L &&
         tolower(bandwidth) %in% valid_rules))) {
    abort("`risk.score.plot.bandwidth` must be a positive number or one of ",
          format_names(c("nrd0", "nrd", "ucv", "bcv", "SJ")), ".")
  }
  groups <- list(`Non-cases` = data$outcome == 0, Cases = data$outcome == 1)
  bandwidths <- vapply(groups, function(g) {
    density_bandwidth(data$score[g], data$weights[g], bandwidth)
  }, numeric(1))
  from <- min(data$score) - 3 * max(bandwidths)
  to <- max(data$score) + 3 * max(bandwidths)
  curves <- do.call(rbind, lapply(names(groups), function(name) {
    g <- groups[[name]]
    w <- data$weights[g]
    estimate <- stats::density(data$score[g], weights = w / sum(w),
                               bw = bandwidths[[name]], kernel = kernel, n = 512,
                               from = from, to = to)
    data.frame(x = estimate$x, y = estimate$y, group = name)
  }))
  curves$group <- factor(curves$group, levels = names(groups))
  colours <- c(`Non-cases` = plot_colours[["primary"]], Cases = plot_colours[["secondary"]])
  ggplot2::ggplot(curves, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group,
                                       fill = .data$group)) +
    ggplot2::geom_area(position = "identity", alpha = 0.18, colour = NA) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(c(0, 0.2))) +
    ggplot2::labs(title = "Discrimination",
                  subtitle = paste("AUC", format_estimate(data$auc, data$auc_ci, 3)),
                  x = "Risk score (linear predictor)",
                  y = if (data$sampled) "Density (weighted)" else "Density") +
    plot_theme(base.size) +
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(0.98, 0.98),
                   legend.justification.inside = c(1, 1),
                   panel.grid.major.y = ggplot2::element_blank())
}

panel_incidence <- function(data, base.size) {
  incidence <- data$study_incidence
  has_ci <- ncol(incidence) >= 4L
  # Rates of the year of age [a, a + 1) are drawn at its middle, a + 0.5.
  study <- data.frame(age = incidence[[1]] + 0.5, rate = 1e5 * incidence[[2]],
                      lower = if (has_ci) 1e5 * incidence[[3]] else NA_real_,
                      upper = if (has_ci) 1e5 * incidence[[4]] else NA_real_)
  study <- study[is.finite(study$rate) & study$rate > 0, , drop = FALSE]
  population <- NULL
  if (!is.null(data$population)) {
    rates <- normalize_rate_table(data$population, "Population_Incidence_Rate")
    population <- data.frame(age = rates$age + 0.5, rate = 1e5 * rates$rate)
    if (nrow(study) > 0L) {
      population <- population[population$age >= min(study$age) &
                                 population$age <= max(study$age), , drop = FALSE]
    }
    population <- population[population$rate > 0, , drop = FALSE]
  }

  layers <- list()
  sources <- character(0)
  if (!is.null(population) && nrow(population) > 0L) {
    population$source <- "Population"
    layers <- c(layers, list(ggplot2::geom_line(data = population, linewidth = 0.7)))
    sources <- c(sources, "Population")
  }
  if (nrow(study) > 0L) {
    study$source <- "Study"
    layers <- c(layers, list(
      ggplot2::geom_linerange(data = study, ggplot2::aes(ymin = .data$lower,
                                                         ymax = .data$upper),
                              linewidth = 0.4, alpha = 0.55, na.rm = TRUE,
                              show.legend = FALSE),
      ggplot2::geom_point(data = study, size = 1.6)
    ))
    sources <- c(sources, "Study")
  }
  # The axis covers the rates; intervals reaching far beyond them are cut at
  # the edge of the panel.
  rates <- c(study$rate, population$rate)
  limits <- NULL
  if (length(rates) > 0L) {
    upper <- c(rates, pmin(study$upper, 2.5 * max(rates)))
    limits <- c(min(rates) / 2.5, 1.1 * max(upper, na.rm = TRUE))
  }
  colours <- c(Population = plot_colours[["secondary"]], Study = plot_colours[["primary"]])
  key <- list(shape = unname(c(Population = NA, Study = 16)[sources]),
              linetype = unname(c(Population = 1, Study = 0)[sources]))
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$age, y = .data$rate,
                                         colour = .data$source)) +
    layers +
    ggplot2::scale_colour_manual(values = colours, breaks = sources,
                                 guide = ggplot2::guide_legend(override.aes = key)) +
    ggplot2::scale_y_log10(breaks = log_breaks(limits %||% c(1, 1e5), n = 7),
                           labels = plain_numbers) +
    ggplot2::coord_cartesian(ylim = limits) +
    ggplot2::labs(title = "Incidence by attained age",
                  subtitle = "Cases per 100,000 person-years",
                  x = "Attained age (years)", y = "Incidence (log scale)") +
    plot_theme(base.size) +
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(0.02, 0.98),
                   legend.justification.inside = c(0, 1))
}

panel_expected_observed <- function(data, base.size) {
  groups <- data$groups
  if (is.null(groups$eo)) {
    abort("These validation results have no E/O by risk-score group; rerun ",
          "ModelValidation() with this version of iCARE.")
  }
  frame <- data.frame(x = seq_len(nrow(groups)), y = groups$eo,
                      lower = groups$lower_eo, upper = groups$upper_eo)
  frame <- frame[is.finite(frame$y), , drop = FALSE]
  values <- c(frame$y, frame$lower, frame$upper, 1)
  values <- values[is.finite(values) & values > 0]
  half <- max(abs(log(range(values))), log(1.25)) * 1.08
  limits <- exp(c(-half, half))
  ggplot2::ggplot(frame, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = plot_colours[["reference"]], linewidth = 0.45) +
    point_interval_layers() +
    ggplot2::scale_x_continuous(breaks = seq_len(nrow(groups))) +
    ggplot2::scale_y_log10(breaks = log_breaks(limits, n = 7), labels = plain_numbers) +
    ggplot2::coord_cartesian(xlim = c(0.4, nrow(groups) + 0.6), ylim = limits,
                             expand = FALSE) +
    ggplot2::labs(title = "Expected/observed by risk group",
                  subtitle = paste("Overall E/O", format_estimate(data$eo, data$eo_ci, 2)),
                  x = "Risk-score group (lowest to highest)",
                  y = "Expected / observed (log scale)") +
    plot_theme(base.size) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}
