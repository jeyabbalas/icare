# Summarize benchmark results as a Markdown table
# ===============================================
#
# Usage (from the package root, after run-benchmarks.R has written the
# results files):
#
#   Rscript tools/benchmarks/summarize-benchmarks.R [reference-label]
#
# Reads tools/benchmarks/results-*.csv and prints one row per workload with
# the median elapsed seconds of each labelled version, and the speed-up of
# the first version over `reference-label` (default "official-c").

args <- commandArgs(trailingOnly = TRUE)
reference_label <- if (length(args) >= 1L) args[1] else "official-c"

files <- list.files(file.path("tools", "benchmarks"), "^results-.*[.]csv$",
                    full.names = TRUE)
if (length(files) == 0L) stop("No results files; run run-benchmarks.R first.")
results <- do.call(rbind, lapply(files, utils::read.csv))

labels <- unique(results$label)
# This version first, then the others in alphabetical order
labels <- c(intersect("refactored", labels), sort(setdiff(labels, "refactored")))
workloads <- unique(results$workload)

seconds <- vapply(labels, function(l) {
  r <- results[results$label == l, ]
  r$seconds[match(workloads, r$workload)]
}, numeric(length(workloads)))
seconds <- matrix(seconds, nrow = length(workloads), dimnames = list(workloads, labels))

format_seconds <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 1, sprintf("%.3f", x), sprintf("%.2f", x)))
}
version_of <- function(l) {
  r <- results[results$label == l, ][1, ]
  paste0(r$version, if (!is.na(r$use_c_code)) paste0(", use.c.code = ", r$use_c_code))
}

header <- c("Workload", paste0(labels, " (", vapply(labels, version_of, ""), ")"))
rows <- lapply(workloads, function(w) {
  description <- results$description[match(w, results$workload)]
  c(description, format_seconds(seconds[w, ]))
})
if (reference_label %in% labels && labels[1] != reference_label) {
  header <- c(header, paste0("Speed-up vs ", reference_label))
  rows <- lapply(seq_along(workloads), function(i) {
    ratio <- seconds[i, reference_label] / seconds[i, 1]
    c(rows[[i]], sprintf("%.1f×", ratio))
  })
}

cat("| ", paste(header, collapse = " | "), " |\n", sep = "")
cat("|", paste(c(":--", rep("--:", length(header) - 1L)), collapse = "|"), "|\n", sep = "")
for (row in rows) cat("| ", paste(row, collapse = " | "), " |\n", sep = "")
cat("\nMedian elapsed seconds over three runs (one run when longer than 20 s).\n")
