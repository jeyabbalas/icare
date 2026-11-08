# Small internal helpers shared by all modules: condition signalling, input
# coercion and formatting. Nothing here is exported.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Signal an error / warning without the (uninformative) internal call.
abort <- function(...) stop(paste0(...), call. = FALSE)
warn <- function(...) warning(paste0(...), call. = FALSE)

# TRUE when `x` is NULL or a single NA (the package's "not supplied" values).
is_absent <- function(x) is.null(x) || (length(x) == 1L && is.atomic(x) && is.na(x))

# A single TRUE/FALSE flag; 0/1 are accepted for backward compatibility.
as_flag <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !(is.logical(x) || (is.numeric(x) && x %in% c(0, 1)))) {
    abort("`", name, "` must be TRUE or FALSE.")
  }
  as.logical(x)
}

# A single whole number within [lower, upper].
as_count <- function(x, name, lower = 1, upper = Inf) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || x %% 1 != 0 ||
      x < lower || x > upper) {
    range_text <- if (is.finite(upper)) {
      paste0("between ", lower, " and ", upper)
    } else {
      paste0("of at least ", lower)
    }
    abort("`", name, "` must be a single whole number ", range_text, ".")
  }
  as.integer(x)
}

# A plain data frame from a data frame, tibble or matrix. Row names are kept so
# that results can be matched back to the caller's rows.
as_plain_data_frame <- function(x, name) {
  if (is.matrix(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) {
    abort("`", name, "` must be a data frame (a tibble is fine) or a matrix.")
  }
  # Drops tibble / data.table classes, whose `[` methods behave differently.
  as.data.frame(x, stringsAsFactors = FALSE)
}

# "50-54, 60, 62-63": compact description of a set of integer ages.
format_age_ranges <- function(ages, max_ranges = 6L) {
  ages <- sort(unique(ages))
  if (length(ages) == 0L) {
    return("")
  }
  breaks <- c(0L, which(diff(ages) != 1), length(ages))
  ranges <- vapply(seq_len(length(breaks) - 1L), function(k) {
    from <- ages[breaks[k] + 1L]
    to <- ages[breaks[k + 1L]]
    if (from == to) format(from) else paste0(format(from), "-", format(to))
  }, character(1))
  if (length(ranges) > max_ranges) {
    ranges <- c(ranges[seq_len(max_ranges)], "...")
  }
  paste(ranges, collapse = ", ")
}

# Collapse a character vector for messages: "a, b and c".
format_names <- function(x, max_names = 8L) {
  x <- as.character(x)
  if (length(x) > max_names) {
    x <- c(x[seq_len(max_names)], paste0("... (", length(x) - max_names, " more)"))
  }
  if (length(x) <= 1L) {
    return(paste(x, collapse = ""))
  }
  paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
}

# Weighted mean that is exact for (and only for) non-negative weights with a
# positive sum; callers guarantee that.
weighted_average <- function(x, w) sum(x * w) / sum(w)

# Handles the arguments removed from the public functions. `use.c.code` is
# accepted with a deprecation warning and ignored (the compiled code was
# removed, the R implementation is used); anything else is an error, exactly
# as R itself would report an unused argument.
check_removed_arguments <- function(dots, fn) {
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }
  dot_names <- names(dots) %||% rep("", length(dots))
  if (any(dot_names == "")) {
    abort("Unused unnamed argument(s) in ", fn, "().")
  }
  unknown <- setdiff(dot_names, "use.c.code")
  if (length(unknown) > 0L) {
    abort("Unused argument(s) in ", fn, "(): ", format_names(unknown), ".")
  }
  warning("`use.c.code` is deprecated and ignored: this version of iCARE has ",
          "no compiled code and always uses its vectorized R implementation.",
          call. = FALSE)
  invisible(NULL)
}
