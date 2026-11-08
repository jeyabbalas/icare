test_that("two-column tables are sorted yearly rates", {
  table <- normalize_rate_table(data.frame(age = c(52, 50, 51), rate = c(3, 1, 2) / 1000),
                                "rates", consecutive = TRUE)
  expect_equal(table, data.frame(age = c(50, 51, 52), rate = c(1, 2, 3) / 1000))
  # matrices and tibble-like data frames are accepted
  expect_equal(normalize_rate_table(cbind(50:52, c(1, 2, 3) / 1000), "rates"), table)
  tbl <- structure(data.frame(age = 50:52, rate = c(1, 2, 3) / 1000),
                   class = c("tbl_df", "tbl", "data.frame"))
  expect_equal(normalize_rate_table(tbl, "rates"), table)
})

test_that("three-column bands give the band's annual rate at every age", {
  inclusive <- data.frame(start = c(50, 55), end = c(54, 59), rate = c(0.002, 0.003))
  shared <- data.frame(start = c(50, 55), end = c(55, 60), rate = c(0.002, 0.003))
  expected <- data.frame(age = 50:59, rate = rep(c(0.002, 0.003), each = 5))
  expect_equal(normalize_rate_table(inclusive, "rates"), expected)
  expect_equal(normalize_rate_table(shared, "rates"), expected)
  # single-year bands and a single band are inclusive
  expect_equal(normalize_rate_table(data.frame(50:52, 50:52, c(1, 2, 3) / 1000), "rates"),
               data.frame(age = 50:52, rate = c(1, 2, 3) / 1000))
  expect_equal(normalize_rate_table(data.frame(50, 52, 0.001), "rates"),
               data.frame(age = 50:52, rate = rep(0.001, 3)))
})

test_that("missing rates mark ages that the table does not cover", {
  table <- normalize_rate_table(data.frame(age = 50:53, rate = c(0.001, 0.002, NA, NA)),
                                "rates")
  expect_equal(table$age, c(50, 51))
  expect_error(normalize_rate_table(data.frame(age = 50:51, rate = c(NA, NA)), "rates"),
               "no non-missing rates")
})

test_that("invalid rate tables are rejected with clear errors", {
  expect_error(normalize_rate_table(data.frame(a = 1:3), "x"), "two columns")
  expect_error(normalize_rate_table(data.frame(a = 1:2, b = 1:2, c = 1:2, d = 1:2), "x"),
               "it has 4")
  expect_error(normalize_rate_table(data.frame(age = c(50, NA), rate = 0.1), "x"),
               "ages .* cannot be missing")
  expect_error(normalize_rate_table(data.frame(age = 50:51, rate = c(0.1, Inf)), "x"),
               "cannot be infinite")
  expect_error(normalize_rate_table(data.frame(age = 50:51, rate = c(0.1, 2)), "x"),
               "between 0 and 1")
  expect_error(normalize_rate_table(data.frame(age = c(50, 50.5), rate = 0.1), "x"),
               "whole numbers")
  expect_error(normalize_rate_table(data.frame(age = c(50, 50), rate = 0.1), "x"),
               "more than once")
  expect_error(normalize_rate_table(data.frame(age = c(50, 52), rate = 0.1), "x",
                                    consecutive = TRUE), "consecutive ages.*51")
  expect_error(normalize_rate_table(data.frame(c(50, 53), c(54, 58), c(0.1, 0.1)), "x"),
               "contiguous")
  expect_error(normalize_rate_table(data.frame(age = c("a", "b"), rate = 0.1), "x"),
               "numeric")
})

test_that("competing rates default to zero over the incidence ages", {
  tables <- resolve_rate_tables(data.frame(age = 50:52, rate = 0.001), NULL)
  expect_equal(tables$competing, data.frame(age = c(50, 51, 52), rate = 0))
  expect_error(resolve_rate_tables(NULL, NULL), "must be supplied")
})

test_that("rate coverage is checked for the ages a window uses", {
  table <- data.frame(age = 50:59, rate = 0.001)
  expect_true(check_rate_coverage(table, 50:59, "rates"))
  expect_error(check_rate_coverage(table, 58:61, "rates"), "60-61")
  expect_equal(rates_at(table, c(49, 50)), c(NA, 0.001))
  expect_equal(rates_at(table, c(49, 50), fill = 0), c(0, 0.001))
})
