library(misclassCRC)

data <- data.frame(
  source_1 = c(1L, 0L, 0L, 1L),
  source_2 = c(0L, 1L, 1L, 0L),
  source_3 = c(0L, 0L, 1L, 1L)
)
prepare_call <- function(control = NULL, verbose = FALSE) {
  misclassCRC:::prepare_crc(
    data, ~ source_1 + source_2 + source_3,
    ~ source_1 + source_2 + source_3,
    outcome = NULL, outcome_formula = NULL, outcome_dist = "ztnegbin",
    misclass = NULL, latent_classes = 1L,
    control = control, verbose = verbose
  )
}

preparation <- prepare_call()
expect_inherits(preparation$control, "crc_control")
expect_equal(preparation$control$init_alpha, 20)
expect_equal(preparation$control$em, list(max_iter = 1000L, tolerance = 1e-6))
expect_false(preparation$verbose)

custom <- prepare_call(list(init_alpha = 12, em = list(max_iter = 25L),
  outcome = list(relative_tolerance = 1e-7)), verbose = TRUE)
expect_equal(custom$control$init_alpha, 12)
expect_equal(custom$control$em, list(max_iter = 25L, tolerance = 1e-6))
expect_equal(custom$control$outcome$relative_tolerance, 1e-7)
expect_true(custom$verbose)

expect_error(prepare_call(list(unknown = 1)), pattern = "Unknown setting")
expect_error(prepare_call(list(em = list(unknown = 1))), pattern = "control\\$em")
expect_error(prepare_call(list(capture = 1)), pattern = "must be a list")
expect_error(prepare_call(list(em = list(max_iter = 1.5))), pattern = "positive integer")
expect_error(prepare_call(verbose = 1), pattern = "logical value")
