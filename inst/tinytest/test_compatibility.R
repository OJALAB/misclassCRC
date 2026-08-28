library(misclassCRC)

data <- data.frame(
  source_1 = c(1L, 0L, 0L, 1L),
  source_2 = c(0L, 1L, 1L, 0L),
  source_3 = c(0L, 0L, 1L, 1L),
  group = c("B", "A", "B", "A"),
  region = c("north", "south", "north", "south"),
  outcome = c(2, 3, 4, 5)
)
group_matrix <- matrix(
  c(0.8, 0.2, 0.1, 0.9), nrow = 2L, byrow = TRUE,
  dimnames = list(c("A", "B"), c("A", "B"))
)

expect_error(
  misclassCRC:::prepare_crc(
    data = data,
    captures = ~ source_1 + source_2 + source_3,
    capture_formula = ~ source_1 + source_2 + source_3,
    outcome = "outcome",
    outcome_formula = ~ region,
    outcome_dist = "ztnegbin",
    misclass = NULL,
    latent_classes = 1L,
    control = NULL,
    verbose = FALSE
  ),
  pattern = "not included in `capture_formula`"
)

expect_error(
  misclassCRC:::prepare_crc(
    data = data,
    captures = ~ source_1 + source_2 + source_3,
    capture_formula = ~ source_1 + source_2 + source_3,
    outcome = NULL,
    outcome_formula = NULL,
    outcome_dist = "ztnegbin",
    misclass = list(group = list(matrix = group_matrix)),
    latent_classes = 1L,
    control = NULL,
    verbose = FALSE
  ),
  pattern = "not used in either model"
)
