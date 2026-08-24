if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package("misclassCRC", ncpu = 1)
}