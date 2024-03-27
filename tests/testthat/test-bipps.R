library(magrittr)
library(dplyr)
library(tidyr)
test_that("bipps reproducibility across runs", {
  set.seed(2024)
  test_data <- get_test_data()

  set.seed(2024)
  reference_data <- get_test_data()

  expect_equal(test_data$beta_test,reference_data$beta_test)
  expect_equal(test_data$lambda_test,reference_data$lambda_test)
  expect_equal(test_data$theta_test,reference_data$theta_test)
  expect_equal(test_data$w_test,reference_data$w_test)
})
