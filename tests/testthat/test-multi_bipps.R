library(magrittr)
library(dplyr)
library(tidyr)

test_that("Test of reproducibility across runs", {

  set.seed(2024)
  test_data <- get_test_data_multi()

  set.seed(2024)
  reference_data <- get_test_data_multi()
  # reference_data <- readRDS("reference_data_multi.rds")

  expect_equal(test_data$beta_test,reference_data$beta_test)
  expect_equal(test_data$lambda_test,reference_data$lambda_test)
  expect_equal(test_data$theta_test,reference_data$theta_test)
  expect_equal(test_data$w_test,reference_data$w_test)
})

test_that("Test of reproducibility betweem different numbers of threads of runs", {

  # data1 <- readRDS("reference_data_multi.rds")
  set.seed(2024)
  data1 <- get_test_data_multi(n_jobs = 1)

  set.seed(2024)
  data2 <- get_test_data_multi(n_jobs = 2)

  set.seed(2024)
  data3 <- get_test_data_multi(n_jobs = 3)

  expect_equal(data1$beta_test,data2$beta_test)
  expect_equal(data1$lambda_test,data2$lambda_test)
  expect_equal(data1$theta_test,data2$theta_test)
  expect_equal(data1$w_test,data2$w_test)

  expect_equal(data2$beta_test,data3$beta_test)
  expect_equal(data2$lambda_test,data3$lambda_test)
  expect_equal(data2$theta_test,data3$theta_test)
  expect_equal(data2$w_test,data3$w_test)
})
