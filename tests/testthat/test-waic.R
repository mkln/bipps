logSumExp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x-m)))
}

r_waic <- function(log_likelihoods) {
  N <- nrow(log_likelihoods)
  S <- ncol(log_likelihoods)

  # Compute lppd (log point-wise predictive density)
  lppd <- apply(log_likelihoods, 1, function(x) logSumExp(x) - log(S))

  # Compute p_waic (penalty term)
  penalty <- apply(log_likelihoods, 1, var)

  # Calculate WAIC
  waic <- -2 * sum(lppd) + 2 * sum(penalty)
  return(waic)
}

test_that("Calculate WAIC", {
  N <- 100
  S <- 100
  log_likelihoods <- matrix(dnorm(rnorm(N * S, mean = 0, sd = 1),log=TRUE), nrow = N, ncol = S)

  streaming_waic <- waic_test_runner(log_likelihoods)

  # out <- loo::waic(t(log_likelihoods))
  # out$estimates
  #
  # LaplacesDemon::WAIC(log_likelihoods)

  simple_waic <- r_waic(log_likelihoods)

  expect_equal(streaming_waic, simple_waic)
})
