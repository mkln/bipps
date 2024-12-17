#include <RcppArmadillo.h>
#include "utils_waic.h"

// [[Rcpp::export]]
double waic_test_runner(const arma::mat& log_likelihood) {
    int N = log_likelihood.n_rows;
		StreamingWAIC waic_calculator(N);
		int S = log_likelihood.n_cols;
		for (int i = 0; i < S; i++) {
			waic_calculator.update(log_likelihood.col(i));
		}
		return waic_calculator.calculate_waic();
}
