#include "multi_bipps.h"

using namespace std;

void MultiBipps::sample_hmc_BetaLambdaTau(bool sample, bool sample_beta, bool sample_lambda, bool sample_tau){
  if(verbose & debug){
    Rcpp::Rcout << "[sample_hmc_BetaLambdaTau] starting\n";
  }
  start = std::chrono::steady_clock::now();

  double mat_sums = 0;
  std::vector<arma::vec> sampleds;
  std::vector<int> indices;

  Rcpp::RNGScope scope;
  arma::mat rnorm_precalc = mrstdnorm(q, k+p);
  arma::vec lambda_runif = vrunif(q);
  arma::vec lambda_runif2 = vrunif(q);
  
  arma::vec tau_rnorm_precalc = mrstdnorm(q, 1);
  arma::vec tau_runif_precalc = vrunif(q);
  
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(auto j=0; j<q; j++){
    
    ///
    /// ** Beta & Lambda update **
    ///

    
    // build W
    // filter: choose value of spatial processes at locations of Yj that are available


    arma::mat XW_joined;
    arma::vec offsets_joined;
    arma::uvec subcols = arma::find(Lambda_mask.row(j) == 1);

    for(Bipps& bipps : multi_bipps) {
      arma::vec offsets_obs = bipps.offsets(bipps.ix_by_q_a(j), oneuv * j);
      arma::vec y_obs = bipps.y(bipps.ix_by_q_a(j), oneuv * j);

      arma::mat WWj = bipps.w.submat(bipps.ix_by_q_a(j), subcols); // acts as X
      if(!sample) { apply2sd(WWj); }

      arma::mat XW = arma::join_horiz(bipps.X.rows(bipps.ix_by_q_a(j)), WWj);

      XW_joined = arma::join_vert(XW_joined, XW);
      offsets_joined = arma::join_vert(offsets_joined, offsets_obs);
    }
    arma::mat BL_Vi = arma::eye( XW_joined.n_cols, XW_joined.n_cols );

    BL_Vi.submat(0, 0, p-1, p-1) = Vi; // prior precision for beta
    arma::vec BL_Vim = arma::zeros(XW_joined.n_cols);
    
    NodeDataB& lambda_block = lambda_node.at(j);
    lambda_block.update_mv(offsets_joined, 1.0/tausq_inv(j), BL_Vim, BL_Vi);
    lambda_block.X = XW_joined;

    arma::vec curLrow = arma::join_vert(
      multi_Beta.col(j),
      arma::trans(multi_Lambda.submat(oneuv*j, subcols)));
    
    arma::mat rnorm_row = arma::trans(rnorm_precalc.row(j).head(curLrow.n_elem));
    
    AdaptE& lambda_adapt = lambda_hmc_adapt.at(j);

    // nongaussian
    lambda_adapt.step();

    // may want to simplfy sampler selection
    if((lambda_hmc_started(j) == 0) && (lambda_adapt.i == 10)){
      // wait a few iterations before starting adaptation
      // Rcpp::Rcout << "reasonable stepsize " << endl;
      
      double lambda_eps = find_reasonable_stepsize(curLrow, lambda_block, rnorm_row);
      
      int n_params = curLrow.n_elem;
      AdaptE new_adapting_scheme;
      new_adapting_scheme.init(lambda_eps, n_params, which_hmc, 1e4);
      lambda_adapt = new_adapting_scheme;
      lambda_hmc_started(j) = 1;
      // Rcpp::Rcout << "done initiating adapting scheme" << endl;
    }

    arma::vec sampled;
    if(which_hmc == 0){
      // some form of manifold mala
      sampled = simpa_cpp(curLrow, lambda_block, lambda_adapt, 
                              rnorm_row, lambda_runif(j), lambda_runif2(j), 
                              debug);
    }
    if(which_hmc == 1){
      // mala
      sampled = mala_cpp(curLrow, lambda_block, lambda_adapt, 
                         rnorm_row, lambda_runif(j), debug);
    }
    if(which_hmc == 2){
      Rcpp::stop("HMC algorithm 2 not implemented");
    }
    
    if(which_hmc == 3){
      // some form of manifold mala
      sampled = smmala_cpp(curLrow, lambda_block, lambda_adapt, 
                              rnorm_row, lambda_runif(j), debug);
    }
    if(which_hmc == 6){
      sampled = hmc_cpp(curLrow, lambda_block, lambda_adapt, 
                        rnorm_row, lambda_runif(j), 0.1, debug);
    }
    if(which_hmc == 7){
      Rcpp::stop("HMC algorithm 7 not implemented");
    }
    if(sample_lambda){
      multi_Lambda.submat(oneuv*j, subcols) = arma::trans(sampled.tail(subcols.n_elem));
    }
    if(sample_beta){
      multi_Beta.col(j) = sampled.head(p);
    } 
    for(Bipps& bipps: multi_bipps) {
      if(sample_beta){
        bipps.Bcoeff.col(j) = sampled.head(p);
      } 
      if(sample_lambda){
        bipps.Lambda.submat(oneuv*j, subcols) = arma::trans(sampled.tail(subcols.n_elem));
      }
      bipps.XB.col(j) = bipps.X * bipps.Bcoeff.col(j); 
      bipps.LambdaHw.col(j) = bipps.w * arma::trans(bipps.Lambda.row(j));
    }
}

  if(sample_lambda) {
    // ensure positive diag, must be done outside of the OMP loop
    multi_Lambda = multi_Lambda * arma::diagmat(arma::sign(multi_Lambda.diag()));
  }
  
  if(verbose & debug){
    Rcpp::Rcout << "[sample_hmc_BetaLambdaTau] XW_joined samples\n";
  }
  

  // refreshing density happens in the 'logpost_refresh_after_gibbs' function
  if(verbose & debug){
    Rcpp::Rcout << "[sample_hmc_Lambda] done\n";
  }
  
}
