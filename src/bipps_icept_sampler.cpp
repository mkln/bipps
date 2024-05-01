#include "bipps.h"

using namespace std;


void Bipps::sample_hmc_icept(){
  if(verbose & debug){
    Rcpp::Rcout << "[sample_hmc_icept] starting\n";
  }
  start = std::chrono::steady_clock::now();
  
  Rcpp::RNGScope scope;
  arma::mat rnorm_precalc = mrstdnorm(q, 1);
  arma::vec icept_runif = vrunif(q);
  arma::vec icept_runif2 = vrunif(q);
  
  for(unsigned int j=0; j<q; j++){
    
    ///
    /// ** Intercept update **
    ///
    
    arma::vec offsets_obs = offsets(ix_by_q_a(j), oneuv * j) + 
      XB(ix_by_q_a(j), oneuv * j) + 
      LambdaHw(ix_by_q_a(j), oneuv * j);
    
    arma::mat XW = arma::ones(ix_by_q_a(j).n_elem, 1);
    
    arma::mat BL_Vi = 1.0/1000.0 * arma::eye( XW.n_cols, XW.n_cols );
    arma::vec BL_Vim = arma::zeros(XW.n_cols);
    
    lambda_node.at(j).update_mv(offsets_obs, 1.0/tausq_inv(j), BL_Vim, BL_Vi);
    lambda_node.at(j).X = XW;
    
    arma::vec curLrow = arma::ones(1) * icept(j);
    arma::mat rnorm_row = arma::ones(1,1) * rnorm_precalc(j, 0); //arma::trans(rnorm_precalc.row(j).head(curLrow.n_elem));
    arma::vec sampled;
    
    // nongaussian
    icept_hmc_adapt.at(j).step();
    if((icept_hmc_started(j) == 0) && (icept_hmc_adapt.at(j).i == 10)){
      // wait a few iterations before starting adaptation
      //Rcpp::Rcout << "reasonable stepsize " << endl;
      
      double icept_eps = find_reasonable_stepsize(curLrow, icept_node.at(j), rnorm_row);
      
      int n_params = curLrow.n_elem;
      AdaptE new_adapting_scheme;
      new_adapting_scheme.init(icept_eps, n_params, which_hmc, 1e4);
      icept_hmc_adapt.at(j) = new_adapting_scheme;
      icept_hmc_started(j) = 1;
      //Rcpp::Rcout << "done initiating adapting scheme" << endl;
    }
    
    if(which_hmc == 0){
      // some form of manifold mala
      sampled = simpa_cpp(curLrow, icept_node.at(j), icept_hmc_adapt.at(j), 
                          rnorm_row, icept_runif(j), icept_runif2(j), 
                          debug);
    }
    if(which_hmc == 1){
      // mala
      sampled = mala_cpp(curLrow, icept_node.at(j), icept_hmc_adapt.at(j), 
                         rnorm_row, icept_runif(j), debug);
    }
    if(which_hmc == 2){
      Rcpp::stop("HMC algorithm 2 not implemented");
    }
    
    if(which_hmc == 3){
      // some form of manifold mala
      sampled = smmala_cpp(curLrow, icept_node.at(j), icept_hmc_adapt.at(j), 
                           rnorm_row, icept_runif(j), debug);
    }
    if(which_hmc == 6){
      sampled = hmc_cpp(curLrow, icept_node.at(j), icept_hmc_adapt.at(j), 
                        rnorm_row, icept_runif(j), 0.1, debug);
    }
    if(which_hmc == 7){
      Rcpp::stop("HMC algorithm 7 not implemented");
    }
    
    icept(j) = sampled(0);
  }
  
  // refreshing density happens in the 'logpost_refresh_after_gibbs' function
  if(verbose & debug){
    Rcpp::Rcout << "[sample_hmc_icept] done\n";
  }
  
}

