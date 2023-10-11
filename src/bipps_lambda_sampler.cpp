#include "bipps.h"

using namespace std;

void Bipps::deal_with_Lambda(BippsDataLMC& data){
  sample_hmc_Lambda();
}

void Bipps::sample_hmc_Lambda(){
  if(verbose & debug){
    Rcpp::Rcout << "[sample_hmc_Lambda] starting\n";
  }
  
  start = std::chrono::steady_clock::now();
  
  arma::vec lambda_runif = vrunif(q);
  arma::vec lambda_runif2 = vrunif(q);
  for(unsigned int j=0; j<q; j++){
    arma::uvec subcols = arma::find(Lambda_mask.row(j) == 1);
    //Rcpp::Rcout << "j: " << j << endl;
    arma::vec offsets_obs = offsets(ix_by_q_a(j), oneuv * j);
    arma::vec xb_obs = XB(ix_by_q_a(j), oneuv * j);
    arma::vec offsets_for_lambda = offsets_obs + xb_obs;
    
    // build W
    // filter: choose value of spatial processes at locations of Yj that are available
    arma::mat WWj = w.submat(ix_by_q_a(j), subcols); // acts as X //*********
    //wmean.submat(ix_by_q_a(j), subcols); // acts as X
    
    arma::mat Wcrossprod = WWj.t() * WWj; 
    
    arma::mat Vi = 1 * arma::eye(WWj.n_cols, WWj.n_cols); 
    arma::vec Vim = arma::zeros(WWj.n_cols);
    
    lambda_node.at(j).update_mv(offsets_for_lambda, 1.0/tausq_inv(j), Vim, Vi);
    lambda_node.at(j).X = WWj;
    lambda_node.at(j).XtX = Wcrossprod;
    
    arma::vec curLrow = arma::trans(Lambda.submat(oneuv*j, subcols));
    arma::mat rnorm_row = mrstdnorm(curLrow.n_elem, 1);
    
    // nongaussian
    //Rcpp::Rcout << "step " << endl;
    lambda_hmc_adapt.at(j).step();
    if((lambda_hmc_started(j) == 0) && (lambda_hmc_adapt.at(j).i == 10)){
      // wait a few iterations before starting adaptation
      //Rcpp::Rcout << "reasonable stepsize " << endl;
      
      
      double lambda_eps = find_reasonable_stepsize(curLrow, lambda_node.at(j), rnorm_row);
      //Rcpp::Rcout << "adapting scheme starting " << endl;
      AdaptE new_adapting_scheme;
      new_adapting_scheme.init(lambda_eps, k, which_hmc, 1e4);
      lambda_hmc_adapt.at(j) = new_adapting_scheme;
      lambda_hmc_started(j) = 1;
      //Rcpp::Rcout << "done initiating adapting scheme" << endl;
    }
    
    arma::vec sampled;
    
    if(which_hmc == 0){
      // some form of manifold mala
      sampled = simpa_cpp(curLrow, lambda_node.at(j), lambda_hmc_adapt.at(j), 
                          rnorm_row, lambda_runif(j), lambda_runif2(j), debug);
    }
    if(which_hmc == 1){
      // mala
      sampled = mala_cpp(curLrow, lambda_node.at(j), lambda_hmc_adapt.at(j), 
                         rnorm_row, lambda_runif(j), debug);
    }
    if(which_hmc == 2){
      Rcpp::stop("HMC algorithm 2 not implemented");
    }
    
    if(which_hmc == 3){
      // some form of manifold mala
      sampled = smmala_cpp(curLrow, lambda_node.at(j), lambda_hmc_adapt.at(j), 
                 rnorm_row, lambda_runif(j), debug);
    }
    if(which_hmc == 6){
      sampled = hmc_cpp(curLrow, lambda_node.at(j), lambda_hmc_adapt.at(j), 
                 rnorm_row, lambda_runif(j), 
                 0.1, debug);
    }
    if(which_hmc == 7){
      Rcpp::stop("HMC algorithm 7 not implemented");
    }
    
    Lambda.submat(oneuv*j, subcols) = sampled.t();
    
    //Rcpp::Rcout << sampled.t();
  } 
  
  LambdaHw = w * Lambda.t();
  
  // refreshing density happens in the 'logpost_refresh_after_gibbs' function
  if(verbose & debug){
    Rcpp::Rcout << "[sample_hmc_Lambda] done\n";
  }
  
}
