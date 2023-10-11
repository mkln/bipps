#include "bipps.h"

using namespace std;

void Bipps::deal_with_beta(bool sample){
  hmc_sample_beta(sample);
}

void Bipps::hmc_sample_beta(bool sample){
  if(verbose & debug){
    Rcpp::Rcout << "[hmc_sample_beta]\n";
  }
  
  // choose between NUTS or Gibbs
  start = std::chrono::steady_clock::now();
  
  arma::mat bmat = mrstdnorm(p, q);
  arma::vec bunifv = vrunif(q);
  arma::vec bunifv2 = vrunif(q);
  
  arma::mat LHW = w * Lambda.t();
  
  for(unsigned int j=0; j<q; j++){
    //Rcpp::Rcout << "j: " << j << endl;
    arma::vec offsets_obs = offsets(ix_by_q_a(j), oneuv * j);
    arma::vec lw_obs = LHW(ix_by_q_a(j), oneuv * j);
    arma::vec offsets_for_beta = offsets_obs + lw_obs;
    
    beta_node.at(j).update_mv(offsets_for_beta, 1.0/tausq_inv(j), Vim, Vi);

    // nongaussian
    beta_hmc_adapt.at(j).step();
    if(sample & (beta_hmc_started(j) == 0)){
      // wait a few iterations before starting adaptation
      //Rcpp::Rcout << "reasonable stepsize " << endl;
      double beta_eps = find_reasonable_stepsize(Bcoeff.col(j), beta_node.at(j), bmat.cols(oneuv * j));
      //Rcpp::Rcout << "adapting scheme starting " << endl;
      AdaptE new_adapting_scheme;
      new_adapting_scheme.init(beta_eps, p, which_hmc, 1e4);
      beta_hmc_adapt.at(j) = new_adapting_scheme;
      beta_hmc_started(j) = 1;
      //Rcpp::Rcout << "done initiating adapting scheme" << endl;
    }
    
    if(which_hmc == 0){
      // some form of manifold mala
      Bcoeff.col(j) = simpa_cpp(Bcoeff.col(j), beta_node.at(j), 
                          beta_hmc_adapt.at(j), bmat.cols(oneuv * j), bunifv(j), 
                          bunifv2(j), debug);
    }
    if(which_hmc == 1){
      // mala
      Bcoeff.col(j) = mala_cpp(Bcoeff.col(j), beta_node.at(j), 
                         beta_hmc_adapt.at(j), bmat.cols(oneuv * j), bunifv(j), 
                         debug);
    }
    if(which_hmc == 2){
      Rcpp::stop("HMC algorithm 2 not implemented");
    }
    
    if(which_hmc == 3){
      // some form of manifold mala
      Bcoeff.col(j) = smmala_cpp(Bcoeff.col(j), beta_node.at(j), 
                           beta_hmc_adapt.at(j), bmat.cols(oneuv * j), 
                           bunifv(j),  debug);
    }
    if(which_hmc == 6){
      Bcoeff.col(j) = hmc_cpp(Bcoeff.col(j), beta_node.at(j), 
                        beta_hmc_adapt.at(j), bmat.cols(oneuv * j), bunifv(j), 
                        0.1, debug);
    }
    if(which_hmc == 7){
      Rcpp::stop("HMC algorithm 7 not implemented");
    }
  
    XB.col(j) = X * Bcoeff.col(j);
  }
  
  if(verbose){
    end = std::chrono::steady_clock::now();
    Rcpp::Rcout << "[hmc_sample_beta] "
                << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
                << "us.\n";
  }
}


