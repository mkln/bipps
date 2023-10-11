#include "bipps.h"

using namespace std;

void Bipps::w_prior_sample(BippsDataLMC& data){
  if(verbose & debug){
    Rcpp::Rcout << "[w_prior_sample] " << "\n";
  }
  //Rcpp::Rcout << "Lambda from:  " << Lambda_orig(0, 0) << " to  " << Lambda(0, 0) << endl;
  
  start_overall = std::chrono::steady_clock::now();
  
  //int ns = coords.n_rows;
  
  bool acceptable = refresh_cache(data);
  if(!acceptable){
    Rcpp::stop("Something went wrong went getting the conditional Gaussians. Try different theta? ");
  }
  
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(unsigned int i = 0; i<n_blocks; i++){
    int u = block_names(i)-1;
    update_block_covpars(u, data);
  }
  
  // assuming that the ordering in block_names is the ordering of the product of conditional densities
  for(unsigned int i=0; i<n_blocks; i++){
    
    int u = block_names(i) - 1;
    
    // recompute conditional mean
    arma::mat Sigi_tot = (*data.w_cond_prec_ptr.at(u)).slice(0); //
    
    arma::mat w_mean = arma::zeros(indexing(u).n_elem);
    if(parents(u).n_elem > 0){
      w_mean = (*data.w_cond_mean_K_ptr.at(u)).slice(0) * w.rows(parents_indexing(u));
    } 
    arma::mat Sigi_chol = arma::inv(arma::trimatl(arma::chol(Sigi_tot, "lower")));
    
    // sample
    arma::vec rnvec = arma::randn(indexing(u).n_elem);
    arma::vec wtemp = w_mean + Sigi_chol.t() * rnvec;
    
    w.rows(indexing(u)) = wtemp;
    
  }
  
  if(verbose & debug){
    end_overall = std::chrono::steady_clock::now();
    
    Rcpp::Rcout << "[w_prior_sample] loops "
                << std::chrono::duration_cast<std::chrono::microseconds>(end_overall - start_overall).count()
                << "us. " << "\n";
  }
}

