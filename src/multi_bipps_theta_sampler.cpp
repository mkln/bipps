#include "multi_bipps.h"

using namespace std;

void MultiBipps::metrop_theta(){
  if(verbose & debug){
    Rcpp::Rcout << "[metrop_theta] start\n";
  }
  
  theta_adapt.count_proposal();
  
  arma::vec param = arma::vectorise(multi_theta);
  arma::vec new_param = arma::vectorise(multi_theta);
  
  Rcpp::RNGScope scope;
  arma::vec U_update = mrstdnorm(new_param.n_elem, 1);
  
  // theta
  new_param = par_huvtransf_back(par_huvtransf_fwd(param, theta_unif_bounds) + 
    theta_adapt.paramsd * U_update, theta_unif_bounds);
  
  arma::mat theta_proposal = 
    arma::mat(new_param.memptr(), new_param.n_elem/k, k);

  std::vector<bool> acceptable_joined;

  double current_loglik = 0;
  double new_loglik = 0;

  // run get_loglik_comps_w on first bipps object, then copy relevant data to all others.

  Bipps &first_bipps = multi_bipps[0];
  first_bipps.alter_data.theta = theta_proposal;
  bool first_acceptable = first_bipps.get_loglik_comps_w( first_bipps.alter_data );

  if(first_acceptable) {
    for(Bipps &bipps: multi_bipps) {
      bipps.alter_data.CC_cache = first_bipps.alter_data.CC_cache;
      bipps.alter_data.Kxxi_cache = first_bipps.alter_data.Kxxi_cache;
      bipps.alter_data.Ri_chol_logdet = first_bipps.alter_data.Ri_chol_logdet;
      bipps.alter_data.Ri_cache = first_bipps.alter_data.Ri_cache;
      bipps.alter_data.H_cache = first_bipps.alter_data.H_cache;
      bipps.alter_data.Kppi_cache = first_bipps.alter_data.Kppi_cache;

      bool acceptable = bipps.calc_ywlogdens(bipps.alter_data);
      acceptable_joined.push_back(acceptable);

      new_loglik += bipps.alter_data.loglik_w;
      current_loglik += bipps.param_data.loglik_w;
    }
  }
  
  // loop over all logliks and add together?
  
  bool accepted = false;
  double logaccept = 0;
  double prior_logratio = 0;
  double jacobian = 0;

  bool acceptable = std::all_of(acceptable_joined.begin(), acceptable_joined.end(), [](bool v) { return v; });
  
  if(acceptable){ 
    // stay the same
    prior_logratio = calc_prior_logratio(
        theta_proposal.tail_rows(1).t(), multi_theta.tail_rows(1).t(), 2, 1); // sigmasq
    
    if(multi_theta.n_rows > 5){
      for(auto i=0; i<multi_theta.n_rows-2; i++){
        prior_logratio += arma::accu( -theta_proposal.row(i) +multi_theta.row(i) ); // exp
      }
    }
    
    // stay the same
    jacobian  = calc_jacobian(new_param, param, theta_unif_bounds);
    logaccept = new_loglik - current_loglik + 
      prior_logratio +
      jacobian;
    
    
    accepted = do_I_accept(logaccept);
    
  } else {
    accepted = false;
    //num_chol_fails ++;
    if(verbose & debug){
      Rcpp::Rcout << "[warning] numerical failure at MH proposal -- auto rejected\n";
    }
  }
  
  if(accepted){
    theta_adapt.count_accepted();
    
    accept_make_change();
    multi_theta = theta_proposal;
    for(Bipps &bipps: multi_bipps) {
      bipps.param_data.theta = theta_proposal;
    }
    
    if(debug & verbose){
      Rcpp::Rcout << "[theta] accepted (log accept. " << logaccept << " : " << new_loglik << " " << current_loglik << 
        " " << prior_logratio << " " << jacobian << ")\n";
    }
  } else {
    if(debug & verbose){
      Rcpp::Rcout << "[theta] rejected (log accept. " << logaccept << " : " << new_loglik << " " << current_loglik << 
        " " << prior_logratio << " " << jacobian << ")\n";
    }
  }
  
  
  theta_adapt.update_ratios();
  
  if(theta_adapt_active){
    theta_adapt.adapt(U_update, acceptable*exp(logaccept), theta_mcmc_counter); 
  }
  theta_mcmc_counter++;
  if(verbose & debug){
    Rcpp::Rcout << "[metrop_theta] end\n";
  }
}
