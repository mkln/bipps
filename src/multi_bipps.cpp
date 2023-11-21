#include "multi_bipps.h"
using namespace std;

MultiBipps::MultiBipps(
  const std::vector<Bipps>& multi_bipps_in,
  const arma::vec& tausq_inv_in,
  const arma::mat& theta_in,
  const arma::mat& metrop_theta_bounds,
  const arma::mat& metrop_theta_sd,
  const arma::mat& coords_in,
  bool adapting_theta,
  bool use_ps=true,
  int num_threads=1,
  int matern_twonu_in
  ) {
	  multi_bipps = multi_bipps_in;

    mb_size = multi_bipps.size();
    lambda_hmc_started = arma::zeros<arma::uvec>(q);
    tausq_inv = tausq_inv_in;

    coords = coords_in;

    dd = coords.n_cols;

    theta_unif_bounds = metrop_theta_bounds;
    theta_metrop_sd = metrop_theta_sd;
    theta_adapt = RAMAdapt(theta_in.n_elem, theta_metrop_sd, 0.24);
    theta_adapt_active = adapting_theta;
    theta_mcmc_counter = 0;

    predicting = false;

    w_joined_rows = 0;
    for(Bipps& bipps : multi_bipps) {
      w_joined_rows += bipps.w.n_rows;
    }

    y_joined_rows = 0;
    for(Bipps& bipps : multi_bipps) {
      y_joined_rows += bipps.y.n_rows;
    }

    init_matern(num_threads, matern_twonu_in, use_ps);

    for(unsigned int j=0; j<q; j++){
      int family = familyid(j);
      arma::mat Xj_joined;
      arma::vec offsets_joined;
      arma::vec yj_joined;
      for(Bipps& bipps : multi_bipps) {

        arma::mat LHW = bipps.w * bipps.Lambda.t();

        arma::vec yj_obs = bipps.y( bipps.ix_by_q_a(j), oneuv * j );
        arma::mat X_obs = bipps.X.rows(bipps.ix_by_q_a(j));
        arma::mat offsets_obs = bipps.offsets(bipps.ix_by_q_a(j), oneuv * j);
        arma::vec lw_obs = LHW(bipps.ix_by_q_a(j), oneuv * j);
        
        arma::vec offsets_for_beta = offsets_obs + lw_obs;

        Xj_joined = arma::join_vert(Xj_joined, X_obs);
        offsets_joined = arma::join_vert(offsets_joined, offsets_obs);
        yj_joined = arma::join_vert(yj_joined, yj_obs);
      }
        
        // Lambda
        NodeDataB new_lambda_block(yj_joined, offsets_joined, Xj_joined, family);
        lambda_node.push_back(new_lambda_block);
        
        // *** sampling beta and lambda together so we use p+k here
        arma::uvec subcols = arma::find(Lambda_mask.row(j) == 1);
        int n_lambdas = subcols.n_elem;
        AdaptE new_lambda_adapt;
        new_lambda_adapt.init(.05, p+n_lambdas, which_hmc);
        lambda_hmc_adapt.push_back(new_lambda_adapt);
    }
  };

void MultiBipps::accept_make_change(){
  for(Bipps& bipps : multi_bipps) {
    bipps.accept_make_change();
  }
}

void MultiBipps::init_matern(int num_threads, int matern_twonu_in=1, bool use_ps=true){
  nThreads = num_threads;
  
  int bessel_ws_inc = 5;
  matern.bessel_ws = (double *) R_alloc(nThreads*bessel_ws_inc, sizeof(double));
  matern.twonu = matern_twonu_in;
  matern.using_ps = use_ps;
  matern.estimating_nu = (dd == 2) & (multi_theta.n_rows == 3);
  
}
