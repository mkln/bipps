
#include <RcppArmadillo.h>
#include "bipps.h"

class MultiBipps {
public:
  std::vector<Bipps> multi_bipps;
  arma::mat multi_Lambda;
  arma::mat multi_Beta;
  arma::mat multi_theta;

  bool verbose;
  bool debug;

  std::chrono::steady_clock::time_point start;

  arma::vec tausq_inv;

  int mb_size;

  int q;
  int p;
  int k;
  int dd;
  int which_hmc;

  arma::mat coords;

  arma::vec y_joined;
  arma::mat Lambda_mask;
  arma::uvec oneuv;
  arma::uvec familyid;
  arma::mat Vi;

  std::vector<NodeDataB> lambda_node; // std::vector
  std::vector<AdaptE> lambda_hmc_adapt; // std::vector
  arma::uvec lambda_hmc_started;

  arma::mat theta_metrop_sd;

  RAMAdapt theta_adapt;
  arma::mat theta_unif_bounds;

  int nThreads;
  MaternParams matern;
  void init_matern(int num_threads, int matern_twonu_in, bool use_ps);

  void sample_hmc_BetaLambdaTau(bool sample, bool sample_beta, bool sample_lambda, bool sample_tau);
  void metrop_theta();

  void accept_make_change();
  bool theta_adapt_active;
  int theta_mcmc_counter;

  bool predicting;

  int w_joined_rows;
  int y_joined_rows;
  // constructors

  MultiBipps(){};
  MultiBipps(
    const std::vector<Bipps>& multi_bipps_in,
    const arma::vec& tausq_inv_in,
    const arma::mat& theta_in,
    const arma::mat& metrop_theta_bounds,
    const arma::mat& metrop_theta_sd,
    const arma::mat& coords_in,
    bool adapting_theta,
    bool use_ps,
    int num_threads,
    int matern_twonu_in
    );
};
