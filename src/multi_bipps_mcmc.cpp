

#include "utils_lmc.h"
#include "utils_interrupt_handler.h"
#include "utils_parametrize.h"
#include "bipps.h"
#include "multi_bipps.h"

//[[Rcpp::export]]
Rcpp::List multi_bipps_mcmc(
    const arma::field<arma::mat>& y_list, // list of responses
    const arma::uvec& family,

    const arma::field<arma::mat>& X_list, // list of covariates

    const arma::mat& coords, // assuming that all coords are the same across samples (and all other derived features)

    int k,

    const arma::field<arma::uvec>& parents, // similar to coords, but for the graph. same for all samples.
    const arma::field<arma::uvec>& children,

    const arma::vec& layer_names,
    const arma::vec& layer_gibbs_group,


    const arma::field<arma::uvec>& indexing,

    const arma::mat& set_unif_bounds_in,
    const arma::mat& beta_Vi,

    int matern_twonu,

    const arma::field<arma::mat>& start_w,
    const arma::mat& lambda,
    const arma::umat& lambda_mask,
    const arma::mat& theta,
    const arma::mat& beta,
    const arma::vec& tausq,

    const arma::mat& mcmcsd,

    int mcmc_keep = 100,
    int mcmc_burn = 100,
    int mcmc_thin = 1,

    int mcmc_startfrom = 0,

    int num_threads = 1,

    int which_hmc=0,
    bool adapting=false,

    bool use_ps=true,

    bool verbose=false,
    bool debug=false,
    int print_every=false,
    bool low_mem=false,

    bool sample_beta=true,
    bool sample_tausq=true,
    bool sample_lambda=true,
    bool sample_theta=true,
    bool sample_w=true){

  if(verbose & debug){
    Rcpp::Rcout << "Initializing.\n";
  }

  
  bool sample_intercept = true;

#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#endif

  // timers
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point start_all = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point end_all = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point start_mcmc = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point end_mcmc = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point tick_mcmc = std::chrono::steady_clock::now();
  // ------

  bool printall = print_every == 1;
  bool verbose_mcmc = printall;

  double tempr = 1;

  unsigned int q  = y_list[0].n_cols;
  unsigned int n = y_list[0].n_rows;

  if(verbose & debug){
    Rcpp::Rcout << "Limits to MCMC search for theta:\n";
    Rcpp::Rcout << set_unif_bounds_in << endl;
  }
  // adaptive params
  int mcmc = mcmc_thin*mcmc_keep + mcmc_burn;

  arma::mat start_lambda = lambda; // *
    //ps_forward(theta, d, matern_twonu, use_ps);

  arma::mat start_theta = theta;

  arma::mat plus_icept_mcmc = arma::zeros(q, mcmc_thin*mcmc_keep);
  arma::cube b_mcmc = arma::zeros(X_list[0].n_cols, q, mcmc_thin*mcmc_keep);
  arma::mat tausq_mcmc = arma::zeros(q, mcmc_thin*mcmc_keep);
  arma::cube theta_mcmc = arma::zeros(theta.n_rows, k, mcmc_thin*mcmc_keep);
  arma::cube lambda_mcmc = arma::zeros(q, k, mcmc_thin*mcmc_keep);
  arma::cube vcov_mcmc = arma::zeros(k, k, mcmc_thin*mcmc_keep);

  int num_samples = y_list.n_elem;
  std::vector<Bipps> bipps_samples;
  bipps_samples.reserve(num_samples);

  for(unsigned int i = 0; i < num_samples; ++i) {
    Bipps bipps(
      y_list[i],
      family,
      X_list[i],
       coords, 
       i,
       k,
              parents, children, layer_names, layer_gibbs_group,

              indexing,

              matern_twonu,
              start_w[i], beta, start_lambda, lambda_mask, start_theta, 1.0/tausq,
              beta_Vi,

              which_hmc,
              adapting,
              mcmcsd,
              set_unif_bounds_in,
              use_ps,
              verbose, debug, num_threads);

    bipps_samples.push_back(bipps);
  }

  MultiBipps msp(
    bipps_samples,
    family,
    beta,
    start_lambda,
    lambda_mask, 
    start_theta, 
    1.0/tausq, 
    set_unif_bounds_in, 
    mcmcsd, 
    coords,
    beta_Vi,
    which_hmc,
    adapting,
    matern_twonu,
    use_ps,
    num_threads,
    verbose,
    debug
  );

  arma::cube icept_ss_mcmc = arma::zeros(msp.mb_size, msp.q, mcmc_thin*mcmc_keep);

  // what to do here?
  Rcpp::List caching_info;
  // caching_info["coords"] = msp.coords_caching.n_elem;
  // caching_info["hrmats"] = msp.kr_caching.n_elem;

  // save bipps index for long multi_bipps vectors
  arma::field<arma::uvec> image_ix(msp.mb_size);
  int Nn = msp.mb_size * n; 
  for(int i=0; i<msp.mb_size; i++){
    image_ix(i) = arma::regspace<arma::uvec>(n*i, n*(i+1)-1);
  }
  arma::mat v_temp = arma::zeros(Nn, k);
  // --
  
  
  arma::vec param = arma::vectorise(msp.multi_theta);

  arma::vec logaccept_mcmc = arma::zeros(mcmc);

  arma::uvec mcmc_ix = arma::zeros<arma::uvec>(mcmc_keep);
  arma::vec llsave = arma::zeros(mcmc_thin*mcmc_keep);
  arma::vec wllsave = arma::zeros(mcmc_thin*mcmc_keep);

  Rcpp::List v_mcmc;
  Rcpp::List w_mcmc;
  Rcpp::List lp_mcmc;
  Rcpp::List yhat_mcmc;

  for(int i=0; i<mcmc_keep; i++){
    if(!low_mem){
      std::string iname = std::to_string(i);
      v_mcmc[iname] = Rcpp::wrap(arma::zeros(msp.w_joined_rows, k));
      yhat_mcmc[iname] = Rcpp::wrap(arma::zeros(msp.y_joined_rows, q));
      w_mcmc[iname] = Rcpp::wrap(arma::zeros(msp.w_joined_rows, q));
      lp_mcmc[iname] = Rcpp::wrap(arma::zeros(msp.y_joined_rows, q));
    }
  }

  if(mcmc > 0){
    // what is going on here? just using this for the side effects I think.
    for(Bipps &bipps: msp.multi_bipps) {
      bipps.get_loglik_comps_w( bipps.param_data );
      bipps.get_loglik_comps_w( bipps.alter_data );
    }
  }



  double current_loglik = 0;
  for(Bipps &bipps: msp.multi_bipps) {
    current_loglik += tempr*bipps.param_data.loglik_w;
  }
  if(verbose & debug){
    Rcpp::Rcout << "Starting from logdens: " << current_loglik << endl;
  }

  double logaccept = 0;

  bool interrupted = false;

  if(verbose){
    Rcpp::Rcout << "Running MCMC for " << mcmc << " iterations.\n\n";
  }


  start_all = std::chrono::steady_clock::now();
  int m=0; int mx=0; int num_chol_fails=0;
  int mcmc_saved = 0; int w_saved = 0;

  try {

    for(m=0; (m<mcmc) & (!interrupted); m++){

      msp.predicting = false;
      mx = m-mcmc_burn;
      if(mx >= 0){
        if((mx % mcmc_thin) == 0){
          msp.predicting = true;
        }
      }

      if(printall){
        tick_mcmc = std::chrono::steady_clock::now();
      }

      if(sample_theta){
        start = std::chrono::steady_clock::now();
        msp.metrop_theta();
        end = std::chrono::steady_clock::now();
        if(verbose_mcmc & verbose){
          Rcpp::Rcout << "[theta] "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us.\n";
        }
      }

      if(sample_w){
        start = std::chrono::steady_clock::now();
        for(Bipps &bipps: msp.multi_bipps) {
          bipps.deal_with_w(bipps.param_data);
        }
        end = std::chrono::steady_clock::now();
        if(verbose_mcmc & verbose){
          Rcpp::Rcout << "[w] "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us.\n";
        }
        if(msp.predicting){
          start = std::chrono::steady_clock::now();
          for(Bipps &bipps: msp.multi_bipps) {
            bipps.predict();
          }
          end = std::chrono::steady_clock::now();
          if(verbose_mcmc & verbose){
            Rcpp::Rcout << "[predict] "
                        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us.\n";
          }
        }
      }
      
      if(sample_intercept){
        start = std::chrono::steady_clock::now();
        msp.sample_icept();
        end = std::chrono::steady_clock::now();
        if(verbose_mcmc & verbose){
          Rcpp::Rcout << "[sample_hmc_icept] "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us.\n";
        }
      }
      
      if(sample_lambda+sample_beta+sample_tausq){
        start = std::chrono::steady_clock::now();
        msp.sample_hmc_BetaLambdaTau(true, sample_beta, sample_lambda, sample_tausq); // true = sample
        end = std::chrono::steady_clock::now();
        if(verbose_mcmc & verbose){
          Rcpp::Rcout << "[BetaLambdaTau] "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us.\n";
        }
      }

      if(sample_tausq || sample_beta || sample_w || sample_lambda){
        start = std::chrono::steady_clock::now();
        for(Bipps &bipps: msp.multi_bipps) {
          bipps.logpost_refresh_after_gibbs(bipps.param_data);
        }
        end = std::chrono::steady_clock::now();
        if(verbose_mcmc & verbose){
          Rcpp::Rcout << "[logpost_refresh_after_gibbs] "
                      << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us.\n";
        }
      }

      //save
      logaccept_mcmc(m) = logaccept > 0 ? 0 : logaccept;

      double ll_joined = 0;
      double wll_joined = 0;
      for(Bipps &bipps: msp.multi_bipps) {
        ll_joined += bipps.logpost;
        wll_joined += bipps.param_data.loglik_w;
      }

      if(mx >= 0){
        tausq_mcmc.col(w_saved) = 1.0 / msp.tausq_inv;
        theta_mcmc.slice(w_saved) = msp.multi_theta;
        
        // calculate mean and cov of v
        // save identifiable parts
        // multi_ss_icept stores the subject-specific intercepts
        arma::mat multi_ss_icept = msp.multi_ss_icept;
        
        for(int im=0; im<msp.mb_size; im++){
          arma::vec plus_subj_icept = arma::trans(arma::mean(msp.multi_bipps[im].w, 0));
          for(int j=0; j<k; j++){
            v_temp(image_ix(im), msp.oneuv * j) = msp.multi_bipps[im].w.col(j) - plus_subj_icept(j); // subtract plus_subj_icept [1]
          }
          multi_ss_icept.row(im) += arma::trans(msp.multi_Lambda * plus_subj_icept); // add plus_subj_icept [1]
        } 
        arma::vec multi_ss_icept_means = arma::trans(arma::mean(multi_ss_icept, 0));
        for(int j=0; j<q; j++){
          multi_ss_icept.col(j) -= multi_ss_icept_means(j); // subtract multi_ss_icept_means [2]
        }
        
        arma::vec vmean = arma::trans(arma::mean(v_temp, 0));
        for(int j=0; j<k; j++){
          v_temp.col(j) -= vmean(j); // subtract vmean (global) [3]
        }
        
        arma::mat vcov = arma::cov(v_temp);
        vcov_mcmc.slice(w_saved) = vcov;
        
        //arma::mat U = arma::chol(vcov, "upper");
        //arma::mat Ui = arma::inv(arma::trimatu(U));
        //lambda_mcmc.slice(w_saved) = msp.multi_Lambda * U.t();
        
        lambda_mcmc.slice(w_saved) = msp.multi_Lambda;
        
        // global intercept might be here but we leave b_mcmc as is
        // the plus_icept_mcmc object stores the operation "add multi_ss_icept_means [2], add vmean (global) [3]"
        // which we can do in postprocessing if Beta includes the intercept
        b_mcmc.slice(w_saved) = msp.multi_Beta; 
        plus_icept_mcmc.col(w_saved) = msp.multi_Lambda * vmean + multi_ss_icept_means; 
        
        icept_ss_mcmc.slice(w_saved) = multi_ss_icept;
        
        llsave(w_saved) = ll_joined;
        wllsave(w_saved) = wll_joined;
        w_saved++;
        

        if(mx % mcmc_thin == 0){
            std::string iname = std::to_string(mcmc_saved);

            //arma::mat v_mcmc_joined; // replaced by v_temp
            arma::mat yhat_mcmc_joined;
            for(Bipps &bipps: msp.multi_bipps) {
              //v_mcmc_joined = arma::join_vert(v_mcmc_joined, bipps.w * H); 

              Rcpp::RNGScope scope;
              bipps.predicty();
              yhat_mcmc_joined = arma::join_vert(yhat_mcmc_joined, bipps.yhat);
            }

            v_mcmc[iname] = Rcpp::wrap(v_temp);
            yhat_mcmc[iname] = Rcpp::wrap(yhat_mcmc_joined);

          if(!low_mem){
            arma::mat w_mcmc_joined;
            arma::mat lp_mcmc_joined;
            for(Bipps &bipps: msp.multi_bipps) {
              w_mcmc_joined = arma::join_vert(w_mcmc_joined, bipps.LambdaHw); //** LambdaHw used here.
              lp_mcmc_joined = arma::join_vert(lp_mcmc_joined, bipps.linear_predictor);
            }
            w_mcmc[iname] = Rcpp::wrap(w_mcmc_joined);
            lp_mcmc[iname] = Rcpp::wrap(lp_mcmc_joined);
          }

          mcmc_ix(mcmc_saved) = w_saved;

          mcmc_saved++;
        }
      }

      interrupted = checkInterrupt();
      if(interrupted){
        Rcpp::stop("Interrupted by the user.");
      }

      if((m>0) & (mcmc > 100)){

        bool print_condition = (print_every>0);
        if(print_condition){
          print_condition = print_condition & (!(m % print_every));
        };

        if(print_condition){
          end_mcmc = std::chrono::steady_clock::now();

          int time_tick = std::chrono::duration_cast<std::chrono::milliseconds>(end_mcmc - tick_mcmc).count();
          int time_mcmc = std::chrono::duration_cast<std::chrono::milliseconds>(end_mcmc - start_mcmc).count();
          msp.theta_adapt.print_summary(time_tick, time_mcmc, m, mcmc);

          tick_mcmc = std::chrono::steady_clock::now();
          if(verbose & debug){
            Rprintf("  p(w|theta) = %.2f    p(y|...) = %.2f  \n ", wll_joined, ll_joined);
          }
          unsigned int printlimit = 10;

          msp.theta_adapt.print_acceptance();
          Rprintf("  theta = ");
          unsigned int n_theta = msp.multi_theta.n_elem;
          unsigned int n_print_theta = min(printlimit, n_theta);
          for(unsigned int pp=0; pp<n_print_theta; pp++){
            Rprintf("%.3f ", msp.multi_theta(pp));
          }


          if(arma::any(msp.familyid == 0)){
            Rprintf("\n  tausq = ");
            unsigned int n_print_tsq = min(printlimit, q);
            for(unsigned int pp=0; pp<n_print_tsq; pp++){
              if(msp.familyid(pp) == 0){
                Rprintf("(%d) %.6f ", pp+1, 1.0/msp.tausq_inv(pp));
              }
            }
          }
          if(arma::any(msp.familyid == 3)){
            Rprintf("\n  tau (betareg) = ");
            unsigned int n_print_tsq = min(printlimit, q);
            for(unsigned int pp=0; pp<n_print_tsq; pp++){
              if(msp.familyid(pp) == 3){
                Rprintf("(%d) %.4f ", pp+1, 1.0/msp.tausq_inv(pp));
              }
            }
          }
          if(arma::any(msp.familyid == 4)){
            Rprintf("\n  tau (negbinom) = ");
            unsigned int n_print_tsq = min(printlimit, q);
            for(unsigned int pp=0; pp<n_print_tsq; pp++){
              if(msp.familyid(pp) == 4){
                Rprintf("(%d) %.4f ", pp+1, 1.0/msp.tausq_inv(pp));
              }
            }
          }
          if(use_ps || q > 1){
            arma::vec lvec = arma::vectorise(msp.multi_Lambda);
            unsigned int n_lambda = lvec.n_elem;
            unsigned int n_print_lambda = min(printlimit, n_lambda);
            if(debug){
              Rprintf("\n  lambdastar = ");
              for(unsigned int pp=0; pp<n_print_lambda; pp++){
                Rprintf("%.3f ", lvec(pp));
              }
            }
            lvec = arma::vectorise(msp.multi_Lambda);
            Rprintf("\n  lambda = ");
            for(unsigned int pp=0; pp<n_print_lambda; pp++){
              Rprintf("%.3f ", lvec(pp));
            }
          }

          Rprintf("\n\n");
        }
      } else {
        tick_mcmc = std::chrono::steady_clock::now();
      }

    }

    end_all = std::chrono::steady_clock::now();
    double mcmc_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count();
    if(print_every>0){
      Rcpp::Rcout << "MCMC done [" << mcmc_time/1000.0 <<  "s]\n";
    }

    return Rcpp::List::create(
      Rcpp::Named("yhat_mcmc") = yhat_mcmc,
      Rcpp::Named("v_mcmc") = v_mcmc,
      Rcpp::Named("w_mcmc") = w_mcmc,
      Rcpp::Named("lp_mcmc") = lp_mcmc,
      Rcpp::Named("plus_icept_mcmc") = plus_icept_mcmc,
      Rcpp::Named("icept_ss_mcmc") = icept_ss_mcmc,
      Rcpp::Named("vcov_mcmc") = vcov_mcmc,
      Rcpp::Named("beta_mcmc") = b_mcmc,
      Rcpp::Named("tausq_mcmc") = tausq_mcmc,
      Rcpp::Named("theta_mcmc") = theta_mcmc,
      Rcpp::Named("lambda_mcmc") = lambda_mcmc,
      Rcpp::Named("paramsd") = msp.theta_adapt.paramsd,
      Rcpp::Named("mcmc") = mcmc,
      Rcpp::Named("mcmc_time") = mcmc_time/1000.0,
      Rcpp::Named("proposal_failures") = num_chol_fails,
      Rcpp::Named("caching_info") = caching_info,
      Rcpp::Named("mcmc_ix") = mcmc_ix,
      Rcpp::Named("success") = true
    );

  } catch(const std::exception& e) {
    Rcpp::Rcout << "Caught exception \"" << e.what() << "\"\n";

    end_all = std::chrono::steady_clock::now();

    double mcmc_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count();
    Rcpp::warning("MCMC has been interrupted. Returning partial saved results if any.\n");

    return Rcpp::List::create(
      Rcpp::Named("yhat_mcmc") = yhat_mcmc,
      Rcpp::Named("v_mcmc") = v_mcmc,
      Rcpp::Named("w_mcmc") = w_mcmc,
      Rcpp::Named("lp_mcmc") = lp_mcmc,
      Rcpp::Named("plus_icept_mcmc") = plus_icept_mcmc,
      Rcpp::Named("icept_ss_mcmc") = icept_ss_mcmc,
      Rcpp::Named("vcov_mcmc") = vcov_mcmc,
      Rcpp::Named("beta_mcmc") = b_mcmc,
      Rcpp::Named("tausq_mcmc") = tausq_mcmc,
      Rcpp::Named("theta_mcmc") = theta_mcmc,
      Rcpp::Named("lambda_mcmc") = lambda_mcmc,
      Rcpp::Named("paramsd") = msp.theta_adapt.paramsd,
      Rcpp::Named("mcmc") = mcmc,
      Rcpp::Named("mcmc_time") = mcmc_time/1000.0,
      Rcpp::Named("proposal_failures") = num_chol_fails,
      Rcpp::Named("caching_info") = caching_info,
      Rcpp::Named("mcmc_ix") = mcmc_ix,
      Rcpp::Named("success") = false
    );
  }
}

