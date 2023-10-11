#include "bipps.h"

using namespace std;

void Bipps::deal_with_w(BippsDataLMC& data, bool sample){
  if(sample){
    // ensure this is independent on the number of threads being used
    Rcpp::RNGScope scope;
    rand_norm_mat = mrstdnorm(coords.n_rows, k);
    rand_unif = vrunif(n_blocks);
    rand_unif2 = vrunif(n_blocks);
  }
  
  nongaussian_w(data, sample);
}
  
void Bipps::update_block_w_cache(int u, BippsDataLMC& data){
  // 
  arma::mat Sigi_tot = build_block_diagonal_ptr(data.w_cond_prec_ptr.at(u));
  arma::mat Smu_tot = arma::zeros(k*indexing(u).n_elem, 1); // replace with fill(0)
  
  for(unsigned int c=0; c<children(u).n_elem; c++){
    int child = children(u)(c);
    arma::cube AK_u = cube_cols_ptr(data.w_cond_mean_K_ptr.at(child), u_is_which_col_f(u)(c)(0));
    AKuT_x_R_ptr(data.AK_uP(u)(c), AK_u, data.w_cond_prec_ptr.at(child)); 
    add_AK_AKu_multiply_(Sigi_tot, data.AK_uP(u)(c), AK_u);
  }

  arma::mat u_tau_inv = arma::zeros(indexing(u).n_elem, q);
  arma::mat ytilde = arma::zeros(indexing(u).n_elem, q);
  
  for(unsigned int j=0; j<q; j++){
    for(unsigned int ix=0; ix<indexing(u).n_elem; ix++){
      if(na_mat(indexing(u)(ix), j) == 1){
        u_tau_inv(ix, j) = pow(tausq_inv(j), .5);
        ytilde(ix, j) = (y(indexing(u)(ix), j) - XB(indexing(u)(ix), j))*u_tau_inv(ix, j);
      }
    }
    // dont erase:
    //Sigi_tot += arma::kron( arma::trans(Lambda.row(j)) * Lambda.row(j), arma::diagmat(u_tau_inv%u_tau_inv));
    arma::mat LjtLj = arma::trans(Lambda.row(j)) * Lambda.row(j);
    arma::vec u_tausq_inv = u_tau_inv.col(j) % u_tau_inv.col(j);
    add_LtLxD(Sigi_tot, LjtLj, u_tausq_inv);
    
    Smu_tot += arma::vectorise(arma::diagmat(u_tau_inv.col(j)) * ytilde.col(j) * Lambda.row(j));
  }
  
  data.Smu_start(u) = Smu_tot;
  data.Sigi_chol(u) = Sigi_tot;
  
  data.Sigi_chol(u) = arma::inv(arma::trimatl(arma::chol( arma::symmatu( Sigi_tot ), "lower")));
}

void Bipps::refresh_w_cache(BippsDataLMC& data){
  if(verbose & debug){
    Rcpp::Rcout << "[refresh_w_cache] \n";
  }
  start_overall = std::chrono::steady_clock::now();
  for(unsigned int i=0; i<n_blocks; i++){
    int u=block_names(i)-1;
    update_block_w_cache(u, data);
  }
  
  if(verbose & debug){
    end_overall = std::chrono::steady_clock::now();
    Rcpp::Rcout << "[refresh_w_cache] "
                << std::chrono::duration_cast<std::chrono::microseconds>(end_overall - start_overall).count()
                << "us. " << "\n";
  }
}

void Bipps::nongaussian_w(BippsDataLMC& data, bool sample){
  if(verbose & debug){
    Rcpp::Rcout << "[hmc_sample_w] " << endl;
  }
  
  start_overall = std::chrono::steady_clock::now();
  //double part1 = 0;
  //double part2 = 0;
  
  int mala_timer = 0;
  
  arma::mat offset_for_w = offsets + XB;
  
  for(int g=0; g<n_gibbs_groups; g++){
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(unsigned int i=0; i<u_by_block_groups(g).n_elem; i++){
      int u = u_by_block_groups(g)(i);
      
      if((block_ct_obs(u) > 0)){
        
        start = std::chrono::steady_clock::now();
        //Rcpp::Rcout << "u :  " << u << endl;
        w_node.at(u).update_mv(offset_for_w, 1.0/tausq_inv, Lambda);
        
        if(parents(u).n_elem > 0){
          //w_node.at(u).Kxxi = (*data.w_cond_prec_parents_ptr.at(u));
          w_node.at(u).Kcx = data.w_cond_mean_K_ptr.at(u);
          //w_node.at(u).w_parents = w.rows(parents_indexing(u));
        }
        
        w_node.at(u).Ri = data.w_cond_prec_ptr.at(u);
        
        if(parents_indexing(u).n_rows > 0){
          w_node.at(u).Kcxpar = arma::zeros((*w_node.at(u).Kcx).n_rows, k);
          for(unsigned int j=0; j<k; j++){
            w_node.at(u).Kcxpar.col(j) = (*w_node.at(u).Kcx).slice(j) * w.submat(parents_indexing(u), oneuv*j);//w_node.at(u).w_parents.col(j);
          }
        } else {
          w_node.at(u).Kcxpar = arma::zeros(0,0);
        }
          
        for(unsigned int c=0; c<w_node.at(u).num_children; c++ ){
          int child = children(u)(c);
          //Rcpp::Rcout << "child [" << child << "]\n";
          arma::uvec c_ix = indexing(child);
          arma::uvec pofc_ix = parents_indexing(child);
          
          arma::uvec pofc_ix_x = u_is_which_col_f(u)(c)(0);
          arma::uvec pofc_ix_other = u_is_which_col_f(u)(c)(1);
          
          arma::mat w_childs_parents = w.rows(pofc_ix);
          arma::mat w_otherparents = w_childs_parents.rows(pofc_ix_other);
          
          w_node.at(u).w_child(c) = w.rows(c_ix);
          w_node.at(u).Ri_of_child(c) = data.w_cond_prec_ptr.at(child);
          
          for(unsigned int r=0; r<k; r++){
            w_node.at(u).Kcx_x(c).slice(r) = (*data.w_cond_mean_K_ptr.at(child)).slice(r).cols(pofc_ix_x);
          }
          
          //Rcpp::Rcout << "hmc_sample_w " << arma::size(w_node.at(u).Kcx_x(c)) << endl;
          
          if(w_otherparents.n_rows > 0){
            //arma::cube Kcx_other = //(*param_data.w_cond_mean_K.at(child)).cols(pofc_ix_other);
            //  cube_cols_ptr(data.w_cond_mean_K_ptr.at(child), pofc_ix_other);
            
            w_node.at(u).Kco_wo(c) = arma::zeros(c_ix.n_elem, k); //cube_times_mat(Kcx_other, w_otherparents);
            for(unsigned int r=0; r<k; r++){
              w_node.at(u).Kco_wo(c).col(r) = (*data.w_cond_mean_K_ptr.at(child)).slice(r).cols(pofc_ix_other) * w_otherparents.col(r);
            }
            
          } else {
            w_node.at(u).Kco_wo(c) = arma::zeros(0,0);
          }
          //Rcpp::Rcout << "child done " << endl;
        } 
        
        // -----------------------------------------------
          
        arma::mat w_current = w.rows(indexing(u));
        
        // adapting eps
        hmc_eps_adapt.at(u).step();
        
        if((hmc_eps_started_adapting(u) == 0) && (hmc_eps_adapt.at(u).i==10) & sample){
          // wait a few iterations before starting adaptation
          //message("find reasonable");
          hmc_eps(u) = find_reasonable_stepsize(w_current, w_node.at(u), rand_norm_mat.rows(indexing(u)));
          //message("found reasonable");
          int blocksize = indexing(u).n_elem * k;
          AdaptE new_adapting_scheme;
          new_adapting_scheme.init(hmc_eps(u), blocksize, which_hmc, 1e4);
          
          hmc_eps_adapt.at(u) = new_adapting_scheme;
          hmc_eps_started_adapting(u) = 1;
        }
        
        
        bool do_gibbs = arma::all(familyid == 0);
        arma::mat w_temp = w_current;
       
        if(which_hmc == 0){
          w_temp = simpa_cpp(w_current, w_node.at(u), hmc_eps_adapt.at(u),
                                 rand_norm_mat.rows(indexing(u)),
                                 rand_unif(u), rand_unif2(u),
                                 debug);
        }
        if(which_hmc == 1){
          // mala
          w_temp = mala_cpp(w_current, w_node.at(u), hmc_eps_adapt.at(u),
                                       rand_norm_mat.rows(indexing(u)), rand_unif(u),
                                       debug);
        }
        if(which_hmc == 2){
          Rcpp::stop("HMC algorithm 2 not implemented");
        }
        if((which_hmc == 3) || (which_hmc == 4)){
          // some form of manifold mala
          w_temp = smmala_cpp(w_current, w_node.at(u), hmc_eps_adapt.at(u),
                                       rand_norm_mat.rows(indexing(u)),
                                       rand_unif(u), debug);
        }
        if(which_hmc == 5){
          Rcpp::stop("ESS not implemented");
        }
        if(which_hmc == 6){
          w_temp = hmc_cpp(w_current, w_node.at(u), hmc_eps_adapt.at(u),
                           rand_norm_mat.rows(indexing(u)),
                           rand_unif(u), 0.1, debug);
        }
        if(which_hmc == 7){
          Rcpp::stop("HMC algorithm 7 not implemented");
        }
      
      
        end = std::chrono::steady_clock::now();
        
        hmc_eps(u) = hmc_eps_adapt.at(u).eps;
        
        w.rows(indexing(u)) = w_temp;//arma::trans(arma::mat(w_temp.memptr(), q, w_temp.n_elem/q));
        //wU.rows(indexing(u)) = w.rows(indexing(u));

        mala_timer += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //Rcpp::Rcout << "done sampling "<< endl;

      }
    }
  }
  
  LambdaHw = w * Lambda.t();
  
  if(verbose & debug){
    end_overall = std::chrono::steady_clock::now();
    Rcpp::Rcout << "[hmc_sample_w] "
                << std::chrono::duration_cast<std::chrono::microseconds>(end_overall - start_overall).count()
                << "us. " << "\n";
    //Rcpp::Rcout << "mala: " << mala_timer << endl;
  }
  
}


void Bipps::predict(bool sample){
  start_overall = std::chrono::steady_clock::now();
  if(predict_group_exists == 1){
    if(verbose & debug){
      Rcpp::Rcout << "[predict] start \n";
    }
    
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(unsigned int i=0; i<u_predicts.n_elem; i++){ //*** subset to blocks with NA
      int u = u_predicts(i);// u_predicts(i);
      // only predictions at this block. 
      arma::uvec predict_parent_indexing, cx;
      arma::cube Kxxi_parents;
      
      // no observed locations, use line of sight
      predict_parent_indexing = parents_indexing(u);
      CviaKron_HRj_chol_bdiag(Hpred(i), Rcholpred(i), Kxxi_parents,
                              na_1_blocks(u),
                              coords, indexing(u), predict_parent_indexing, k, param_data.theta, matern);
    
      
      //Rcpp::Rcout << "step 1 "<< endl;
      arma::mat wpars = w.rows(predict_parent_indexing);
      
      for(unsigned int ix=0; ix<indexing(u).n_elem; ix++){
        if(na_1_blocks(u)(ix) == 0){
          arma::rowvec wtemp = arma::sum(arma::trans(Hpred(i).slice(ix)) % wpars, 0);
          
          if(sample){
            wtemp += arma::trans(Rcholpred(i).col(ix)) % rand_norm_mat.row(indexing(u)(ix));
          }
          
          w.row(indexing(u)(ix)) = wtemp;
          //wU.row(indexing(u)(ix)) = wtemp;
          
          LambdaHw.row(indexing(u)(ix)) = w.row(indexing(u)(ix)) * Lambda.t();
        }
      }
      
      //Rcpp::Rcout << "done" << endl;
      
      
    }
    
    if(verbose & debug){
      end_overall = std::chrono::steady_clock::now();
      Rcpp::Rcout << "[predict] "
                  << std::chrono::duration_cast<std::chrono::microseconds>(end_overall - start_overall).count()
                  << "us. ";
    }
  }
}

void Bipps::predicty(bool sample){  
  int n = XB.n_rows;
  yhat.fill(0);
  Rcpp::RNGScope scope;
  arma::mat Lw = w*Lambda.t();
  for(unsigned int j=0; j<q; j++){
    linear_predictor.col(j) = XB.col(j) + Lw.col(j);
    if(familyid(j) == 0){
      // gaussian
      yhat.col(j) = linear_predictor.col(j);
      if(sample){
        yhat.col(j) += pow(1.0/tausq_inv(j), .5) * mrstdnorm(n, 1);
      }
    } else if(familyid(j) == 1){
      // poisson
      if(sample){
        yhat.col(j) = vrpois(exp(linear_predictor.col(j)));
      } else {
        yhat.col(j) = exp(linear_predictor.col(j));
      }
    } else if(familyid(j) == 2){
      // binomial
      if(sample){
        yhat.col(j) = vrbern(1.0/(1.0 + exp(-linear_predictor.col(j))));
      } else {
        yhat.col(j) = 1.0/(1.0 + exp(-linear_predictor.col(j)));
      }
    } else if(familyid(j) == 3){
      arma::vec mu =  1.0/ (1.0 + exp(-linear_predictor.col(j)));
      arma::vec aa = tausq_inv(j) * mu;
      arma::vec bb = tausq_inv(j) * (1.0-mu);
      yhat.col(j) = vrbeta(aa, bb);
    } else if(familyid(j) == 4){
      arma::vec muvec = exp(linear_predictor.col(j));
      double alpha = 1.0/tausq_inv(j); 
      yhat.col(j) = vrnb(muvec, alpha);
    }
  }
}


