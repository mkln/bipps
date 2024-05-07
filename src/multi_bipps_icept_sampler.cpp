#include "multi_bipps.h"
using namespace std;


void MultiBipps::sample_icept(){
  for(int i=0; i<mb_size; i++) {
    multi_bipps[i].sample_hmc_icept();
    multi_ss_icept.row(i) = arma::trans(multi_bipps[i].icept);
  }
}
