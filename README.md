## Bayesian Inference of Point Patterns in Space 

<!-- badges: start -->
[![R-CMD-check](https://github.com/mkln/bipps/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mkln/bipps/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

### Install: 

`devtools::install_github("mkln/bipps")` installs from this repository.

#### Tips for best performance:

 - `bipps` works best with OpenMP and OpenBLAS or Intel MKL. 
 - [Dirk Eddelbuettel has a great guide on installing Intel MKL on Debian/Ubuntu systems](http://dirk.eddelbuettel.com/blog/2018/04/15/#018_mkl_for_debian_ubuntu). In that case it is important to add `MKL_THREADING_LAYER=GNU` to `~/.Renviron`. 
 - On systems with AMD CPUs, it may be best to install `intel-mkl-2019.5-075` and then also add the line `MKL_DEBUG_CPU_TYPE=5` to `~/.Renviron`. I have not tested more recent versions of Intel MKL.
 - If using OpenBLAS, it might be important to let OpenMP do *all* the parallelization when running `bipps`. I think this can be done with the [RhpcBLASctl](https://CRAN.R-project.org/package=RhpcBLASctl) package. YMMV.

