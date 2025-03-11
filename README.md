# BIPPS: Bayesian Inference for Point Patterns in Space

<!-- badges: start -->
[![R-CMD-check](https://github.com/mkln/bipps/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mkln/bipps/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

BIPPS (Bayesian Inference for Point Patterns in Space) is a spatial statistical modeling package for analyzing point patterns or count data across a spatial domain. It uses latent Gaussian processes to model the underlying intensity of a point process, making it particularly useful for applications in:

- Ecological data analysis
- Disease mapping
- Astronomical observations
- Any spatial point pattern or count data analysis

The package implements a fully Bayesian approach with efficient MCMC sampling, allowing for inference on complex spatial dependencies through a latent factor model structure.

## Features

- Model spatial point processes using latent Gaussian processes
- Support for multiple response types (multivariate spatial data)
- Configurable number of latent factors
- Support for multiple data replicates/images
- Parallel processing capabilities
- Various spatial correlation structures
- Support for different distributional families (Poisson, Negative Binomial)
- Implementation of spatial cross-correlation functions

## Installation

You can install the development version of BIPPS from GitHub:

```r
devtools::install_github("mkln/bipps")
```

### Performance Optimization

For best performance, consider the following:

- **OpenMP and BLAS Libraries**: BIPPS works best with OpenMP and optimized BLAS libraries like OpenBLAS or Intel MKL.

- **Intel MKL**: [Dirk Eddelbuettel has a guide on installing Intel MKL on Debian/Ubuntu systems](http://dirk.eddelbuettel.com/blog/2018/04/15/#018_mkl_for_debian_ubuntu). Add `MKL_THREADING_LAYER=GNU` to your `~/.Renviron` file.

- **AMD CPUs**: When using AMD CPUs, it may be best to install `intel-mkl-2019.5-075` and add `MKL_DEBUG_CPU_TYPE=5` to your `~/.Renviron`. More recent versions of Intel MKL have not been extensively tested.

- **OpenBLAS**: If using OpenBLAS, it's recommended to let OpenMP handle all parallelization when running BIPPS. This can be managed with the [RhpcBLASctl](https://CRAN.R-project.org/package=RhpcBLASctl) package.

## Getting Started

Here's a basic example of how to use BIPPS:

```r
library(bipps)
library(dplyr)
library(ggplot2)

# Set random seed for reproducibility
set.seed(42)

# Create spatial grid
nx <- 30
ny <- 30
coords <- expand.grid(
  x = seq(0, 1, length.out = nx),
  y = seq(0, 1, length.out = ny)
)

# Example usage with synthetic data
# (See vignettes for complete examples)
results <- multi_bipps(
  y_list = y_list,         # Your response data (point counts)
  x_list = x_list,         # Your covariates
  coords = coords,         # Spatial coordinates
  k = 3,                   # Number of latent factors
  family = "poisson",      # Response distribution
  n_samples = 1000,        # MCMC samples to collect
  n_burn = 500,            # Burn-in iterations
  n_threads = 2,           # Parallel threads
  verbose = 10             # Progress reporting interval
)

# Access results
summary(results)
```

## Documentation

For more detailed examples and usage information, please refer to the package vignettes:

```r
# Install package with vignettes
devtools::install_github("mkln/bipps", build_vignettes = TRUE)

# View available vignettes
vignette(package = "bipps")
```

## Model Description

BIPPS implements spatial models where:

- The data follow a point process with spatially varying intensity
- Spatial dependence is modeled through latent Gaussian processes
- A factor model structure allows for dimensionality reduction
- Spatial cross-correlations describe relationships between different response types

The spatial cross-correlation function between response types r and s is defined as:

$$
\rho_{rs}(\mathbf{h}; \mathbf{\theta}) = \frac{C_{rs}(\mathbf{h}; \mathbf{\theta})}{\sqrt{C_{rr}(\mathbf{0}; \mathbf{\theta})} \sqrt{C_{ss}(\mathbf{0}; \mathbf{\theta})}}
$$

where:

$$C_{rs}(\mathbf{h}) = \sum_{j=1}^k \lambda_{rj}\lambda_{sj} \rho(\mathbf{h}; \varphi_j)$$

## Applications

BIPPS is particularly useful for:

- Analyzing spatial point patterns in ecology
- Modeling disease incidence across geographical regions
- Studying stellar or galactic distributions
- Any application involving spatially distributed count data

## Contributing

Contributions to the package are welcome. Please feel free to submit issues or pull requests on GitHub.

## Citation

If you use BIPPS in your research, please cite:

[Citation information]
