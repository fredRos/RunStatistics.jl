# RunStatistics.jl
[![Documentation for stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://bat.github.io/RunStatistics.jl/stable)
[![Documentation for development version](https://img.shields.io/badge/docs-dev-blue.svg)](https://bat.github.io/RunStatistics.jl/dev)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://github.com/bat/RunStatistics.jl/workflows/CI/badge.svg?branch=main)](https://github.com/bat/RunStatistics.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/bat/RunStatistics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bat/RunStatistics.jl)

This package implements the cumulative distribution function of the weighted-runs statistic originally defined in 

Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. [doi:10.1016/j.jspi.2011.04.022](http://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](http://arxiv.org/abs/1005.3233)

The authors further derived an approximation to be able to compute the cumulative also for large numbers of observations in

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](http://arxiv.org/abs/1710.06642)

where they renamed the weighted-runs statistic to the SQUARES statistic.

This code is based on the [original implementation](https://github.com/fredRos/runs) by Frederik Beaujean in c++ and mathematica.


## Installation
-------

To install RunStatistics.jl, start Julia and run 

```Julia
julia> using Pkg
julia> pkg"add RunStatistics"
```

## SQUARES statistic
-------

When calculating the p value P(T >= Tobs | N) for a sequence of `N` observations, first the value of the SQUARES test statistic `Tobs` needs to be calculated. It denotes the largest $\chi^2$ of any run of consecutive successes (above expectation) in a sequence of `N` independent trials with Gaussian uncertainty

For the Squares statistic to be calculable, the observed data must satisfy following conditions:

> 1. All observations {$\sf{X_{i}}$} are independent. 
> 2. Each observation is normally distributed, Xi ∼ N( $\sf{µ_{i}}$
, $\sf{σ^2_{i}}$ ). 
> 3. Mean $\sf{µ_{i}}$ and variance $\sf{σ^2_{i}}$ are known.

Calculating `Tobs` for the observed data {$\sf{X_{i}}$} is done with the `tobs()` function:

```Julia
Tobs = tobs(X::AbstractArray, μ::Real, σ2::Real)
```

Where `X` is a vector containing the observations, and `μ` and `σ2` are their mean and variance.

If the ovservations don't all have the same mean and variance, use `tobs_ind()`:

```Julia
Tobs = tobs_ind(X::AbstractArray, μ::AbstractArray, σ2::AbstractArray)
```

with the i-th elements of `μ` and `σ2` are the mean and variance of the i-th element of `X`.

The cumulative distribution `P(T < Tobs | N)` and the p value `P(T >= Tobs | N)` are calculated by:

```Julia 
julia> cumulative(Tobs::Float64, N::Int)

julia> pvalue(Tobs::Float64, N::Int)
```

### Approximation for large N
-------

For large `N`, the number of therms in the exact expression scales like `exp(N^1/2)/N` and quickly grows too large. An approximate formula is implemented here for `n*N`, where the cumulative for example `N = 100` is computed exactly and `n` may be 1 or >> `N` and need not even be an integer. 

The cumulative distribution P(T < Tobs | n\*N) and the p value P(T >= Tobs | n\*N) are approximated by:

```Julia 
julia> approx_cumulative(Tobs::Float64, N::Int, n::Float64, epsrel::Float64, epsabs::Float64)

julia> approx_pvalue(Tobs::Float64, N::Int, n::Float64, epsrel::Float64, epsabs::Float64)
```

The approximation involves a 1D numerical integration whose relative and absolute target precision are `epsrel` and `epsabs`. 
## Documentation

* [Documentation for stable version](https://bat.github.io/RunStatistics.jl/stable)
* [Documentation for development version](https://bat.github.io/RunStatistics.jl/dev)


## Citing RunStatistics.jl

When using RunStatistics.jl for research, teaching or similar, please cite the original authors' work:

```
@article{beaujean2011test,
title={A test statistic for weighted runs},
author={Beaujean, Frederik and Caldwell, Allen},
journal={Journal of Statistical Planning and Inference},
volume={141},
number={11},
pages={3437--3446},
year={2011},
publisher={Elsevier}
}

@article{Beaujean:2017eyq,
  author         = "Beaujean, Frederik and Caldwell, Allen and Reimann, Olaf",
  title          = "{Is the bump significant? An axion-search example}",
  year           = "2017",
  eprint         = "1710.06642",
  archivePrefix  = "arXiv",
  primaryClass   = "hep-ex",
  SLACcitation   = "%%CITATION = ARXIV:1710.06642;%%"
}
```