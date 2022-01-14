# RunStatistics.jl
[![Documentation for stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://bat.github.io/RunStatistics.jl/stable)
[![Documentation for development version](https://img.shields.io/badge/docs-dev-blue.svg)](https://bat.github.io/RunStatistics.jl/dev)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://github.com/bat/RunStatistics.jl/workflows/CI/badge.svg?branch=main)](https://github.com/bat/RunStatistics.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/bat/RunStatistics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bat/RunStatistics.jl)

This package implements the evaluation of the cumulative distribution function P(T < T_obs | N) of the weighted-runs statistic  originally defined in 

Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. [doi:10.1016/j.jspi.2011.04.022](http://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](http://arxiv.org/abs/1005.3233)

The authors further derived an approximation to be able to compute the cumulative also for large numbers of observations in

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](http://arxiv.org/abs/1710.06642)

where they renamed the weighted-runs statistic to the SQUARES statistic.

This code is based on the [original implementation](https://github.com/fredRos/runs) by Frederik Beaujean in c++ and mathematica.


## Installation

To install `RunStatistics.jl`, start Julia and run 

```Julia
julia> using Pkg
julia> pkg"add RunStatistics"
```

## SQUARES statistic
To use the RunStatistics.jl package after installing it, do:

```Julia
julia> using RunStatistics
```
When calculating the p value P(T >= T_obs | N) for a sequence of `N` independent observations with gaussian uncertainty, first the observed value of the `Squares test statistic` `T_obs` needs to be calculated. 

`T_obs` denotes the largest `χ^2` of any run of consecutive successes (above expectation) in this sequence of observations.

For the Squares statistic to be calculable, the observed data must satisfy following conditions:

> 1. All observations {X_i} are independent. 
> 2. Each observation is normally distributed, X_i ∼ N( µ_i, σ^2_i ). 
> 3. Mean µ_i and variance σ^2_i are known.

Calculating `T_obs` for the observed data X_i is done with the `t_obs()` function:

```Julia
T_obs = t_obs(X::AbstractArray, μ::Real, σ2::Real)
```

Where `X` is a vector containing the observations, and `μ` and `σ2` are their mean and variance.

If the obvservations don't all have the same mean and variance, use:

```Julia
T_obs = T_obs(X::AbstractArray, μ::AbstractArray, σ2::AbstractArray)
```

with the i-th elements of `μ` and `σ2` being the mean and variance of the i-th element of `X`.

The cumulative distribution `P(T < T_obs | N)` and the p value `P(T >= T_obs | N)` at the observed value `T_obs` are calculated with:

```Julia 
julia> cumulative(T_obs::Float64, N::Int)

julia> pvalue(T_obs::Float64, N::Int)
```

### Approximation for large N

For large `N`, the number of terms in the exact expression scales like `exp(N^1/2)/N` and quickly grows too large. An approximate formula is implemented here for cases when the total number of data points is `n*N`.

So for example, if 50 000 data points were collected, choose `N = 50` and `n = 1000`.

For this example the approximation the cumulative for `N = 50` is computed exactly and then further processed to obtain the desired approximation, via equation (17) from 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](http://arxiv.org/abs/1710.06642)

As a rule of thumb, `N` should not exceed `100`. 

The cumulative distribution P(T < T_obs | n\*N) and the p value P(T >= T_obs | n\*N) are approximated by:

```Julia 
julia> approx_cumulative(T_obs::Float64, N::Int, n::Float64, [epsrel::Float64, epsabs::Float64])

julia> approx_pvalue(T_obs::Float64, N::Int, n::Float64, [epsrel::Float64, epsabs::Float64])
```

The approximation involves a 1D numerical integration whose relative and absolute target precision are `epsrel` and `epsabs`; these are `optional arguments` in the above functions. 
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