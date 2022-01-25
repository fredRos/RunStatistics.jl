# RunStatistics.jl
[![Documentation for stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://bat.github.io/RunStatistics.jl/stable)
[![Documentation for development version](https://img.shields.io/badge/docs-dev-blue.svg)](https://bat.github.io/RunStatistics.jl/dev)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://github.com/bat/RunStatistics.jl/workflows/CI/badge.svg?branch=main)](https://github.com/bat/RunStatistics.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/bat/RunStatistics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bat/RunStatistics.jl)

This package implements the evaluation of the cumulative distribution function `P(T < T_obs | N)` of the weighted-runs statistic  originally defined in 

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
and 
```Julia
julia> using RunStatistics
```
to use the functions the package provides.

For an explanation on the theory behind the package and how to use it for your data, see the [Documentation for stable version](https://bat.github.io/RunStatistics.jl/stable).
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