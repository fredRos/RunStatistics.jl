# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).
__precompile__()

"""
This package implements the evaluation of the cumulative distribution function of the `Squares test statistic` originally defined in

Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs.*

https://arxiv.org/abs/1005.3233

The authors further derived an approximation to be able to compute the cumulative also for large numbers of observations in

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* 

https://arxiv.org/abs/1710.06642 

Where they renamed the weighted-runs statistic to the `*SQUARES statistic*`.
This code is based on the original implementation by Frederik Beaujean in c++ and mathematica:

https://github.com/fredRos/runs
"""
module RunStatistics


using Distributions
using QuadGK
using ArgCheck
using Distributions

include("t_obs.jl")
include("partitions.jl")
include("squares.jl")
include("squares_approx.jl")

end # module
