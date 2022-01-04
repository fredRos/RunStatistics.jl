# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).
__precompile__()

"""
This package implements the cumulative distribution function of the weighted-runs statistic originally defined in

    Frederik Beaujean and Allen Caldwell. A Test Statistic for Weighted Runs. Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. doi:10.1016/j.jspi.2011.04.022 https://arxiv.org/abs/1005.3233
    The authors further derived an approximation to be able to compute the cumulative also for large numbers of observations in
    Frederik Beaujean and Allen Caldwell. Is the bump significant? An axion-search example https://arxiv.org/abs/1710.06642 
    Where they renamed the weighted-runs statistic to the SQUARES statistic.
    This code is based on the original implementation by Frederik Beaujean in c++ and mathematica.

"""
module RunStatistics

using ArgCheck
using Distributions
using QuadGK


include("tobs.jl")
include("partitions.jl")
include("squares.jl")
include("squares_approx.jl")

end # module
