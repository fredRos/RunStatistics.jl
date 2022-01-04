# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).
__precompile__()

module RunStatistics

export pvalue, cumulative

include("tobs.jl")
include("partitions.jl")
include("squares.jl")
include("squares_approx.jl")
end # module
