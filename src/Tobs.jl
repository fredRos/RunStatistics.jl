# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

# Calculates the Test statistic T = max_j χ^2_{run,j} for a sequence of observations that fulfill the assumptions:
# 1 All observations are independent
# 2 Each observation is normally distributed
# 3 Mean μ and variance σ^2 are known


#Hello world
#Think about implementing a special case with μ and σ as vectors, so each obsevation X[i] can have an individual normal distribution

function tobs(X::AbstractArray, μ::Real, σ2::Real)
    χ2 = []
    χi = 0

    for i in range(1, length(X))
        X[i] > μ ?  χi += (X[i] - μ)^2 / σ2 : (append!(χ2, χi); χi = 0)
    end 
    
    return maximum(χ2)
end

function tobs_ind(X::AbstractArray, μ::AbstractArray, σ2::AbstractArray)
    χ2 = []
    χi = 0

    for i in range(1, length(X))
        X[i] > μ[i] ?  χi += (X[i] - μ[i])^2 / σ2[i] : (append!(χ2, χi); χi = 0)
    end 

    return maximum(χ2)
end


#=

example of docstring? 

"""
    RunStatistics.hello_world()

Prints "Hello, World!" and returns 42.

```jldoctest
using RunStatistics

RunStatistics.hello_world()

# output

Hello, World!
42
```
"""

=#