# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


"""
    tobs(X::AbstractArray, μ::Float64, σ2::Float64)

Cumpute the value of the squares test statistic `Tobs` i.e. the largest \\chi^2 of any run of consecutive successes (above expectation)
in a sequence of `N` independent trials with Gaussian uncertainty. With μ and σ2 being the expectation and variance of the observations.

For the Squares statistic to be calculable, the observed data must satisfy following conditions:

        All observations {X_i} are independent.
        Each observation is normally distributed, Xi ∼ N( µ_i, σ^2_i ).
        Mean µ_i and variance σ^2_i are known.
"""
function tobs(X::AbstractArray, μ::Float64, σ2::Float64)
    χ2 = []
    χi = 0
    for i in range(1, length(X); step = 1)
        X[i] > μ ?  χi += (X[i] - μ)^2 / σ2 : (append!(χ2, χi); χi = 0)
    end 
    
    return maximum(χ2)
end

function tobs(X::AbstractArray, μ::AbstractArray, σ2::AbstractArray)
    χ2 = []
    χi = 0

    for i in range(1, length(X))
        X[i] > μ[i] ?  χi += (X[i] - μ[i])^2 / σ2[i] : (append!(χ2, χi); χi = 0)
    end 

    return maximum(χ2)
end

