# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

"""
    t_obs(X::AbstractArray, μ::Real, σ2::Real)

Compute the value of the *Squares test statistic* `T_obs` i.e. the largest `χ2` of any run of consecutive successes (above expectation)
in a sequence of `N` independent trials with Gaussian uncertainty. `μ` and `σ2` are the expectation and variance of the observations.

For the Squares statistic to be calculable, the observed data must satisfy following conditions:

        All observations {X_i} are independent. 
        Each observation is normally distributed, X_i ∼ N( µ_i, σ^2_i ).
        Mean µ_i and variance σ^2_i are known.

In case the observations {X_i} have individiual expectations and variances, use:

    t_obs(X::AbstractArray, μ::AbstractArray, σ2::AbstractArray)

With `μ[i]` and `σ2[i]` being the mean and variance of the i-th element of `X`. 

See:

Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs.*
Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. 

https://www.sciencedirect.com/science/article/abs/pii/S0378375811001935?via%3Dihub

https://arxiv.org/abs/1005.3233
"""
function t_obs(X::AbstractArray, μ::Real, σ2::Real)

    T = float(promote_type(eltype(X), eltype(μ), eltype(σ2)))
    χ2 = T[]
    χi = zero(T)

    @inbounds for i in eachindex(X)
        X[i] > μ ? χi += (X[i] - μ)^2 / σ2 : (append!(χ2, χi); χi = 0)
    end

    return maximum(χ2)
end

function t_obs(X::AbstractArray, μ::AbstractArray, σ2::AbstractArray)

    @argcheck axes(X) == axes(μ) == axes(σ2)
    T = float(promote_type(eltype(X), eltype(μ), eltype(σ2)))

    χ2 = T[]
    χi = zero(T)

    @inbounds for i in eachindex(X)
        X[i] > μ[i] ? χi += (X[i] - μ[i])^2 / σ2[i] : (append!(χ2, χi); χi = 0)
    end

    return maximum(χ2)
end
