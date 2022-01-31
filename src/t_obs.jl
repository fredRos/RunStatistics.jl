# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

"""
    t_obs(X::AbstractArray, μ::Real, σ2::Real)

Compute the value of the *Squares test statistic* `T_obs` i.e. the largest `χ2` of any run of consecutive successes (above expectation)
in a sequence of `N` independent trials with Gaussian uncertainty. `μ` and `σ2` are the expectation and variance of the observations.

Find the location(s) of the run(s) that produces `T_obs`.

Returns a tuple containing `T_obs` and one or more arrays containing the indices of the runs that produce `T_obs`.

For the Squares statistic to be calculable, the observed data must satisfy following conditions:

        All observations {X_i} are independent. 
        Each observation is normally distributed, X_i ∼ N( µ_i, σ^2_i ).
        Mean µ_i and variance σ^2_i are known.

In case the observations {X_i} have individual expectations and variances, use:

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

    locs_cache = Array{Array}(undef,0)
    locs_run = Array{Int64}(undef,0)
    locs = Array{Union{Array, Int64}}(undef,0)

    run_at_end = false

    @inbounds for i in eachindex(X) 

        if X[i] > μ
            χi += (X[i] - μ)^2 / σ2
            append!(locs_run, i)
            run_at_end = true

        else
            χi > 0 && append!(χ2, χi) 
            isempty(locs_run) || append!(locs_cache, [locs_run])
            χi = zero(T)
            locs_run = Array{Int64}(undef,0)
            run_at_end = false
        end 

    end

    run_at_end && (append!(χ2, χi); append!(locs_cache, [locs_run]))

    locs_ids = findall(x -> x == maximum(χ2), χ2)

    if length(locs_ids) == 1 
        return maximum(χ2), locs_cache[locs_ids[1]]

    else
        for i in locs_ids
            append!(locs, [locs_cache[i]])
        end
    end 

    return maximum(χ2), locs
end

function t_obs(X::AbstractArray, μ::AbstractArray, σ2::AbstractArray)

    @argcheck axes(X) == axes(μ) == axes(σ2)
    T = float(promote_type(eltype(X), eltype(μ), eltype(σ2)))

    χ2 = T[]
    χi = zero(T)

    locs_cache = Array{Array}(undef,0)
    locs_run = Array{Int64}(undef,0)
    locs = Array{Union{Array, Int64}}(undef,0)

    run_at_end = false

    @inbounds for i in eachindex(X)

        if  X[i] > μ[i]

            χi += (X[i] - μ[i])^2 / σ2[i]
            append!(locs_run, i)
            run_at_end = true
        else
            χi > 0 && append!(χ2, χi)
            isempty(locs_run) || append!(locs_cache, [locs_run])
            χi = 0
            locs_run = Integer[]
            run_at_end = false
        end
    end

    run_at_end && (append!(χ2, χi); append!(locs_cache, [locs_run]))

    locs_ids = findall(x -> x == maximum(χ2), χ2)

    if length(locs_ids) == 1 

        return maximum(χ2), locs_cache[locs_ids[1]]

    else
        for i in locs_ids

            append!(locs, [locs_cache[i]])
        end
    end 

    return maximum(χ2), locs
end