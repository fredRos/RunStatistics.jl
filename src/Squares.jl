# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


using ArgCheck
using Distributions


# include("partitions.jl")

log_factorial = []

function CacheFactorials(N)
    if N < length(log_factorial)
        return length(log_factorial)
    end

    sizehint!(log_factorial, N)

    if isempty(log_factorial)
        push!(log_factorial, 0.0)
    end 

    for i in range(length(log_factorial), N)
        push!(log_factorial, last(log_factorial) + log(i))
    end

    return length(log_factorial)
end

function CacheChi2(Tobs::Float64, N::Int)
    
    @argcheck 0 < N

    res = zeros(N + 1)
    res[1] = NaN # ? does this fulfill the intended purpose?

    Threads.@threads for i in range(2, N + 1)
        res[i] = logcdf(Chisq(i - 1), Tobs)
    end

    return res
end


"""
    cumulative(Tobs::Float64, N::Int)

Compute the cumulative distribution of the weighted-runs, or SQUARES, statistic `T`.

`Tobs` is the value of the test statistic for the observed data set; i.e., the largest \chi^2 of 
any run of consecutive observed values above the expectation. `N` is the total number of data points.

The calculation implements Eqns. (16) and (17) from

Frederik Beaujean and Allen Caldwell. “A Test Statistic for
Weighted Runs.” Journal of Statistical Planning and Inference 141,
no. 11 (November 2011): 3437–46. doi:10.1016/j.jspi.2011.04.022
http://arxiv.org/abs/1005.3233.

"""
function cumulative(Tobs::Float64, N::Int)
    
    CacheFactorials(N)

    log_cumulative = CacheChi2(Tobs, N)

    poch = 0.0

    logpow2N1 = (N <= 63) ? log((1 << N) - 1) : N * log(2) 

    p = 0.0

    Threads.@threads for r in range(1, N)
        Mmax = min(r, N - r + 1)
        poch = 0

        for M in range(1, Mmax)
            poch += log(N - r + 2 - M)
            scale = poch - logpow2N1
            ppi = 0.0

            g = partition(r, M)
            n = g.c
            y = g.y

            check = true
            while check
                h = g.h

                ppartition = 0
                for l in range(2, h + 1)
                    ppartition += n[l] * log_cumulative[y[l] + 1] - log_factorial[n[l] + 1]
                end 
                ppi += exp(ppartition)
                
                check = !final_partition(g)
                iterate!(g)
            end

            p += exp(scale + log(ppi))
        end
    end

    @assert p < 1
    return p
end

# println(cumulative(3.4, 10), " ", cumulative(6.5, 20), " ", cumulative(9.0, 6), " ", cumulative(5.4, 4), " ", cumulative(3.9, 11)) 
"""
    pvalue(Tobs::Float64, N::Int)

Compute the p value P(T >= `Tobs` | `N`) with `Tobs` being the value of the squares test statistic,
i.e. the larges chi^2 of any run of consecutive successes (above expectation) in a sequence of `N` 
independent trials with Gaussian uncertainty.
"""
function pvalue(Tobs::Float64, N::Int)
    return 1 - cumulative(Tobs, N)    
end

