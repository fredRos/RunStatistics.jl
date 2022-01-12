# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


export pvalue, cumulative

const log_factorial = Vector{Float64}(undef, 0)

function cachefactorials(N::Integer)

    if N < length(log_factorial)
        return length(log_factorial)
    end

    sizehint!(log_factorial, N)

    if isempty(log_factorial)
        push!(log_factorial, 0.0)
    end

    for i = length(log_factorial):N
        push!(log_factorial, last(log_factorial) + log(i))
    end

    return length(log_factorial)
end


function cachechi2(T_obs::Real, N::Int)

    @argcheck 0 < N

    res = zeros(N + 1)
    res[1] = NaN

    for i = firstindex(res)+1:lastindex(res)
        
        res[i] = logcdf(Chisq(i - 1), T_obs)
    end

    return res
end

"""
    cumulative(T_obs::Float64, N::Int)

Compute the cumulative distribution of the weighted-runs, or SQUARES, statistic `T`.

`T_obs` is the value of the test statistic for the observed data set; i.e., the largest chi^2 of 
any run of consecutive observed values above the expectation. `N` is the total number of data points.

The calculation implements Eqns. (16) and (17) from Frederik Beaujean and Allen Caldwell. “A Test Statistic for
Weighted Runs.” Journal of Statistical Planning and Inference 141,
no. 11 (November 2011): 3437–46. doi:10.1016/j.jspi.2011.04.022
http://arxiv.org/abs/1005.3233.

"""
function cumulative(T_obs::Real, N::Int)

    cachefactorials(N)

    log_cumulative = cachechi2(T_obs, N)

    logpow2N1 = (N <= 63) ? log((1 << N) - 1) : N * log(2)

    # p = Threads.Atomic{float(T)}(0)
    Ps = zeros(Real, numthreads)

    #TODO: think about creating a Partition() object for each thread and then initiate it within the thread
    Threads.@threads for r = 1:N

        Mmax = min(r, N - r + 1)
        poch = 0.0

        for M = 1:Mmax

            poch += log(N - r + 2 - M)
            scale = poch - logpow2N1
            ppi = 0.0

            g = Partition(r, M)
            n = g.c
            y = g.y

            done = false

            while !done

                h = g.h
                ppartition = 0.0

                for l = 2:(h + 1)
                    ppartition += n[l] * log_cumulative[y[l]+1] - log_factorial[n[l]+1]
                end

                ppi += exp(ppartition)

                done = next_partition!(g)
            end

            Ps[threadid()] += exp(scale + log(ppi))
        end
    end
    p = sum(Ps)

    @assert p < 1
    return p
end

"""
    pvalue(T_obs::Float64, N::Int)

    Compute the p value P(T >= `T_obs` | `N`) with `T_obs` being the value of the squares test statistic,
    i.e. the larges chi^2 of any run of consecutive successes (above expectation) in a sequence of `N` 
    independent trials with Gaussian uncertainty.
"""
function pvalue(T_obs::Real, N::Int)

    return 1 - cumulative(T_obs, N)

end

runs