# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


export pvalue, cumulative

const log_factorial = Vector{Float64}(undef, 0)

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
    cumulative(T_obs::Real, N::Integer)

Compute P(T < T_obs | N), the value of the cumulative distribution of the Squares statistic `T` at the 
value `T_obs` observed in the data.

`T_obs` is the value of the test statistic for the observed data set; i.e., the largest `χ^2` of 
any run of consecutive observed values above the expectation in a sequence of `N` independent trials 
with Gaussian uncertainty. 

`N` is the total number of data points.

The calculation implements equations (16) and (17) from 

Frederik Beaujean and Allen Caldwell. *A Test Statistic for
Weighted Runs.* Journal of Statistical Planning and Inference 141,
no. 11 (November 2011): 3437–46. doi:10.1016/j.jspi.2011.04.022

https://arxiv.org/abs/1005.3233.

"""
function cumulative(T_obs::Real, N::Integer)

    T = promote_type(typeof(T_obs), typeof(N))

    if N >= length(log_factorial)

        sizehint!(log_factorial, N)

        if isempty(log_factorial)
            push!(log_factorial, 0.0)
        end

        for i = length(log_factorial):N
            push!(log_factorial, last(log_factorial) + log(i))
        end

    end


    log_cumulative = cachechi2(T_obs, N)

    logpow2N1 = (N <= 63) ? log((1 << N) - 1) : N * log(2)

    p = zero(T)

    # p = Threads.Atomic{float(T)}(0)
    # Ps = zeros(Real, numthreads)

    #TODO: think about creating a Partition() object for each thread and then initiate it within the thread
    #Threads.@threads 
    for r = 1:N

        Mmax = min(r, N - r + 1)
        poch = 0.0

        for M = 1:Mmax

            poch += log(N - r + 2 - M)
            scale = poch - logpow2N1
            ppi = 0.0

            g = init_partition(r, M)
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

            #Ps[threadid()] += exp(scale + log(ppi))
            p += exp(scale + log(ppi))
        end
    end
    #p = sum(Ps)

    @assert p < 1
    return p
end

"""
    pvalue(T_obs::Real, N::Integer)

Compute P(T >= `T_obs` | `N`), the p value for the Squares test statistic `T` being larger or equal to `T_obs`, 
the value of the Squares statistic observed in the data. 

`N` is the total number of data points.

The Squares statistic `T` denotes the largest `χ^2` of any run of consecutive successes (above expectation) in a 
sequence of `N` independent trials with Gaussian uncertainty.

`T_obs` is the value of the test statistic for the observed data, `N` is the total number of data points.

Via `cumulative()` this function inmplements equations (16) and (17) from

Frederik Beaujean and Allen Caldwell. *A Test Statistic for
Weighted Runs.* Journal of Statistical Planning and Inference 141,
no. 11 (November 2011): 3437–46. doi:10.1016/j.jspi.2011.04.022

https://arxiv.org/abs/1005.3233.
  
"""
function pvalue(T_obs::Real, N::Integer)

    return 1 - cumulative(T_obs, N)

end

