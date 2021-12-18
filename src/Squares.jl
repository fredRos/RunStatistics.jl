# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


using ArgCheck
using Distributions


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
    
    # @argcheck 0 < N

    res = zeros(N + 1)
    res[1] = NaN # ? does this fulfill the intended purpose?

    #create team of threads somehow?

    for i in range(2, N + 1)
        res[i] = logcdf(Chisq(i - 1), Tobs)
    end

    return res
end




function cumulative(Tobs::Float64, N::Int)
    
    CacheFactorials(N)

    log_cumulative = CacheChi2(Tobs, N)

    poch = 0.0

    logpow2N1 = (N <= 63) ? log((1 << N) - 1) : N * log(2) 

    p = 0.0

    #again create team of threads

    for r in range(1, N)
        Mmax = min(r, N - r + 1)
        poch = 0

        for M in range(1, Mmax)
            poch += log(N - r + 2 - M)
            scale = poch - logpow2N1
            ppi = 0.0

            g = partition(r, M)
            n = g.c
            y = g.y

            check = true # to assure that the initial partition is also counted 
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

    # assert p<1
    return p
end

# println(cumulative(3.4, 10), " ", cumulative(6.5, 20), " ", cumulative(9.0, 6), " ", cumulative(5.4, 4), " ", cumulative(3.9, 11)) 

function pvalue(Tobs::Float64, N::Int)
    return 1 - cumulative(Tobs, N)    
end

