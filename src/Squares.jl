# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).


using ArgCheck
using QNaNs


log_factorial = [] ::AbstractArray{Float64}



function CacheFactorials(N::Int)::Unsigned
    if N < length(log_factorial)
        return length(log_factorial)
    end

    sizehint!(log_factorial, N)

    if isempty(log_factorial)
        push!(log_factorial, 0.0) # make sure this is an ok substitute for 0L in c++
    end 

    for i in length(log_factorial) : N
        push!(log_factorial, last(log_factorial) + log(Float64(i)))
    end

    return length(log_factorial)
end


function CacheChi2(Tobs::Float64, N::Int)::AbstractArray{Float64}
    
    @argcheck 0 < N

    res = zeros(N + 1)
    res[1] = qnan(0) # ? does this fulfill the intended purpose?

    #create team of threads somehow?

    for i in range(2, N)
        res[i] = log(cquantile(Chisq(Tobs), i))
    end

    return res
end


function cumulative(Tobs::Float64, N:Int)::Float64
    
    CacheFactorials(N)

    log_cumulative = CacheChi2(Tons, N) # find analogon for 'auto' from c++ in julia 

    poch = 0.0::Float64 

    logpow2N1 = (N <= 63) ? log((1 << N) - 1) : N * log(2) # find analogon for 'const' from c++ in julia; also does the bitshift work properly?

    p = 0.0::Float64

    #again create team of threads

    for r in range(1, N)
        Mmax = min(r, N - r + 1)
        poch = 0

        for M in range(1, Mmax)
            poch += log(Float64(N - r + 2 - M))
            scale = poch - logpow2N1::Float64
            ppi = 0.0::Float64

            g(r, M) = e

        end
    end

    return p
end


function pvalue(Tobs::Float64, N:Int)::Float64
    return 1 - cumulative(Tobs, N)    
end

