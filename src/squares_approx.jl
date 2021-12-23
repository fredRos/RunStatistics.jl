# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using Distributions
#using GSL
#using HCubature
#using StaticArrays
using QuadGK

#include("squares.jl")

# in many cases the Float64s here could probably be of a smaller type e.g. Float34 or something. check if it  makes a difference

# fred often works with pointers to storage locations of values of variables instead of variables themselves. Check how to do that properly here

mutable struct IntegrandData 
    Tobs::Float64
    Nl::Int
    Nr::Int
end



mutable struct CubaIntegrandData
    Tobs::Float64 # check if needs to be Float64, or something less suffices. also with IntegrandData
    Nl::Int 
    Nr::Int
    counter::Int
    spline::Any
    acc::Any
end

d = CubaIntegrandData(1.0, 1, 1, 0, 0, 0)


function h(chisq::Float64, N::Int)

    res = 0
    weight = 0.5

    for i in range(1, N)
        if (i < N)
            weight *= 0.5
        end
        res += weight * pdf(Chisq(i), chisq) 
    end

    # println("resh: ", res)
    return res
end


function H(a::Float64, b::Float64, N::Int)
    res = 0
    weight = 0.5
    
    for i in range(1, N)
        if (i < N)
            weight *= 0.5
        end
        res += weight * (cdf(Chisq(i), b) - cdf(Chisq(i), a))
    end
    # println("resH: ", res)
    return res
end


function (integrand::IntegrandData)(x::Float64)
    # println("reshH: ", h(x, e.Nl) * H(e.Tobs - x, e.Tobs, e.Nr))
    return h(x, integrand.Nl) * H(integrand.Tobs - x, integrand.Tobs, integrand.Nr)
end 


function Delta(Tobs::Float64, Nl::Int, Nr::Int, epsrel::Float64, epsabs::Float64, maxevals=10^3) # thinks about making error args optional
    # this doesn't support passing on additional parameters to the integrand function :( is this a problem?
    # could maybe also work with cubature.jl. is it sensible to only use one package?
    F = IntegrandData(Tobs, Nl, Nr)

    return quadgk(F, 0, Tobs, rtol=epsrel, atol=epsabs, maxevals=maxevals, order=10)  # maybe delete F afterwards?
end

#=

function (c_int::CubaIntegrandData)(uv::AbstractArray) # analog to cubature_integrand in fred's code. investigate the storage pointers he uses 
    c_int.counter += 1 

    u = uv[1]
    v = uv[2]
    x = c_int.Tobs * u * (1 - v) + c_int.Tobs * (1 - u)
    y = c_int.Tobs * u * v + c_int.Tobs * (1 - u)
    jac = c_int.Tobs * c_int.Tobs * u

    fval = jac * h(x, c_int.Nl) * h(y, c_int.Nr)
    fval *= (c_int.spline && c_int.acc) ? spline_eval(c_int.spline, x + y, c_int.acc) : cumulative(c_int.Tobs, c_int.Nl + c_int.Nr)

    return fval
end 


function full_correction(Tobs::Float64, Nl::Int, Nr::Int, epsrel::Float64, epsabs::Float64, ninterp::Int) # what is this even exactly supposed to do??
    

    #wonky at best

    data = CubaIntegrandData(Tobs, Nl, Nr, 0)
    
    #find out how to emulate the stuff with null pointers in the c++ code


    if (ninterp >= 2)
        acc = interp_accel_alloc()
        spline = spline_alloc(interp_linear, ninterp)

        x = zeros(ninterp)
        y = zeros(ninterp)

        for i in range(1, ninterp)
            x[i] = Tobs + i * Tobs / (ninterp - 1)
            y[i] = cumulative(x[i], Nl + Nr)
        end
    
        spline_init(spline, x, y, ninterp)
    end

    data.acc = acc
    data.spline = spline

    #fdim = 1
    #dim = 2
    uvmin = [0, 0]
    uvmax = [1, 1]
    maxEval = 10000

    #=
    if hcubature

    throw runtime_error

    end
    =#

    res = hcubature(cubature_integrand, uvmin, uvmax, norm = norm, epsrel, epsabs, maxEval)

    # consider wheter using StaticArrays for performance is necessary; see HCubature.jl package wiki 

    spline_free(spline)
    interp_accel_free(acc)
    
    x = []
    y = [] # decide if this is enough, or x,y should somehow be deleted as in the c++ code
    
    return res

end
=#


# find sensible default values for epsrel and epsabs.

function approx_cumulative(Tobs::Float64, N::Int, n::Float64, epsrel::Float64, epsabs::Float64)

    F = cumulative(Tobs, N)
    Fn1 = (F / (1 + Delta(Tobs, N, N, epsrel, epsabs)[1])) ^ (n - 1)
    return F * Fn1
end

function approx_pvalue(Tobs::Float64, N::Int, n::Float64, epsrel::Float64, epsabs::Float64)

    return 1 - approx_cumulative(Tobs, N, n, epsrel, epsabs)
end

#println(approx_pvalue(30.0, 50, 100.0, 0.1, 0.001))