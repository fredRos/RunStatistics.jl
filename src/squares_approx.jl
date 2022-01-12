# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

export approx_cumulative, approx_pvalue

# NOTE: in many cases the Float64s here could probably be of a smaller type e.g. Float34 or something. check if it  makes a difference
# NOTE: is this a good way to implement methods with different argument types?


mutable struct IntegrandData
    T_obs::Float64
    Nl::Int
    Nr::Int
end

#=

mutable struct CubaIntegrandData

    T_obs::Real # check if needs to be Real, or something less suffices. also with IntegrandData
    Nl::Int 
    Nr::Int
    counter::Int
    spline::Any
    acc::Any

end

d = CubaIntegrandData(1.0, 1, 1, 0, 0, 0)

=#

function h(chisq::Float64, N::Int)
    res = 0
    weight = 0.5

    for i = 1:N
        if (i < N)
            weight *= 0.5
        end
        res += weight * pdf(Chisq(i), chisq)
    end

    return res
end


function H(a::Float64, b::Float64, N::Int)

    res = 0
    weight = 0.5

    for i = 1:N
        if (i < N)
            weight *= 0.5
        end

        res += weight * (cdf(Chisq(i), b) - cdf(Chisq(i), a))
    end

    return res
end


function (integrand::IntegrandData)(x::Float64)
    return h(x, integrand.Nl) * H(integrand.T_obs - x, integrand.T_obs, integrand.Nr)  
end


function Delta(T_obs::Real, Nl::Int, Nr::Int, epsrel::Real, epsabs::Real)

    #NOTE:  this doesn't support passing on additional parameters to the integrand function :( is this a problem?
    #       could maybe also work with cubature.jl. is it sensible to only use one package?
    F = IntegrandData(T_obs, Nl, Nr)

    return quadgk(F, 0, T_obs, rtol = epsrel, atol = epsabs, order = 10)
end


function approx_cumulative(T_obs::Real, N::Int, n::Real, epsrel::Real = nothing, epsabs::Real = nothing)

    F = cumulative(T_obs, N)
    Fn1 = (F / (1 + Delta(T_obs, N, N, epsrel, epsabs)[1]))^(n - 1)
    return F * Fn1
end


function approx_pvalue(T_obs::Real, N::Int, n::Real, epsrel::Real, epsabs::Real)

    return 1 - approx_cumulative(T_obs, N, n, epsrel, epsabs)
end


#=

function (c_int::CubaIntegrandData)(uv::AbstractArray) # analog to cubature_integrand in fred's code. investigate the storage pointers he uses 
    c_int.counter += 1 

    u = uv[1]
    v = uv[2]
    x = c_int.T_obs * u * (1 - v) + c_int.T_obs * (1 - u)
    y = c_int.T_obs * u * v + c_int.T_obs * (1 - u)
    jac = c_int.T_obs * c_int.T_obs * u

    fval = jac * h(x, c_int.Nl) * h(y, c_int.Nr)
    fval *= (c_int.spline && c_int.acc) ? spline_eval(c_int.spline, x + y, c_int.acc) : cumulative(c_int.T_obs, c_int.Nl + c_int.Nr)

    return fval
end 


function full_correction(T_obs::Real, Nl::Int, Nr::Int, epsrel::Real, epsabs::Real, ninterp::Int) # what is this even exactly supposed to do??


    #wonky at best

    data = CubaIntegrandData(T_obs, Nl, Nr, 0)

    #find out how to emulate the stuff with null pointers in the c++ code


    if (ninterp >= 2)
        acc = interp_accel_alloc()
        spline = spline_alloc(interp_linear, ninterp)

        x = zeros(ninterp)
        y = zeros(ninterp)

        for i in range(1, ninterp)
            x[i] = T_obs + i * T_obs / (ninterp - 1)
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
