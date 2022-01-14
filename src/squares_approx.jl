# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

export approx_cumulative, approx_pvalue

"""
    IntegrandData

Represent the parameters needed for the 1D numerical integration performed in `Delta()`.

`T_obs` is the value for the Squares statistic observed in the data, `Nl` the left-hand lenght and 
`Nr` the right-hand lenght of a boundary spanning run, as defined in section II.A. in

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
"""
mutable struct IntegrandData
    T_obs::Float64
    Nl::Int
    Nr::Int
end

"""
    h(chisq::Real, N::Integer)

Compute the probability density h(χ2 | Nr) for the right-hand side of a boundary spanning run to be above expectation; 
as explained in section II.A. in the paper below.

Caclulate it as the sum of probability densities for runs of different length times the χ2 probability for that number of degrees of freedom.

Implements the term defined in equation (8) in 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
"""
function h(chisq::Real, N::Integer)
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


"""
    H(a::Real, b::Real, N::Integer)

Compute the cumulative of `h()` as defined in section II.A. in

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
"""
function H(a::Real, b::Real, N::Integer)

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

"""
    (integrand::IntegrandData)(x::Real)

Compute the integrand in the Δ(T_obs | N_l, N_r) term defined in equation (13) in 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
"""
function (integrand::IntegrandData)(x::Real)
    return h(x, integrand.Nl) * H(integrand.T_obs - x, integrand.T_obs, integrand.Nr)  
end

 
"""
    Delta(T_obs::Real, Nl::Integer, Nr::Integer, epsrel::Real, epsabs::Real)

Compute the Δ(T_obs | N_l, N_r) term defined in equation (13) in 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642

The calculation involves a 1D numerical integration using the `quadgk()` function with the 
relative and absolute target precision `epsrel` and `epsabs`. If not specified, the default values of `quadgk()` are used.
See https://juliamath.github.io/QuadGK.jl/stable/ for documentation.
"""
function Delta(T_obs::Real, Nl::Integer, Nr::Integer, epsrel::Real, epsabs::Real)

    F = IntegrandData(T_obs, Nl, Nr)
    return quadgk(F, 0, T_obs, rtol = epsrel, atol = epsabs, order = 10)
end

"""
approx_cumulative(T_obs::Real, N::Integer, n::Real, [epsrel::Real, epsabs::Real])

Compute an approximation of P(T < `T_obs` | `n * N`), the value of the cumulative distribution function for the Squares test statistic at `T_obs`, 
the value of the Squares statistic observed in the data. 
The total number of datapoints is `n * N`

The approximation involves a 1D numerical integration whose relative and absolute target precision are `epsrel` and `epsabs`. 
These are optional arguments here and if left empty, the default values of the `quadgk()` function are used. See https://juliamath.github.io/QuadGK.jl/stable/

This function implements equation (17) from:

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
  
"""
function approx_cumulative(T_obs::Real, N::Integer, n::Real, epsrel::Real = nothing, epsabs::Real = nothing)

    F = cumulative(T_obs, N)
    Fn1 = (F / (1 + Delta(T_obs, N, N, epsrel, epsabs)[1]))^(n - 1)
    return F * Fn1
end


"""
    approx_pvalue(T_obs::Real, N::Integer, n::Real, [epsrel::Real, epsabs::Real])

Compute an approximation of P(T >= `T_obs` | `n * N`), the p value for the Squares test statistic T being larger or equal to `T_obs`, 
the value of the Squares statistic observed in the data. 
The total number of datapoints is `n * N`.

The approximation involves a 1D numerical integration whose relative and absolute target precision are `epsrel` and `epsabs`. 
These are optional arguments here and if left empty, the default values of the `quadgk()` function are used. See https://juliamath.github.io/QuadGK.jl/stable/

Via `approx_cumulative` this function implements equation (17) from:

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
  
"""
function approx_pvalue(T_obs::Real, N::Integer, n::Real, epsrel::Real, epsabs::Real)

    return 1 - approx_cumulative(T_obs, N, n, epsrel, epsabs)
end



