# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

export squares_cdf_approx, squares_pvalue_approx

"""
    IntegrandData

Represent the parameters needed for the 1D numerical integration performed in `Delta()`.

`T_obs` is the value for the Squares statistic observed in the data, `Nl` the left-hand length and 
`Nr` the right-hand length of a boundary spanning run, as defined in section II.A. in

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

Compute the probability density h(χ2 | Nr) for the right-hand side of a boundary spanning 
run to be above expectation; 
as explained in section II.A. in the paper below.

Calculate it as the sum of probability densities for runs of different length times the χ2 
    probability for that number of degrees of freedom.

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
    Delta(T_obs::Real, Nl::Integer, Nr::Integer, [epsrel::Nothing, epsabs::::Union{Real, Nothing}])

Compute the Δ(T_obs | N_l, N_r) term defined in equation (13) in 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642

The calculation involves a 1D numerical integration using the `quadgk()` function with the 
relative and absolute target precision `epsrel` and `epsabs`. If not specified, the default 
values of `quadgk()` are used.
See https://juliamath.github.io/QuadGK.jl/stable/ for documentation.
"""
function Delta(T_obs::Real, Nl::Integer, Nr::Integer, epsrel::Nothing = nothing, epsabs::Union{Real, Nothing} = nothing)

    F = IntegrandData(T_obs, Nl, Nr)
    return quadgk(F, 0, T_obs, rtol = epsrel, atol = epsabs, order = 10)
end

"""
squares_cdf_approx(T_obs::Real, L::Integer, [epsp::Real])

Compute an approximation of P(T < `T_obs` | `L = n * N`), the value of the cumulative distribution 
function for the Squares test statistic at `T_obs`, 
the value of the Squares statistic observed in the data. 
The total number of datapoints is `L = n * N`, if not defined otherwise, the function chooses the default values `N = 80` and `n = L / N`.

To specify a certain choice for `N` and `n`, do:

    squares_cdf_approx(T_obs::Real,  Ns::AbstractArray, epsp::Real = 0)

With `Ns` being an array holding `N::Integer` and `n::Real` as its first and second element: Ns = [N, n].

The accuracy's lower bound is `10^(-14)`, a desired accuracy up to this boundary can be specified with the optional `epsp` argument.
See documentation on Accuracy.

This function implements equation (17) from:

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
  
"""
function squares_cdf_approx(T_obs::Real, L::Integer, epsp::Real = 0)
    N = 80
    n = L / N

    @argcheck (epsp == 0 || epsp / n >= 10^(-14)) error("The desired accuracy is too high. See documentation on Accuracy.")

    if epsp != 0

        epsabs = (epsp / n) * 0.1
        F = squares_cdf(T_obs, N)
        Fn1 = (F / (1 + Delta(T_obs, N, N, nothing, epsabs)[1]))^(n - 1)
        return F * Fn1
    end
        
    F = squares_cdf(T_obs, N)
    Fn1 = (F / (1 + Delta(T_obs, N, N)[1]))^(n - 1)
    return F * Fn1
end

function squares_cdf_approx(T_obs::Real, Ns::AbstractArray, epsp::Real = 0)

    N = Integer(Ns[1])
    n = Ns[2]

    @argcheck (epsp == 0 || epsp / n >= 10^(-14)) error("The desired accuracy is too high. See documentation on Accuracy.")

    if epsp != 0

        epsabs = (epsp / n) * 0.1
        F = squares_cdf(T_obs, N)
        Fn1 = (F / (1 + Delta(T_obs, N, N, nothing, epsabs)[1]))^(n - 1)
        return F * Fn1
    end
        
    F = squares_cdf(T_obs, N)
    Fn1 = (F / (1 + Delta(T_obs, N, N)[1]))^(n - 1)
    return F * Fn1
end

"""
    squares_pvalue_approx(T_obs::Real, L::Integer, [epsp::Real])

Compute an approximation of P(T >= `T_obs` | `L`), the p value for the Squares test 
statistic T being larger or equal to `T_obs`, 
the value of the Squares statistic observed in the data. 
The total number of datapoints is `L = n * N`, if not defined otherwise, the function chooses the default values `N = 80` and `n = L / N`.

To specify a certain choice for `N` and `n`, do:

    squares_pvalue_approx(T_obs::Real, Ns::AbstractArray,  [epsp::Real])
 
With `Ns` being an array holding `N::Integer` and `n::Real` as its first and second element: Ns = [N, n]

The accuracy's lower bound is `10^(-14)`, a desired accuracy up to this boundary can be specified with the optional `epsp` argument.
See documentation on Accuracy.

Via `squares_cdf_approx()` this function implements equation (17) from:

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example*

https://arxiv.org/abs/1710.06642
  
"""
function squares_pvalue_approx(T_obs::Real, L::Integer, epsp::Real = 0)

    return 1 - squares_cdf_approx(T_obs, L, epsp)
end

function squares_pvalue_approx(T_obs::Real, Ns::AbstractArray,  epsp::Real = 0)

    return 1 - squares_cdf_approx(T_obs, Ns, epsp)
end