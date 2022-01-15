# Guide
---
The following guide provides information on how to install and use the `RunStatistics.jl` package and the functions it provides and on some of their the implementation. 

```@contents
Pages = ["guide.md"]
Depth = 3
```

## Installation
---

To install `RunStatistics.jl`, start Julia and run 

```Julia
julia> using Pkg
julia> pkg"add RunStatistics"
```

To use `RunStatistics.jl` after installation, run 

```Julia
julia> using RunStatistics
```
to gain access to the functions provided in the package.


## #TODO: Explain usage of package
---


## Details of computation
---

In the following some more in-depth information on the calculations performed with `RunStatistics` is given.

As during the derivation of the p-value for ``T`` in [^1], the quantity that is actually being computed here is 

```math
P(T < T_{obs} | N)
```

This is the value of the [cumulative distribution function](https://en.wikipedia.org/wiki/Cumulative_distribution_function) of ``T`` at the value ``T_{obs}`` observed in a sequence of ``N`` datapoints. 

The p-value then is obtained as ``p = 1- P(T < T_{obs} | N)``.

The central calculation in this package implements Equation (16) from [^1]:

```math
P(T < T_{obs} | N) = \sum_{r = 1}^{N}\sum_{M = 1}^{M_{max}}\sum_{\pi} X(T_{obs}, N)
```

The full derivation of ``P(T < T_{obs} | N)`` and an explanation for the parameters in the above equation can be found in section 2. of [^1]. 

For this manual, suffice it to say that the only input parameters that need to be known are ``T_{obs}`` and ``N``, the total number of observed data points. The other parameters are then calculated from them.

The computationally expensive operation is the sum over ``\pi``. Here ``\pi`` denotes the set of inequivalent sequences of success and failure runs. Given a total number ``r`` of successes in a sequence of ``N`` observations with ``M`` success runs, ``\pi`` is in one-to-one correspondence with the set of [*integer partitions*](https://en.wikipedia.org/wiki/Partition_(number_theory)) of ``r`` into ``M`` summands [^2].

When the cumulative of ``T`` is evaluated at ``T_{obs}`` for ``N`` data points, the relevant partitions need to be calculated.

### Partitions
---

In this package, a partition of an integer ``n`` into ``k`` parts is represented with a `Partition()` object. It holds the fields:

    n::Int   
    k::Int
    h::Int          Number of *distinct* parts
    c::Vector{Int}  Multiplicities of parts
    y::Vector{Int}  parts

With these parameters, a partition can be represented in the *multiplicity representation*:
```math
\begin{align}
n = \sum_{i = 2}^{h + 1} c_i \cdot y_i
\end{align}
```
It is important to note the indexing of the arrays holding the parts and their multiplicities. Due to the implementation of the algorithm used to generate partitions, the first elements of `c` and `y` hold a buffer value, and the arrays always have a length of `1 + the maximum possible number of parts`.

So when reading a partition, always use the above equation: ignore the first element of `c` and `y` and do not read beyond `c[h + 1]` and `y[h + 1]`.

To save memory during the evaluation of the cumulative, a partition object is initated and updated in place during the summation over the set of possible partitions.

The function that updates a partition is `next_partition!()` it implements a modified version of *Algorithm Z* from 

*A. Zoghbi: Algorithms for generating integer partitions, Ottawa (1993)*. 
### Approximation for large numbers of data
---
#TODO: discuss choice of n and N

For the approximation of the p-value of the Squares statistic for large numbers of data, this package implements equation (17) from [^2]:

```math
\begin{align}
F(T_{obs} | nN) = \frac{F(T_{obs} | N)^n}{(1 + \Delta(T_{obs}))^{n-1}} \quad \text{for} \quad n\ge 2
\end{align}
```

``F(T_{obs} | L) \equiv P(T < T_{obs} | N)`` denotes the value of the cumulative of ``T`` for a sequence of ``L`` observations. 

So if a total number of ``L`` data points have been observed, choose ``n`` and ``N`` so that ``n \cdot N = L``. The exact value of P(T < T_obs | N) is then calculated and further processed in accordance with the above equation.

The approximate p-value for the data set then is:

```math
\begin{align}
p = 1 - F(T_{obs} | nN)
\end{align}
```

``\Delta(T_{obs})`` is a correction term (see equation (13) in [^2]) whose computation involves a 1D numerical integration. This is performed with the `quadgk()` function from the [`QuadGK.jl`](https://juliapackages.com/p/quadgk) package. When calling the `approx_cumulative()` or `approx_pvalue()` functions, the optional arguments`epsrel` and `epsabs` are the relative and absolute target precision of the integration performed in `quadgk()`. If not specified, the default values of `quadgk()` are used. See [documentation](https://juliamath.github.io/QuadGK.jl/stable/).



[^1]: Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. [doi:10.1016/j.jspi.2011.04.022](https://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](https://arxiv.org/abs/1005.3233)

[^2]: Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](https://arxiv.org/abs/1710.06642)
