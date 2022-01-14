# RunStatistics.jl

A package implementing the exact evaluation of the cumulative distribution function of the `Squares test statistic` ``T`` as defined in 

Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. [doi:10.1016/j.jspi.2011.04.022](https://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](https://arxiv.org/abs/1005.3233)

This package also includes an implementation of an approximation of this cumulative for the more general case of large numbers of observations, as derived in 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](https://arxiv.org/abs/1710.06642)

This code is based on the [original implementation](https://github.com/fredRos/runs) by Frederik Beaujean in c++ and mathematica.

Much of the following explanations is taken from the above two papers.
## Introduction 
---

One of the most common tasks in scientific inference is comparing observations and model predictions. Based on this comparison, the hypothesized model may be either accepted or rejected. In the latter case usually an improved model is sought. The comparison between observations and the new model is then repeated until a satisfactory model has been constructed.

In model validation the goal is to provide quantitative test procedures.
The standard approach consists of defining a scalar function of the data ``D``, called [*test statistic*](https://en.wikipedia.org/wiki/Test_statistic) ``T (D)``, such that a large value of ``T`` indicates a large deviation of the data from the expectations under the hypothesized model ``\mathcal{H}``. Correspondingly, small ``T`` is seen as good agreement.
## The Squares test statistic
---

The `Squares test statistic` or `Squares statistic` for short, in the following denoted with ``T``, is a test statistic sensitive to local deviations of the data from expectations within an ordered data set[^1].

It supplements the classic [``\chi^2`` test](https://en.wikipedia.org/wiki/Chi-squared_test) which ignores the ordering of observations and provides additional sensitivity to local deviations from expectations. 

The `Squares statistic` ``T`` can be defined for data that follows any symmetric distribution[^1], but in this package only data with a [gaussian probability distribution](https://en.wikipedia.org/wiki/Normal_distribution) is considered:

```math
\begin{align}
X_i \sim \mathcal{N}(\mu_i, \sigma_i^2)
\end{align}
```
The hypothesis ``\mathcal{H}`` for the data is:

- All observations ``\{X_i\}`` are independent. 
- Each observation is normally distributed, ``X_i \sim \mathcal{N}(\mu_i, \sigma^2_i)``
- Mean ``\mu_i`` and variance ``\sigma^2_i`` are known.

``T`` is based on `runs` of weighted deviations from a mean value, observed in samples ``X_i`` from independent normal distributions. 

A `run` in this context refers to a sequence of observations that share a common attribute commonly called a `success`. Here an observation is called a *success*, ``S``, if the observed value exceeds the expected value. Similarly, an expected value exceeding the observation is considered a *failure*, ``F``.

``T`` is formally defined via:

---
-  Split the data ``{X_i}`` into runs. Keep the success runs and ignore the
    failure runs. Denote by ``A_j = \{X_{j_1} ,X_{j_2}, ...\}`` the set of 
    observations in the ``j``-th success run.

-  Associate a weight ``\omega(A_j)`` with each success run:

```math
\begin{align}
\omega(A_j) \equiv \chi_{run,j}^2 = \displaystyle\sum_i\frac{(X_i-\mu_i)^2}{\sigma_i^2}
\end{align}
```
- Choose ``T`` as the largest weight of any run in the entire sequence of observed data:

```math
\begin{align}
T \equiv \max_j \omega(A_j)
\end{align}
```
---

Note that the choice for the weight ``\omega(A_j)`` implemented in this package is only one of the most significant ones, more general options are available (see section 1. of [^1]).

Consider for example an observed data sequence of: 

```math
SSSFFSFFFSSF 
```
where ``S`` denotes a *success*, a value above the expected value, and ``F`` a *failure*. In accordance with the above steps, only the success runs are considered:

```math
\underbrace{\mathbf{SSS}}_{A_1}~~FF~\underbrace{\mathbf{S}}_{A_2}~FFF~\underbrace{\mathbf{SS}}_{A_3}~F 
```
For each of the three success runs ``A_j`` observed in this example, the weight ``\omega(A_j)=\chi_{run,j}^2`` is calculated. The value ``T_{obs}``, denoting the observed value of ``T`` in this sequence of data is then the maximum of the three weights ``T_{obs} = \max_j \omega(A_j)``.


## Interpreting the Squares statistic
---

To facilitate the interpretation of ``T`` (how large is too large?), it is useful to introduce the [*p-value*](https://en.wikipedia.org/wiki/P-value) ``p``. Assuming ``\mathcal{H}``, the p-value is defined as the tail area probability to randomly sample a value of ``T`` larger than or equal to ``T_{obs}``, the value of ``T`` observed in the data:
```math
\begin{align}
p \equiv P (T ≥ T_{obs} ~|~ \mathcal{H}) 
\end{align}
```
If ``\mathcal{H}`` is correct and all parameters are fixed, then ``p`` is a random variable with uniform distribution on ``[0, 1]``. An incorrect model will typically yield smaller values of ``p``.

### Approximation for large numbers of observations
---

The cost for calculating the exact `p-value` for the Squares statistic as described in the initial paper[^1], scales with the number ``N`` of observations in a sequence of data like ``\exp[N^{\frac{1}{2}}]/N`` and quickly grows to large for ``N \gtrsim 80``. 

The authors derived an approximation for large numbers of data in the follow-up paper[^2].

The underlying principle is to split the (long) total sequence of observed data into shorter sequences, for which the p-value can be computed exactly. An approximate p-value ``p`` for the entire observed data sequence can then be extrapolated. The approximation is constructed to have high accuracy in the region of interest, for small values for ``p``. 

For a comprehensive explanation see section II. of [^2].
## Using RunStatistics.jl
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

#TODO: Explain usage of package




## Details of computation
---

In the following some more in-depth information on the calculations performed with `RunStatistics` is given.

As during the derivation of the p-value for ``T`` in [^1], the quantity that is actually being computed here is ``P(T < T_{obs} | N)``. This is the value of the cumulative distribution function of ``T`` at the value ``T_{obs}`` observed in a sequence of ``N`` datapoints. 

The p-value then is obtained as ``p = 1- P(T < T_{obs} | N)``.

The central calculation in this package implements Equation (16) from [^1]:

```math
P(T < T_{obs} | N) = \sum_{r = 1}^{N}\sum_{M = 1}^{M_{max}}\sum_{\pi} X(T_{obs}, N)
```

The full derivation of ``P(T < T_{obs} | N)`` and an explanation for the parameters in the above equation can be found in section 2. of [^1]. 

For this manual, suffice it to say that the only input parameters that need to be known are ``T_{obs}`` and ``N``, the total number of observed data points. The other parameters are then calculated from them.

The computationally expensive operation is the sum over ``\pi``. Here ``\pi`` denotes the set of inequivalent sequences of success and failure runs. Given a total number ``r`` of successes in a sequence of ``N`` observations with ``M`` success runs, ``\pi`` is in one-to-one correspondence with the set of [*integer partitions*](https://en.wikipedia.org/wiki/Partition_(number_theory)) of ``r`` into ``M`` summands.

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

*A. Zoghbi: Algorithms for generating integer partitions, Ottawa (1993)* find it . 
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
