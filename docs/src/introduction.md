# Introduction 
---
This page provides a summarized explanation for the *Squares test statistic* ``T`` and when and in what context it can be used. For a comprehensive derivation and explanation see [^1] and [^2].

```@contents
Pages = ["introduction.md"]
Depth = 3
```

## Motivation
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

The cost for calculating the exact `p-value` for the Squares statistic as described in the initial paper[^1], scales with the number ``N`` of observations in a sequence of data like ``\exp[N^{\frac{1}{2}}]/N`` and quickly grows too large for ``N \gtrsim 80``. 

The authors derived an approximation for large numbers of data in the follow-up paper[^2].

The underlying principle is to split the (long) total sequence of observed data into shorter sequences, for which the p-value can be computed exactly. An approximate p-value ``p`` for the entire observed data sequence can then be extrapolated. The approximation is constructed to have high accuracy in the region of interest, for small values for ``p``. 

For a comprehensive explanation see section II. of [^2].


[^1]: Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. [doi:10.1016/j.jspi.2011.04.022](https://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](https://arxiv.org/abs/1005.3233)

[^2]: Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](https://arxiv.org/abs/1710.06642)
