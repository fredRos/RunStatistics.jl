# RunStatistics.jl

A package implementing the exact evaluation of the cumulative distribution function of the `Squares test statistic` `T` as defined in 

Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. [doi:10.1016/j.jspi.2011.04.022](http://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](http://arxiv.org/abs/1005.3233)

This package also includes an implementation of an approximation of this cumulative for the more general case of large numbers of observations, as derived in 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](http://arxiv.org/abs/1710.06642)

Much of the following explanations is copied verbatim from the above two papers.
## Introduction 

One of the mose common taks in scientific inference is comparing observations and model predictions. Based on this comparison, the hypothesized model may be either accepted or rejected. In the latter case usually an improved model is sought. The comparison between observations and the new model is then repeated until a satisfactory model has been constructed.

In model validation the goal is to provide quantitative test procedures.
The standard approach consists of defining a scalar function of the data $D$, called [*test statistic* $T (D)$](https://en.wikipedia.org/wiki/Test_statistic), such that a large value of $T$ indicates a large deviation of the data from the expectations under the hypothesized model $\mathcal{H}$. Correspondingly, small $T$ is seen as good agreement.

Let $T_{obs}$ denote the value of $T$ observed in the actual data set. In order to facilitate the interpretation of $T$ (how large is too large?), it is useful to introduce the [*p-value*](https://en.wikipedia.org/wiki/P-value). Assuming $\mathcal{H}$, the p-value is defined as the tail area probability to randomly sample a value of $T$ larger than or equal to $T_{obs}$:

$$
\begin{align}
p \equiv P (T ≥ T_{obs} ~|~ \mathcal{H}) 
\end{align}
$$


If $\mathcal{H}$ is correct and all parameters are fixed, then $p$ is a random variable with uniform distribution on $[0, 1]$. An incorrect model will typically yield smaller values of $p$. This is used to guide model selection. For the same data, different models will give different $p$. Similarly, a different choice of the test statistic produces a different $p$ for the same model and data. Why use different statistics? Because one statistic is sensitive to certain, but not to all properties of the model.

The `Squares test statistic`, in the following denoted with `T`, is a test statistic sensitive to local deviations of the data from expectations within an ordered data set.

It supplements the classic [$\chi^2$ test](https://en.wikipedia.org/wiki/Chi-squared_test) which ignores the ordering of observations and provides additional sensitivity to local deviations from expectations. 

The Squares test statistic, or `Squares statistic` for short, is based on `runs` of weighted deviations from a mean value, observed in samples from independent normal distributions. 

A `run` in this context refers to a sequence of observations that share a common attribute commonly called a `success`. Here an observation is called a *success*, $S$, if the observed value exceeds the expected value. Similarly an expected value exceeding the observation is considered a *failure*, $F$.

The `Squares statistic` $T$ is formally defined in three steps:


1.  Split the data ${X_i}$ into runs. Keep the success runs and ignore the
    failure runs. Denote by $A_j = \{X_{j_1} ,X_{j_2}, ...\}$ the set of 
    observations in the $j$-th success run.

2.  Associate a weight with each success run. The weight $\omega(A_j)$ ought 
    to be chosen such that a large weight indicates large discrepancy between 
    model and observations. A natural choice of the weight function is a 
    convenient one-to-one function of the probability (density) of $A_j$ such as 

$$
\begin{align}
\omega(A_j) = [P (A_j | \mathcal{H})]^{-1} \quad \text{or} \quad \omega(A_j) = −~2\log(P (A_j | \mathcal{H}))
\end{align}
$$


3. Choose $T$ as the largest weight:

$$
\begin{align}
T \equiv \max_j \omega(A_j)
\end{align}
$$


