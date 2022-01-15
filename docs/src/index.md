# RunStatistics.jl

A package implementing the exact evaluation of the cumulative distribution function of the *Squares test statistic* ``T`` as defined in 

Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. [doi:10.1016/j.jspi.2011.04.022](https://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](https://arxiv.org/abs/1005.3233)

This package also includes an implementation of an approximation of this cumulative for the more general case of large numbers of observations, as derived in 

Frederik Beaujean and Allen Caldwell. *Is the bump significant? An axion-search example* [arXiv:1710.06642](https://arxiv.org/abs/1710.06642)

This code is based on the [original implementation](https://github.com/fredRos/runs) by Frederik Beaujean in c++ and mathematica.

## Package Features
---

- Calculate the [*p-value*](https://en.wikipedia.org/wiki/P-value) for the *Squares test statistic*[^1] ``T`` observed in ``N`` independent trials of gaussian uncertainty
- Improved performance and readability over the [original implementation](https://github.com/fredRos/runs) thanks to Julia 
- Inbuilt implementation for generating [*integer partitions*](https://en.wikipedia.org/wiki/Partition_(number_theory)) of an integer ``n`` into ``k`` parts

## Manual Outline
---

```@contents
Pages = ["introduction.md","guide.md"]
Depth = 3
```

[^1]: Frederik Beaujean and Allen Caldwell. *A Test Statistic for Weighted Runs*. Journal of Statistical Planning and Inference 141, no. 11 (November 2011): 3437–46. [doi:10.1016/j.jspi.2011.04.022](https://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](https://arxiv.org/abs/1005.3233)

## Citing RunStatistic.jl
---

When using `RunStatistics.jl` for research, teaching or similar, please cite the original authors' work:

```
@article{beaujean2011test,
title={A test statistic for weighted runs},
author={Beaujean, Frederik and Caldwell, Allen},
journal={Journal of Statistical Planning and Inference},
volume={141},
number={11},
pages={3437--3446},
year={2011},
publisher={Elsevier}
}

@article{Beaujean:2017eyq,
  author         = "Beaujean, Frederik and Caldwell, Allen and Reimann, Olaf",
  title          = "{Is the bump significant? An axion-search example}",
  year           = "2017",
  eprint         = "1710.06642",
  archivePrefix  = "arXiv",
  primaryClass   = "hep-ex",
  SLACcitation   = "%%CITATION = ARXIV:1710.06642;%%"
}
```