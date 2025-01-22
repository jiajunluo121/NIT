What is NIT?
======

NIT is a new emperical BAYES estimation procedure that directly incorporates the side information into the estimation procedure. Our approach relies on an improved nonparametric empirical Bayes deconvolution estimator which utilize a nonparametric integrative Tweedie approach.


How to use this repository?
----------

This repository holds the scripts that reproduce the analysis in the paper[1]. In `example.R`, we provide an example about how to use the NIT estimator. `NIT.R` must be available in the current R working directory when using NIT. Meanwhile, the subfolder `simulations` has the R code to reproduce settings 1-4 in  simulation 1 and settings 1-4 in simulation 2 in the paper. In each case, the code `funcs.R` and `NIT.R` must be available in the current R working directory. You will also need `MOSEK 9.3 or higher` and the associated R interface available in the [`Rmosek`](https://docs.mosek.com/latest/rmosek/index.html) package.

[1.] EMPIRICAL BAYES ESTIMATION WITH SIDE INFORMATION: A NONPARAMETRIC INTEGRATIVE TWEEDIE APPROACH. Jiajun Luo, Trambak Banerjee, Gourab Mukherjee and Wenguang Sun

