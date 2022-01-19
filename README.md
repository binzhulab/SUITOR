# SUITOR: selecting the number of mutational signatures through cross-validation
<br/>

### Introduction
For the  _de novo_ mutational signature analysis, estimating the correct number of signatures is the crucial starting point, since it influences all the downstream steps, including extraction of signature profiles, estimation of signature activities and classification of tumors based on the andestimated activities. Here we present an **R** package `SUITOR`, an unsupervised cross-validation tool to select the optimal number of signatures. This tutorial introduces the usage of two main functions `suitor()` and `suitor_Extract_WH()` as follows.

<br/>

### Installation
To install from Github, use the devtools R package:
```r
if (!requireNamespace("devtools", quietly = TRUE))  
	install.packages("devtools")
devtools::install_github("binzhulab/SUITOR/source")
```
Alternatively, download the package and follow the steps below. Download SUITOR_1.0.0.tar.gz (for Unix) or SUITOR_1.0.0.zip (for Windows, R version >= 4.0). To install SUITOR on Unix, enter the command (without the quotes) from a Unix prompt:
```
R CMD INSTALL SUITOR_1.0.0.tar.gz -l path_to_install_package
```
Alternatively, SUITOR_1.0.0.tar.gz (for Unix) or SUITOR_1.0.0.zip (for Windows, R version >= 4.0) from the [Github page](https://github.com/binzhulab/SUITOR) are available and one may use the following commands:
```
install.packages("SUITOR_1.0.0.tar.gz", repose = NULL, type = "source")
install.packages("SUITOR_1.0.0.zip", repose = NULL, type = "win.binary")
```
Once the installation is successful, it can be loaded on **R** by calling 
```
library(SUITOR)
```
<br/>

### Run `suitor()` function
The main function `suitor(data, op)` is to select the number of mutational signatures based on cross-validation. The first argument of the function is `data`. It could be an **R**  `dataframe` or `matrix` containing mutational catalog whose elements are non-negative counts. We set each column of `data` corresponding to a tumor (or sample) while its rows representing a mutation type. For example the dimension of `data` is 96 by _N_ for single base substitution where _N_ is the number of tumors. An example data (`SimData`) with dimension 96 by 300 is available in the paackage for illustrative purpose. 
First, we start with the default option:
```
> data(SimData)
> re <- suitor(SimData)
> re$rank
[1] 8
```
It computes the estimated optimal rank stored at `re$rank` and generates the cross validation error plot. 

#### `suitor()` Options 
Since SUTIOR is based on cross-validation and the Expectation Conditional Maximization (ECM) algorithm, it is necessary to set a list of tuning parameters which control the fitting process.

| Name   |      Description      |  Default Value |
|----------|:-------------|:------|
|`min.value` | Minimum value of matrix before factorizing | 1e-4 |
|`min.rank` | Minimum rank | 1 |
|`max.rank` | Maximum rank | 10 |
|`k.fold` | Number of folds | 10 |
|`em.eps`| EM algorithm stopping tolerance | 1e-5 |
|`max.iter` | Maximum number of iterations in EM algorithm | 2000 |
|`n.seeds` | Number of seeds (starting points) | 30 |
|`n.cores` | Number of cores to use | 1 |
|`get.summary` | 0 or 1 to create summary results | 1 |
|`plot` | 0 or 1 to produce an error plot | 1 |
|`print` | 0 or 1 to print info (0=no printing) | 1 |
|`seeds` | Vector of seeds (takes precedence over n.seeds) | NULL |

`min.value` is a small number added to the `data` matrix for stable computation of non-negative matrix factorization. For a given number of signatures or called rank _r_  (`min.rank` &le; _r_ &le; `max.rank`), the `data` matrix is divided into `k.fold` parts for the cross-validation. The default value of the maximal rank `max.rank` is 10 but it can be changed depending on the cancer type. The default value of the number of fold K (`k.fold`) is 10 and it can be modified depending on the computer resources.

Since the ECM algorithm may converge to a local saddle point, SUITOR tries multiple initial values for _**W**_~0~ and _**H**_~0~. For this purpose, the number of seeds (`n.seeds`) or a vector of seeds (`seeds`) is used. For example, when setting the number of seeds `n.seeds = 30`, 30 seeds are randomly generated while setting a vector of seeds as `seeds = 1:30` defines the seeds as 1, 2, ..., 30. Note that `seeds` takes precedence over `n.seeds`; in fact, if `seeds = 1:30` is used, `n.seeds = 30` will be ignored. Although the default `n.seeds` is set to 30, it can be increased depending on the size of the `data` matrix and/or computational resources.

For the ECM algorithm, the default value of the maximal iteration `max.iter` is set to 2000. It is possible for some cases to reach the maximal iteration, for which the function would produce a warning message. Overall, we recommend a two-stage approach where the user would  run `suitor()` with the default option first and then narrow down the set of plausible ranks (`min.rank` &le; _r_ &le; `max.rank`) with more seeds (`seeds`) and a larger number of maximal iteration (`max.iter`) if necessary.

To ease the computation burden, SUITOR supports the parallel computing for Windows and UNIX machines by specifying `n.cores` greater than 1. To check whether parallel computing is available for the computer,
```
> detectCores()
```
If `detectCores()` returns a value greater than 1, that return value may be used for `n.cores`.


To run the function `suitor()` with different option values, we create a list of options as follows. Note that the name of the option `list` should be matched with elements in the above option table.
```r
> OP <- list(min.rank = 5, max.rank = 13, k.fold = 5, n.seeds = 50, 
	     get.summary = 0)
> re2 <- suitor(data = SimData, op = OP)
> Summary2 <- getSummary(re2$all.results, ncol(SimData))
> Summary2$rank
[1] 8
> plotErrors(Summary2$summary)
```
If `get.summary` is set to 0, `suitor()` only produce a matrix containing all possible results. In that case, as shown below one needs to use the `getSummary(obj, NC, NR)` function to compute the estimated optimal rank (`Summary1$rank`), where `obj` is the matrix containing all results from `suitor()`, `NC` and `NR` are the numbers of columns and rows respectively in the input `data` for `suitor()`. Please note that `NR` has a default value of 96 for single base substitution signature analysis. The `plotErrors()` function could be used to draw a cross validation error plot.
[comment]: # (plotErrors.pdf Here!!)

### Run `suitor_extract_WH()` function
Once the optimal number of signature or called rank is estimated by `suitor()`, we can extract the signature profiles  _**W**_~est~ and activities _**H**_~est~ with the function `suitor_extract_WH(data, rank, op)`. As in the `suitor()` function, the input `data` is a data frame or matrix containing mutational catalog whose elements are non-negative counts. A non-negative integer `rank` is the number of mutational signatures to be extracted. The possible option values are summarized in the following table and they can be used in the same manner as `suitor()`.
| Name   |      Description      |  Default Value |
|----------|:-------------|:------|
|`min.value` | Minimum value of matrix before factorizing | 1e-4 |
|`n.seeds` | Number of seeds (starting points) | 30 |
|`n.cores` | Number of cores to use for parallel computing | 1 |
|`print` | 0 or 1 to print info (0=no printing) | 1 |
|`seeds` | Vector of seeds (takes precedence over n.seeds) | NULL |

To run the `suitor_extract_WH()` function, 
```
> re <- suitor(SimData)
> re$rank
[1] 8
> Extract <- suitor_extract_WH(SimData, re$rank)
seed 30
> head(Extract$W)
	  denovo A     denovo B     denovo C     denovo D    denovo E
A[C>A]A 0.04983920 6.864379e-20 9.210142e-04 5.021053e-03 0.010287996
A[C>A]C 0.03619912 2.909597e-03 6.688993e-04 4.062685e-07 0.007367939
A[C>A]G 0.01907898 1.760262e-03 5.568929e-20 5.822657e-20 0.001869262
A[C>A]T 0.03533304 8.483854e-04 5.568929e-20 2.024296e-03 0.008467782
C[C>A]A 0.08706398 7.221585e-03 5.568929e-20 1.323915e-03 0.010022022
C[C>A]C 0.10698556 6.864379e-20 5.568929e-20 2.885607e-12 0.003802372
	    denovo F     denovo G     denovo H
A[C>A]A 6.236367e-20 1.121825e-03 5.495331e-20
A[C>A]C 6.236367e-20 1.053258e-03 9.980362e-04
A[C>A]G 6.236367e-20 2.401631e-07 1.528276e-10
A[C>A]T 4.031693e-07 5.511649e-20 2.594604e-03
C[C>A]A 1.427441e-04 1.903321e-03 5.495331e-20
C[C>A]C 1.287044e-06 5.511649e-20 1.492892e-03

> Extract$H[,1:3]
	         [,1]         [,2]         [,3]
denovo A 2.803011e+00 2.050670e+00 7.926640e-12
denovo B 1.508570e+01 2.105153e+00 1.815314e+01
denovo C 9.011722e+01 1.227970e+01 2.238425e+00
denovo D 8.467579e-13 6.883381e+01 5.620215e+01
denovo E 2.281342e+01 6.971187e-13 2.699901e+00
denovo F 4.592150e+01 4.630038e+01 1.337674e-04
denovo G 1.431012e+01 2.083340e+01 8.982586e+01
denovo H 8.595290e+01 1.260151e+01 1.788483e+01
```

For more information please refer to the [user guide](link_to_userguide)
