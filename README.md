# gapit2mashr ![alt text][travisbuild]



[travisbuild]: https://travis-ci.org/Alice-MacQueen/gapit2mashr.svg?branch=master

Convert GAPIT output files to input data frames used by mashr.

The most time-consuming part of data analysis is often the first step: getting the data into the right format to run an analysis. To that end, `gapit2mashr` takes output files from [GAPIT](https://www.maizegenetics.net/gapit "Genome Associated Prediction Integrated Tool"), a commonly used software for genome wide association studies (GWAS), and creates a series of dataframes from them suitable for input into [mashr](https://github.com/stephenslab/mashr "multivariate adaptive shrinkage"), a new and exciting software for testing and estimating many effects in many conditions. Once you have run multiple GWAS in R using GAPIT, `gapit2mashr`  allows you to quickly prepare your results for use in mashr.

# Setup

1. Install the [devtools](https://github.com/r-lib/devtools) package, if you haven't already done so:
```R
install.packages("devtools")
```

2. Install `gapit2mashr` with the following commands:
```R
library(devtools)
devtools::install_github("Alice-MacQueen/gapit2mashr")
```
This command should automatically install any missing dependencies that are available on CRAN.

# Usage

1. Determine the number of SNPs you want to select from each file. A quick and dirty way to estimate the number of SNPs you should select is to divide the number of cells you want in your output data frames by the number of conditions you plan to include (the number of columns in your dataframe). Currently [mashr](https://github.com/stephenslab/mashr "multivariate adaptive shrinkage") recommends dataframes of around one million cells as a maximum. So if you have 40 conditions, for example, you likely want to set `numSNPs` to 25000. You can use `gapit_phenotypes_in_folder` to obtain this quantity and to find the names of the phenotypes in your GAPIT results directory:
```R
phenotypes_vector <- gapit_phenotypes_in_folder(path = "path/to/your/GAPIT/Results", model = "CMLM") 
# currently only CMLM is supported
numSNPs <- 1000000 / length(phenotypes_vector)
```

2. To run mashr on every file in a specific directory, run:
```R
gapit2mashr(path = "path/to/your/GAPIT/Results", numSNPs = numSNPs, S_hat = "Hedges' G", saveoutput = TRUE)
```
By default, this function will not save output files. If you specify `saveoutput = TRUE`, then it will save output files as `.rds` files to the path you specified in the function.

Sometimes, GAPIT can not determine standard errors about the effect sizes of SNPs. When GAPIT does not provide standard errors to use for the S_hat matrix, which contains the standard error in the effect size difference between the reference and alternate allele, `gapit2mashr` currently uses Hedges' g (Hedges & Olkin 1985 p. 86) to calculate S_hat. Hedges's G allows the calculation of both the effect size of the alternate allele, and the confidence interval around the effect size. `gapit2mashr` uses the effect sizes provided by GAPIT to compute the confidence interval calculation. The calculations are:
```R
d = 2*r/sqrt(1-r^2)
d_unbiased = (1-(3/(4*(N-2)-1)))*d
sigma_squared_d_i = (n_i^e + n_i^c)/n_i^e*n_i^c + d_i^2 / 2*(n_i^e + n_i^c)
```
 where r is the effect size, scaled between -1 and 1; n's are the sample sizes of the two experimental groups; N is the total sample size.
