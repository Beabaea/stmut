
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stmut

<!-- badges: start -->
<!-- badges: end -->

The goal of stmut is providing a series of functions to visualize which
genes have somatic mutations in which spots when working on sptial
transcriptomics data.

## Installation

You can install the development version of stmut from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("limin321/stmut")
```

## Instruction before using this package.

Have your data ready:

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(stmut)
## basic example code
files <- list.files(path = "./inst/extdata/spotIndex", pattern = ".txt", full.names = TRUE, recursive = FALSE)
data1 <- read.csv("./inst/extdata/filtered_feature_bc.csv", header = TRUE)
data2 <- read.csv("./inst/extdata/Graph-Based.csv", header = TRUE)
path <- "./inst/extdata/tissue_positions_list.csv"
df <- sptBClstRds(files=files,data1=data1, data2=data2, path=path)

head(df[[1]])
#>              barcode    spot TotalRDs
#> 1 AAACAAGTATCTCCCA-1 spot000    22655
#> 2 AGATACTCAAGATCGA-1 spot100    10631
#> 3 AGATCTCAGGTGTGAT-1 spot101      797
#> 4 AGATGACTCGCCCACG-1 spot102      260
#> 5 AGCAACATATCTTATT-1 spot103      180
#> 6 CAGAACTTAGCCCTCT-1 spot212     5019
head(df[[2]])
#>              barcode   cluster    spot TotalRDs array_row array_col
#> 1 AAACAAGTATCTCCCA-1 Cluster 3 spot000    22655        50       102
#> 2 AGATACTCAAGATCGA-1 Cluster 4 spot100    10631        46        82
#> 3 AGATCTCAGGTGTGAT-1 Cluster 2 spot101      797        29        75
#> 4 AGATGACTCGCCCACG-1 Cluster 2 spot102      260        21        55
#> 5 AGCAACATATCTTATT-1 Cluster 2 spot103      180        41        93
#> 6 CAGAACTTAGCCCTCT-1 Cluster 4 spot212     5019        53        85
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
