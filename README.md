
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

This is a basic example of how to run sptBClstRds function:

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

An example on how to run sptMutCt()

``` r
library(stmut)
index <- read.csv("./inst/extdata/spotBC.csv", header = TRUE)
files2 <- list.files(path = "./inst/extdata/mpileup",pattern = "MpileupOutput_RNA.txt", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
df1 <- sptMutCt(index = index, files = files2)

head(df1[[1]])
#>        spot RreadC MreadC Treads groups            barcode
#> 267 spot266      0      5      5     G5 CCTCCTGTTGTGTCGT-1
#> 260 spot259      0      4      4     G4 CCGTAGGGTTGTTTAC-1
#> 663 spot662      0      3      3     G3 TGATGTCAATTAAGTG-1
#> 213 spot212      0      2      2     G2 CAGAACTTAGCCCTCT-1
#> 700 spot699      0      2      2     G2 TTACACGATCTGCGAC-1
#> 1   spot000      0      0      0     G0 AAACAAGTATCTCCCA-1
head(df1[[2]])
#>                  index groups
#> 267 CCTCCTGTTGTGTCGT-1     G5
#> 260 CCGTAGGGTTGTTTAC-1     G4
#> 663 TGATGTCAATTAAGTG-1     G3
#> 213 CAGAACTTAGCCCTCT-1     G2
#> 700 TTACACGATCTGCGAC-1     G2
#> 1   AAACAAGTATCTCCCA-1     G0
```

An example to run nonZeRdCts

``` r
library(stmut)
d1 <- read.csv("./inst/extdata/spotBClster.csv", header = TRUE)
d2 <- read.csv("./inst/extdata/spotRdCtFinal.csv", header = TRUE)

dim(d1)
#> [1] 10  6
dim(d2)
#> [1] 700   6

df3 <- nonZeRdCts(df1 = d1, df2 = d2)
head(df3[[1]])
#>              barcode    spot   cluster TotalRDs array_row array_col RreadC
#> 1 CAGAACTTAGCCCTCT-1 spot212 Cluster 4     5019        53        85      0
#> 2 CCGTAGGGTTGTTTAC-1 spot259 Cluster 3    10758        60       104      0
#> 3 CCTCCTGTTGTGTCGT-1 spot266 Cluster 3    19768        62        96      0
#> 4 TGATGTCAATTAAGTG-1 spot662 Cluster 3     8780        59       103      0
#> 5 TTACACGATCTGCGAC-1 spot699 Cluster 4     3699        62       110      0
#>   MreadC Treads groups
#> 1      2      2     G2
#> 2      4      4     G4
#> 3      5      5     G5
#> 4      3      3     G3
#> 5      2      2     G2
head(df3[[2]])
#>              barcode    spot   cluster TotalRDs array_row array_col RreadC
#> 1 CAGAACTTAGCCCTCT-1 spot212 Cluster 4     5019        53        85      0
#> 2 CCGTAGGGTTGTTTAC-1 spot259 Cluster 3    10758        60       104      0
#> 3 CCTCCTGTTGTGTCGT-1 spot266 Cluster 3    19768        62        96      0
#> 4 TGATGTCAATTAAGTG-1 spot662 Cluster 3     8780        59       103      0
#> 5 TTACACGATCTGCGAC-1 spot699 Cluster 4     3699        62       110      0
#>   MreadC Treads groups
#> 1      2      2     G2
#> 2      4      4     G4
#> 3      5      5     G5
#> 4      3      3     G3
#> 5      2      2     G2
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
