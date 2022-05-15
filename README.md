
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
# ## basic example code
# files <- list.files(path = "./inst/extdata/spotIndex", pattern = ".txt", full.names = TRUE, recursive = FALSE)
# data1 <- read.csv("./inst/extdata/filtered_feature_bc.csv", header = TRUE)
# data2 <- read.csv("./inst/extdata/Graph-Based.csv", header = TRUE)
# path <- "./inst/extdata/tissue_positions_list.csv"
# df <- sptBClstRds(files=files,data1=data1, data2=data2, path=path)
# 
# head(df[[1]])
# head(df[[2]])
?sptBClstRds
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

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
