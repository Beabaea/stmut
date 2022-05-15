
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

## Notes before running the package.

This package has a series of 4 functions: sptBClstRds, sptMutCt,
nonZeRdCts, spotSummary. <br /> You have to run each function
sequentially to obtain the ideal output, or You have your input of each
function strictly match the input format of our examples.

Preparations – 5 files from spaceranger output: <br />
filtered_feature_bc.csv <br /> Graph-Based.csv <br />
possorted_genome_bam.bam <br /> spatial/tissue_positions_list.csv <br />
raw_feature_bc_matrix/barcodes.tsv.gz <br />

spotIndex: one spatial barcode one .txt file, thie file name is the spot
name you choose. <br />

spotBam: <https://github.com/10XGenomics/subset-bam> <br />
subset-bam_linux –bam possorted_genome_bam.bam –cell-barcodes
spot000.txt –out-bam spot000.bam samtools index spot000.bam

Count each spot ref and mut reads using Mpileup_RNA.pl (samtools
mpileup): <br /> perl Mpileup_RNA.pl somaticSNPsMML.txt
spot000/MpileupOutput_RNA.txt

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

An example to run spotSummary

``` r
library(stmut)
df <- read.csv("./inst/extdata/NonZeroRdSpotIndex.csv", header = TRUE)
path1 <-"./inst/extdata/mpileup/"

df4 <- spotSummary(df = df, path1 = path1)

head(df4)
#>              barcode    spot   cluster TotalRDs array_row array_col RreadC
#> 1 TTACACGATCTGCGAC-1 spot699 Cluster 4     3699        62       110      0
#> 2 CAGAACTTAGCCCTCT-1 spot212 Cluster 4     5019        53        85      0
#> 3 CCGTAGGGTTGTTTAC-1 spot259 Cluster 3    10758        60       104      0
#> 4 TGATGTCAATTAAGTG-1 spot662 Cluster 3     8780        59       103      0
#> 5 CCTCCTGTTGTGTCGT-1 spot266 Cluster 3    19768        62        96      0
#>   MreadC Treads groups propM2Total propM2SumRM MoreThan1Mrd    Score
#> 1      2      2     G2    5.406867           1            1 6.406867
#> 2      2      2     G2    3.984858           1            1 4.984858
#> 3      4      4     G4    3.718163           1            1 4.718163
#> 4      3      3     G3    3.416856           1            1 4.416856
#> 5      5      5     G5    2.529340           1            1 3.529340
#>                                                  GenesWithMutRead
#> 1                        WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;
#> 2                        WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;
#> 3                        WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;
#> 4                                           RFNG-17:82049969 G>A;
#> 5 WNT4-1:22120068 C>T;WNT4-1:22120070 C>T;CTNNBIP1-1:9850735 C>T;
```

Column meaning of df4: <br /> column 1: Spatial barcode <br /> column 2:
spot name <br /> column 3: cluster in Loupe.Loupe file <br /> column 4:
Total reads each spot has <br /> column 5: row position in
position_lists_bc.csv <br /> column 6: col position in
position_lists_bc.csv <br /> column 7: Reference read counts <br />
column 8: Mutant read counts <br /> column 9: Sum of Ref read count and
Mut read count (sum of col 7 and col 8) <br /> column 10: The number (
\< 10) in this column indicates how many mutant reads in a spot. For
example, G9 means there are 9 mutant reads in that spotl G10 means there
are more than 10 mutant reads in the spot. <br /> column 11: proportion
of mutant to Total reads (col8/col4) <br /> column 12: proportion of
mutant to sum of ref and mut reads (col8/col9) <br /> column 13: A
binary column, when the spot has \> 1 mutant read, it is 1, otherwise,
it is 0. <br /> column 14: A score calculated based on the ref and
mutant read counts. The higher the score is, the more likely it is a
tumor spot. The table is ordered by this column. <br /> column 15:
Details on which genes containing what mutations in which spots. <br />
