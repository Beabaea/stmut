## code to prepare `DATASET` dataset goes here

# read in external data for package testing
files <- list.files(path = "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotIndex", pattern = ".txt", full.names = TRUE, recursive = FALSE)
data1 <- read.csv("/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/filtered_feature_bc.csv", header = TRUE)
data2 <- read.csv("/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/Graph-Based.csv", header = TRUE)
path <- "/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/tissue_positions_list.csv"

# add test data to the development package
usethis::use_data(files, compress = "xz", overwrite = TRUE)
usethis::use_data(data1, compress = "xz", overwrite = TRUE)
usethis::use_data(data2, compress = "xz", overwrite = TRUE)
usethis::use_data(path, overwrite = TRUE)



# # test data for sptMutCt function
# write.csv(df[[1]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotBC.csv", row.names = FALSE)
# write.csv(df[[2]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotBClster.csv", row.names = FALSE)

index <- read.csv("/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotBC.csv", header = TRUE)
files2 <- list.files(path = "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/mpileup", pattern = "MpileupOutput_RNA.txt", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
usethis::use_data(index, compress = "xz", overwrite = TRUE)
usethis::use_data(files2, compress = "xz", overwrite = TRUE)
# write.csv(df1[[1]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotRdCtFinal.csv", row.names = FALSE)
# write.csv(df1[[2]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotMutReadCount.csv", row.names = FALSE)


# test data for nonZeRdCts function
d1 <- read.csv("/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotBClster.csv", header = TRUE)
d2 <- read.csv("/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotRdCtFinal.csv", header = TRUE)
usethis::use_data(d1, compress = "xz", overwrite = TRUE)
usethis::use_data(d2, compress = "xz", overwrite = TRUE)
# write.csv(df3[[1]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/NonZeroRdSpotIndex.csv", row.names = FALSE)
# write.csv(df3[[1]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/MutantSpotIndex.csv", row.names = FALSE)


# test data for spotSummary function
df <- read.csv("/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/NonZeroRdSpotIndex.csv", header = TRUE)
path1 <- "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/mpileup/"
usethis::use_data(df, compress = "xz", overwrite = TRUE)
usethis::use_data(path1, compress = "xz", overwrite = TRUE)


# test data for wtArmMedianOne function
cnr <- read.table("/Volumes/Bastian/Limin/Ji_data/ST_paper/mergedP6/cnr/spot1_rep1.cnr", header = TRUE)
centm <- read.table("/Volumes/Bastian/Limin/reference/hg38_centromereSimple.bed")
usethis::use_data(cnr, overwrite = TRUE)
usethis::use_data(centm, overwrite = TRUE)

# test data for CtArmGenes function
cdt <- read.table("/Volumes/Bastian/Limin/Ji_data/ST_paper/mergedP6/cdt/p6GrpWtClstRdsMissRmdCtdNonTumMedian.cdt", sep = "\t", header = TRUE)
cdt <- cdt[-c(1,2),-c(1,4)] #
cdt <- cdt[,1:5]
data3 <- read.table("/Volumes/Bastian/Limin/Ji_data/ST_paper/mergedP6/armWtedOne/summary.txt", sep = "\t", header = TRUE)
usethis::use_data(cdt, compress = "xz", overwrite = TRUE)
usethis::use_data(data3, compress = "xz", overwrite = TRUE)


# test data for cdt_filt_sort function
data4 <- read.csv("/Volumes/Bastian/Limin/Ji_data/ST_paper/mergedP6/cdt/CtArmGenSummary.csv", header = TRUE)
usethis::use_data(data4, compress = "xz", overwrite = TRUE)

