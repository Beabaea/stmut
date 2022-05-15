## code to prepare `DATASET` dataset goes here

# read in external data for package testing
files <- list.files(path = "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotIndex", pattern = ".txt",full.names = TRUE, recursive = FALSE)
data1 <- read.csv("/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/filtered_feature_bc.csv", header = TRUE)
data2 <- read.csv("/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/Graph-Based.csv",header = TRUE)
path <- "/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/tissue_positions_list.csv"

# add test data to the development package
usethis::use_data(files, compress = "xz", overwrite = TRUE)
usethis::use_data(data1, compress = "xz", overwrite = TRUE)
usethis::use_data(data2, compress = "xz", overwrite = TRUE)
usethis::use_data(path, overwrite = TRUE)



# test data for sptMutCt function
write.csv(df[[1]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotBC.csv", row.names = FALSE)
write.csv(df[[2]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotBClster.csv", row.names = FALSE)

index <- read.csv("/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotBC.csv", header = TRUE)
files2 <- list.files(path = "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/mpileup",pattern = "MpileupOutput_RNA.txt", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
usethis::use_data(index, compress = "xz", overwrite = TRUE)
usethis::use_data(files2, compress = "xz", overwrite = TRUE)
write.csv(df1[[1]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotRdCtFinal.csv", row.names = FALSE)
write.csv(df1[[2]], "/Users/limin/limin_practice/Rprojects/stmut/inst/extdata/spotMutReadCount.csv", row.names = FALSE)





