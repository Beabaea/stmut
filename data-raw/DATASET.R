## code to prepare `DATASET` dataset goes here

# read in external data for package testing
files <- list.files(path = "/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/spotIndex", pattern = ".txt",full.names = TRUE, recursive = FALSE)
data1 <- read.csv("/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/filtered_feature_bc.csv", header = TRUE)
data2 <- read.csv("/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/Graph-Based.csv",header = TRUE)
path <- "/Volumes/Bastian/Limin/Ji_data/ST_paper/rFunctions/data/tissue_positions_list.csv"

# add test data to the development package
usethis::use_data(files, compress = "xz", overwrite = TRUE)
usethis::use_data(data1, compress = "xz", overwrite = TRUE)
usethis::use_data(data2, compress = "xz", overwrite = TRUE)
usethis::use_data(path, overwrite = TRUE)
