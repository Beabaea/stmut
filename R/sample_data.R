#' Sample Data for sptBClstRds Function
#'
#' @description First input for sptBClstRds
#'
#' @format A vector containing all spotIndex path
#' \describe{
#'   \item{files}{a character}
#' }
"files"


#' Sample Data for sptBClstRds Function
#'
#' @description Second input for sptBClstRds
#'
#' @format A data frame with 36601 rows and 11 variables
#' \describe{ ## document each column/variable; no "" is needed for the variable name, otherwise it will cause warning.
#'   \item{X}{the first colum is the gene name}
#'   \item{AAACAAGTATCTCCCA.1}{spatial barcode of spot1}
#'   \item{AGATACTCAAGATCGA.1}{spatial barcode of spot2}
#'   \item{AGATCTCAGGTGTGAT.1}{spatial barcode of spot3}
#'   \item{AGATGACTCGCCCACG.1}{spatial barcode of spot4}
#'   \item{AGCAACATATCTTATT.1}{spatial barcode of spot5}
#'   \item{CAGAACTTAGCCCTCT.1}{spatial barcode of spot6}
#'   \item{CCGTAGGGTTGTTTAC.1}{spatial barcode of spot7}
#'   \item{CCTCCTGTTGTGTCGT.1}{spatial barcode of spot8}
#'   \item{TGATGTCAATTAAGTG.1}{spatial barcode of spot9}
#'   \item{TTACACGATCTGCGAC.1}{spatial barcode of spot10}
#' }
"data1"

#' Sample Data for sptBClstRds Function
#'
#' @description Third input for sptBClstRds
#'
#' @format A data frame with 10 rows and 2 variables
#' \describe{
#'   \item{Barcode}{the spatial barcode of spots}
#'   \item{Graph.based}{the cluster of 10X Loupe}
#' }
"data2"

#' Sample Data for sptBClstRds Function
#'
#' @description Fourth input for sptBClstRds
#'
#' @format A path to tissue_positions_list.csv
#' \describe{
#'   \item{path}{a character}
#' }
"path"
