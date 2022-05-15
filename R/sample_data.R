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




#' Sample Data for sptMutCt Function
#'
#' @description First input for sptMutCt
#'
#' @format A data frame with 10 rows and 3 variables
#' \describe{
#'   \item{barcode}{the spatial barcode of spots}
#'   \item{spot}{the spot name}
#'   \item{TotalRDs}{the total number reads of each spot}
#' }
"index"


#' Sample Data for sptMutCt Function
#'
#' @description Second input for sptBClstRds
#'
#' @format A vector containing all spots samtools mpileup output path
#' \describe{
#'   \item{files2}{a character}
#' }
"files2"



#' Sample Data for nonZeRdCts Function
#'
#' @description First input for nonZeRdCts
#'
#' @format A data frame with 10 rows and 6 variables
#' \describe{
#'   \item{barcode}{the spatial barcode of spots}
#'   \item{cluster}{the cluster the spot belongs to}
#'   \item{spot}{the spot name}
#'   \item{TotalRDs}{the total number of reads of each spot}
#'   \item{array_row}{spot row position}
#'   \item{array_col}{spot col position}
#' }
"d1"

#' Sample Data for nonZeRdCts Function
#'
#' @description Second input for nonZeRdCts
#'
#' @format A data frame with 5 rows and 6 variables
#' \describe{
#'   \item{spot}{the spot name}
#'   \item{RreadC}{reference read counts}
#'   \item{MreadC}{Mutant reads counts}
#'   \item{Treads}{the sum of RreadC and MreadC}
#'   \item{groups}{if MreadC == 0, the spot will given "G0", if MreadC == 1, the spot will given "G1",
#'   .... if MreadC > 9, the spot will given "G10"}
#'   \item{barcode}{the spatial barcode of spots}
#' }
"d2"





#' Sample Data for spotSummary Function
#'
#' @description First input for spotSummary
#'
#' @format A data frame with 5 rows and 10 variables
#' \describe{
#'   \item{barcode}{the spatial barcode of spots}
#'   \item{spot}{the spot name}
#'   \item{cluster}{the cluster the spot belongs to}
#'   \item{TotalRDs}{the total number of reads of each spot}
#'   \item{array_row}{spot row position}
#'   \item{array_col}{spot col position}
#'   \item{RreadC}{reference read counts}
#'   \item{MreadC}{Mutant reads counts}
#'   \item{Treads}{the sum of RreadC and MreadC}
#'   \item{groups}{if MreadC == 0, the spot will given "G0", if MreadC == 1, the spot will given "G1",
#'   .... if MreadC > 9, the spot will given "G10"}
#' }
"df"

#' Sample Data for spotSummary Function
#'
#' @description Second input for spotSummary
#'
#' @format A path to spot Mpileup file
#' \describe{
#'   \item{path}{a character}
#' }
"path1"
