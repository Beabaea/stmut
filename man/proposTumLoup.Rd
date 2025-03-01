% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptmutations.R
\name{proposTumLoup}
\alias{proposTumLoup}
\title{Generate new spot-groups based on the number of mut and ref reads(the spots having
more mutant reads are more likely to have tumor cells/spots) for loupe visualization.}
\usage{
proposTumLoup(data1, data2)
}
\arguments{
\item{data1}{A dataframe of MutantSpotIndex_scored.csv, which is the output of spotSummary func.}

\item{data2}{A df of barcode, which is the first output of the sptBClstRds func.}
}
\value{
A df for Loupe visualization of AllSptTumPropsed.csv
}
\description{
Generate new spot-groups based on the number of mut and ref reads(the spots having
more mutant reads are more likely to have tumor cells/spots) for loupe visualization.
}
\examples{
\dontrun{
print("check user guider")

}
}
