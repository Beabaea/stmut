% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNVsScorPermutation.R
\name{fdrTumor}
\alias{fdrTumor}
\title{Generate tumor groups based on FDR rate for Loupe visualization.}
\usage{
fdrTumor(a = 0, b = 1.82, c = 4.4, data1, data2, data3, data4)
}
\arguments{
\item{a}{CNVscor cutoff less than a for Non_tumor}

\item{b}{CNVscor cutoff for potential_tumor and Tossup ( b=< Potential_tumor_CNVscor < c;
a =< Tossup_CNVscor_cutoff < b)}

\item{c}{CNVscor cutoff greater than c for Tumor spots}

\item{data1}{Sample/case CNVscore summary, the output of permutScore function}

\item{data2}{Sample grouped spots representative summary}

\item{data3}{Grouped_spots summary, which spots into the same group, which
should be given the same color for loupe visualization}

\item{data4}{All spots summary}
}
\value{
A dataframe of FDR-based tumor groups, used for Loupe visualization
}
\description{
Generate tumor groups based on FDR rate for Loupe visualization.
}
\examples{
\dontrun{
print("check user guider")
}
}
