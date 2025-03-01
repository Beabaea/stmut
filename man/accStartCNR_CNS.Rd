% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNV_Plot.R
\name{accStartCNR_CNS}
\alias{accStartCNR_CNS}
\title{Get accumulated start position for cnr and cns file for arm-level CNV plotting.}
\usage{
accStartCNR_CNS(cnr, seg)
}
\arguments{
\item{cnr}{the dataframe of cnr file from cnvkit output}

\item{seg}{the dataframe of cns file from cnvkit output}
}
\value{
A list of 2 dataframe, the new_ch for cnr file; and the newSeg for cns
}
\description{
Get accumulated start position for cnr and cns file for arm-level CNV plotting.
}
\examples{
\dontrun{
cnr <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cnr",
package = "stmut"), sep="\t", header = TRUE)
cns <- read.table(system.file("extdata/", "P6_T.deduplicated.realign.cns",
package = "stmut"), sep="\t", header = TRUE)
df9 <- accStartCNR_CNS(cnr=cnr, seg=cns)
}
}
