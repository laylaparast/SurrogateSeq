\name{StudyA.aids}
\alias{StudyA.aids}
\docType{data}
\title{
ACTG 320 clinical trial data
}
\description{
Primary outcome and surrogate marker measurements over time from the ACTG 320 clinical trial data
}
\usage{data("StudyA.aids")}
\format{
 A list with 4 elements:
  \describe{
  	  \item{\code{y1}}{the primary outcome in the treatment group in Study A; the primary outcome is defined as -1 times (log of RNA at 40 weeks - log of RNA at baseline) because a DECREASE in RNA is better}
    \item{\code{y0}}{the primary outcome in the control group in Study A}
    \item{\code{s1}}{a dataframe of the surrogate markers at different time points in the treatment group in Study A; the surrogate marker is change in CD4 cell count from baseline to 4 weeks (CD4_4weeks), 8 weeks (CD4_8weeks), 24 weeks (CD4_24weeks), and 40 weeks (CD4_40weeks). Note that higher values indicate increasing CD4 cell count which is "better".}
    \item{\code{s0}}{a dataframe of the surrogate markers at different time points in the control group in Study A}
}
}

\examples{
data(StudyA.aids)
}
\keyword{datasets}
