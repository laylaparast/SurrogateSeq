\name{StudyB.aids}
\alias{StudyB.aids}
\docType{data}
\title{
ACTG 193A clinical trial data
}
\description{
Surrogate marker measurements over time from the ACTG 193A clinical trial data. Note that the time points do not exactly match up to ACTG 320. In the paper, we use Study A surrogate data at 24 weeks to construct the conditional mean function applied to Study B at 16 weeks. Also note that some subjects are missing values of the surrogate at one or more time points. The naive estimate of the treatment effect using the surrogates uses all non-missing data available at each time point.
}
\usage{data("StudyB.aids")}
\format{
 A list with 2 elements:
  \describe{
  	   \item{\code{s1}}{a dataframe of the surrogate markers at different time points in the treatment group in Study B; the surrogate marker is change in CD4 cell count from baseline to 8 weeks (CD4_8weeks), 16 weeks (CD4_16weeks), 24 weeks (CD4_24weeks), and 40 weeks (CD4_40weeks). Note that higher values indicate increasing CD4 cell count which is "better".}
    \item{\code{s0}}{a dataframe of the surrogate markers at different time points in the control group in Study B}
}
}
\examples{
data(StudyB.aids)
}
\keyword{datasets}
