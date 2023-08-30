\name{example.data}
\alias{example.data}
\docType{data}
\title{
Example data
}
\description{
Example data
}
\usage{data("example.data")}
\format{
 A list with 9 elements:
  \describe{
  	  \item{\code{w0.p}}{the baseline covariate in the control group in the prior study (Study A)}
    \item{\code{s0.p}}{the surrogate marker in the control group in the prior study (Study A}
    \item{\code{y0.p}}{the primary outcome in the control group in the prior study (Study A}
  	 \item{\code{w1}}{a baseline covariate in the treatment group in the current study (Study B)}
    \item{\code{w0}}{a baseline covariate in the control group in the current study (Study B)}
     \item{\code{s1}}{the surrogate marker in the treatment group in the current study (Study B)}
    \item{\code{s0}}{the surrogate marker in the control group in the current study (Study B)}
      \item{\code{y1}}{the primary outcome in the treatment group in the current study (Study B)}
    \item{\code{y0}}{the primary outcome in the control group in the current study (Study B)}
}
}
\examples{
data(example.data)
names(example.data)
}

