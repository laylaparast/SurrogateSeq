\name{cov.surr.gs}
\alias{cov.surr.gs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes variances and standardized covariance matrix for the group sequential statistic
}
\description{
Computes variances and standardized covariance matrix for the group sequential statistic
}
\usage{
cov.surr.gs(s0.4.est, s1.4.est, sa.0, ya.0, nb.0, nb.1, full.matrix = TRUE, 
naive = FALSE)
}
\arguments{
  \item{s0.4.est}{
surrogate marker in the control group which is used for estimating means and covariances of S0, S1 in the Study B data.  For designing tests (e.g., finding boundaries) these may come from Study A data, but for analyzing tests these may come from Study B data. Number of columns is the number of stages, number of rows may differ from rows in sa.0}
  \item{s1.4.est}{
surrogate marker in the treated group which is used for estimating means and covariances of S0, S1 in the Study B data.  For designing tests (e.g., finding boundaries) these may come from Study A data, but for analyzing tests these may come from Study B data. Number of columns is the number of stages, number of rows may differ from rows in sa.0}
  \item{sa.0}{
surrogate marker in the control group in Study A
}
  \item{ya.0}{
primary outcome in the control group in Study A
}
  \item{nb.0}{
sample size for the control group in Study B
}
  \item{nb.1}{
sample size for the treated group in Study B
}
  \item{full.matrix}{
if TRUE, the standardized covariance matrix is provided; default is TRUE
}
  \item{naive}{
user should set to TRUE to compute covariance for "cumulative" test statistic, FALSE for naive statistic that only uses study B data from timepoint J at the J-th analysis; default is FALSE
}
}

\value{
Returns a list:
 \item{var.vec.del }{variance vector computed by the delta method}
  \item{cov.stand.del}{if full.matrix = TRUE, covariance matrix of the standardized test statistic computed by the delta method}
    \item{var.vec.samp}{variance vector computed by the sample mean and covariance of s0.4.est and s1.4.est}
      \item{cov.stand.samp}{if full.matrix = TRUE, covariance matrix of the standardized test statistic computed by the sample mean and covariance of s0.4.est and s1.4.est}
}

\author{
Jay Bartroff
}
\keyword{internal}
