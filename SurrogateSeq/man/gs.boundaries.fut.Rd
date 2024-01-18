\name{gs.boundaries.fut}
\alias{gs.boundaries.fut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes group sequential boundaries with futility stopping
}
\description{
Computes group sequential (and naive) boundaries for the nonparametric test for a treatment effect on the primary outcome using surrogate marker information. The boundaries and test statistic borrow information from a prior study (Study A) about the relationship between the surrogate and the primary outcome to test for a treatment effect in the current study (Study B). The group sequential boundaries allow for futility stopping (bounds given).}
\usage{
gs.boundaries.fut(szerop, sonep, yzerop, n.stg, B.norm = 1e+06, alpha = 0.05, 
pp = 0.4, inf.fraction = (1:n.stg)/n.stg, j.star=1, alpha0=(j.star/n.stg)*alpha, 
plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{szerop}{
surrogate marker in the control group in Study A
}
  \item{sonep}{
surrogate marker in the treated group in Study A
}
  \item{yzerop}{
primary outcome in the control group in Study A
}
  \item{n.stg}{
maximum number of analyses
}
  \item{B.norm}{
number of multivariate normal vectors to use in simulation for boundaries; default is 1e+06
}
  \item{alpha}{
desired rejection probability of the test; default is 0.05
}
\item{pp}{
	power parameter for Wang-Tsiatis boundaries; default is 0.4
}
\item{inf.fraction}{
	information fraction vector of the same length as n.stg which reflects the fraction of information at each analysis, should be positive, non-decreasing, and the last entry should be 1; default is (1:n.stg)/n.stg, user may want to specify a different vector for unequal time points
}
\item{j.star}{
	 earliest stage at which futility stopping is allowed. Should be <= n.stg-1 (there is already "futility stopping" at the n.stg-th stage anyway). Default is 1.
}
\item{alpha0}{
	the part of alpha that c1 is chosen to spend in first j.star stages; default is (j.star/n.stg)*alpha
}
\item{plot}{
	TRUE or FALSE if a plot of the boundaries is desired; default is FALSE
}
}
\value{
Returns a list of boundaries:
 \item{Naive}{Naive boundaries }
 \item{Bonf}{Bonferroni boundaries}
 \item{Pocock.futility}{Pocock futility boundaries}
 \item{Pocock.nullrejection}{Pocock null rejection boundaries}
 \item{OBrien_Fleming.futility}{O'Brien-Fleming futility boundaries}
  \item{OBrien_Fleming.nullrejection}{O'Brien-Fleming null rejection boundaries}
 \item{Wang_Tsiatis.futility}{Wang-Tsiatis futility boundaries}
  \item{Wang_Tsiatis.nullrejection}{Wang-Tsiatis null rejection boundaries}
}

\references{
Parast and Bartroff (2023+). Group Sequential Testing of a Treatment Effect Using a Surrogate Marker. Under Review.
}
\author{
Layla Parast and Jay Bartroff
}
\examples{
data(example.data)
data(StudyA.aids)
data(StudyB.aids)
s0.studya = StudyA.aids$s0
s1.studya = StudyA.aids$s1

bound = gs.boundaries.fut(szerop = s0.studya, sonep = s1.studya, 
yzerop=StudyA.aids$y0, n.stg=4, B.norm=1e6, alpha=0.05)

bound
}

\keyword{nonparametric}
\keyword{test}
\keyword{boundaries}
\keyword{futility}