\name{delta.e.estimate}
\alias{delta.e.estimate}
\title{
Tests for a treatment effect on the primary outcome using surrogate marker information}
\description{
Nonparametric test for a treatment effect on the primary outcome using surrogate marker information. This test borrows information from a prior study (Study A) about the relationship between the surrogate and the primary outcome to test for a treatment effect in the current study (Study B).}
\usage{
delta.e.estimate(sone = NULL, szero = NULL, szerop, yzerop, extrapolate = TRUE, 
mat = NULL, n1 = NULL, n0 = NULL)
}
\arguments{
  \item{sone}{
surrogate marker in the treated group in Study B}
  \item{szero}{
surrogate marker in the control group in Study B
}
  \item{szerop}{
surrogate marker in the control group in Study A
}
  \item{yzerop}{
primary outcome in the control group in Study A
}
  \item{extrapolate}{
TRUE or FALSE; extrapolate for values outside of the support in Study A
}
  \item{mat}{
for Study B, the user can either provide sone and szero or can provide a vector, mat, where the first n1 values are the surrogate marker in the treated group in the Study B, and the remaining values are the surrogate marker in the control group in Study B
}
  \item{n1}{
sample size of treated group in Study B; only needed if mat is provided instead of sone and szero
}
  \item{n0}{
sample size of control group in Study B; only needed if mat is provided instead of sone and szero
}
}
\value{
\item{delta.e}{estimated treatment effect using surrogate marker information}
\item{sd.closed}{estimated standard error of treatment effect estimate}
\item{delta.e.z}{test statistic}
\item{delta.e.p}{p-value of test statistic} 
}
\references{
Parast, Cai, and Tian (2023). Using a Surrogate with Heterogeneous Utility to Test for a Treatment Effect. Statistics in Medicine, 42(1): 68-88.

Parast and Bartroff (2024+). Group Sequential Testing of a Treatment Effect Using a Surrogate Marker. Under Review.
}
\author{
Layla Parast
}
\examples{
data(example.data)
delta.e.estimate(sone = example.data$s1, szero = example.data$s0, szerop = example.data$s0.p, 
yzerop = example.data$y0.p)

data(StudyA.aids)
data(StudyB.aids)
s1.studyb = StudyB.aids$s1
s0.studyb = StudyB.aids$s0
s0.studya = StudyA.aids$s0

#24 weeks

delta.e.vec = delta.e.estimate(sone=s1.studyb$CD4_24weeks[!is.na(s1.studyb$CD4_24weeks)],
szero=s0.studyb$CD4_24weeks[!is.na(s0.studyb$CD4_24weeks)], szerop = s0.studya$CD4_24weeks, 
yzerop = StudyA.aids$y0, extrapolate = TRUE)
delta.e.vec
}
\keyword{nonparametric}
\keyword{test}