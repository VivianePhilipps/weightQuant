\name{bootwrq}
\alias{bootwrq}

\title{
Bootstrap procedure for weighted quantile regressions
}
\description{
A subject-level bootstrap method for weighted quantile regressions is
implemented in this function. Quantile regressions are estimated in a
generalized estimating equation framework with independent working
covariance matrix. Weights are estimated using \code{weightsIMD} or
\code{weightsMMD} functions.
}
\usage{
bootwrq(B, form, tau, data, Y, X1 = NULL, X2 = NULL, subject,
death, time, interval.death = NULL, impute = NULL, weight = NULL,
wcompute = 2, seed = NULL, intermittent, file = NULL,
nproc = 1, MPI = FALSE)
}

\arguments{
  \item{B}{
integer, number of bootstrap samples
}
  \item{form}{
formula indicating the quantile regression model to be estimated
}
  \item{tau}{
numeric vector indicating the quantiles to be estimated
}
  \item{data}{
data frame containing the data
}
  \item{Y}{
character indicating the name of the response outcome
}
  \item{X1}{
optional character vector passed to the weight functions
}
  \item{X2}{
optional character vector passed to the weight functions
}
  \item{subject}{
character indicating the name of the subject identifier
}
  \item{death}{
optional character passed to the weight functions
}
  \item{time}{
optional character passed to the weight functions
}
  \item{interval.death}{
optional numeric vector passed to the weight function weightsMMD
}
  \item{impute}{
optional numeric vector passed to the weight function weightsIMD
}
  \item{weight}{
character indicating the name of the weight variable in data
}
  \item{wcompute}{
integer indicating if weights should be estimated in each bootstrap
  sample. If wcompute=0, weights are supposed to be known. If
  wcompute=1, weights are re-estimated in each bootstrp sample. If
  wcompute=2, both results are returned.
}
  \item{seed}{
optional integer vector of length B indicating the seeds.
}
  \item{intermittent}{
logical indicating if data contains intermittent missing data. If
  intermittent=TRUE, the weights are estimated using weightsIMD
  function, if intermittent=FALSE, the weights are estimated using weightsMMD
  function.
}
  \item{file}{
optional character indicating the name of the results file. If
  file=NULL, no results file is created.
}
  \item{nproc}{
number of processors to be used for parallel computing. Default to 1,
  sequential computation.
}
  \item{MPI}{
logical indicating if MPI parallelization should be used. Default to FALSE.
}
}
\value{
a matrix with B columns containing the results on each
bootstrap sample.
}
\author{
Viviane Philipps, Robert Darlin Mba
}

\seealso{
 \code{\link{summary.bootwrq}},  \code{\link{test.bootwrq}}
}
\examples{
\dontrun{
## computation of the weights with intermittent missing data 
w_simdata <- weightsIMD(data=simdata,Y="Y",X1="X",X2=NULL,subject="id",
death="death",time="time",impute=20,name="w_imd")$data

## estimation of the weighted quantile regressions
## for the first quartile and the median
m_simdata <- rq(Y~time*X,data=w_simdata,weights=w_imd,tau=c(0.25,0.5))

## estimation of the standard erros using the bootstrap procedure
boot_simdata <- bootwrq(B=1000, form=Y~time*X, tau=c(0.25,0.5),
data=w_simdata, Y="Y",X1="X",X2=NULL,subject="id",
death="death",time="time",impute=20,wcompute=0,intermittent=TRUE)

## the summary of the results
summary(boot_simdata,m_simdata)

## comparison of the covariate effects
## between the first quartile and the median
test.bootwrq(boot_simdata,m_simdata)
}
}
