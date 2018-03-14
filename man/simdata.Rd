\name{simdata}
\alias{simdata}
\docType{data}
\title{
Simulated dataset
}
\description{
The data were simulated from a linear mixed model. 
Repeated data of the longitudinal outcome were simulated for 500 subjects. 
Death time was simulated depending on the (observed and unobserved)
longitudinal outcome and on the binary covariate.
Missing data before death were simulated using a logistic regression
model including the binary covariate, the outcome at the previous visit and
the observation status at the previous visit.
}
\usage{simdata}
\format{
  A data frame with 2123 observations over 500 different subjects and 7 variables.
  \describe{
    \item{\code{id}}{subject identification number}
    \item{\code{X}}{binary covariate}
    \item{\code{death}}{death time (missing for subjects alive)}
    \item{\code{time}}{measurement time}
    \item{\code{age}}{age at measurement time}
    \item{\code{Y}}{longitudinal outcome}
    \item{\code{Ytrunc}}{longitudinal outcome truncated at the first
      missing value}
  }
}
\keyword{datasets}
