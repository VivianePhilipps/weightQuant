\name{weightsIMD}
\alias{weightsIMD}
\title{
Estimation of observation-specific weights with intermittent missing
data
}
\description{
This function provides stabilized weights for incomplete longitudinal
data selected by death. The procedure allows intermittent missing data
and assumes a missing at random (MAR) mechanism.
Weights are defined as the inverse of the probability of being
observed. These are obtained by pooled logistic regressions.
}
\usage{
weightsIMD(data, Y, X1, X2, subject, death, time, impute = 0, name = "weight")
}
\arguments{
  \item{data}{
data frame containing the observations and all variables named in
\code{Y}, \code{X1}, \code{X2}, \code{subject}, \code{death} and
\code{time} arguments.
}
  \item{Y}{
character indicating the name of the response outcome
}
  \item{X1}{
character vector indicating the name of the covariates with 
interaction with the outcome Y in the logistic regressions
}
  \item{X2}{
character vector indicating the name of the covariates without 
interaction with the outcome Y in the logistic regressions
}
  \item{subject}{
character indicating the name of the subject identifier
}
  \item{death}{
character indicating the time of death variable 
}
  \item{time}{
character indicating the measurement time variable. Time should be 1 for
the first (theoretical) visit, 2 for the second (theoretical) visit, etc.
}
  \item{impute}{
numeric indicating the value to impute if the outcome Y is missing
}
  \item{name}{
character indicating the name of the weight variable that will be added
to the data
}
}
\details{
Denoting T_i the death time, R_{ij} the observation indicator for
subject i and occasion j, t the time, Y the outcome and X1 and X2 the covariates, we propose
weights for intermittent missing data defined as :


 w_{ij} = P(R_{ij} = 1 | T_i > t_{ij}, X1_{ij}, X2_{ij})  / P(R_{ij} = 1 | T_i >
 t_{ij}, X1_{ij}, X2_{ij}, Y_{ij-1}) 

 The numerator corresponds to the conditional probability of being
 observed in the population currently alive under the MCAR assumption.

 The denominator is computed by recurrence :
 
 P(R_{ij} = 1 | T_i > t_{ij}, X1_{ij}, X2_{ij}, Y_{ij-1}) =

 P(R_{ij} = 1 | T_i > t_{ij-1}, X1_{ij}, X2_{ij}, Y_{ij-1}, R_{ij-1} = 0) *
 P(R_{ij-1} = 0 | T_i > t_{ij}, X1_{ij}, X2_{ij}, Y_{ij-1}) + P(R_{ij} =
 1 | T_i > t_{ij-1}, X1_{ij}, X2_{ij}, Y_{ij-1}, R_{ij-1} = 1) *
 P(R_{ij-1} = 1 | T_i > t_{ij}, X1_{ij}, X2_{ij}, Y_{ij-1})

 Under the MAR assumption, the conditional probabilities lambda_{ij} = P(R_{ij} =
 1 | T_i > t_{ij}, X1_{ij}, X2_{ij}, Y_{ij-1}, R_{ij-1}) are obtained
 from the logistic regression :

 logit(lambda_{ij}) = b_{0j} + b_1 X1_{ij} +  b_2 X2_{ij} +  b_3 Y_{i(j-1)} + b_4
 X1_{ij} Y_{i(j-1)} + b_5 (1-R_{ij-1})
}
\value{
  A list containing :
  \item{data}{the data frame with initial data and estimated weights as
    last column}
  \item{coef}{a list containing the estimates of the logistic
    regressions. The first element of coef contains the estimates under the MCAR
  assumption, the second contains the estimates under the MAR
  assumption.}
\item{se}{a list containing the standard erros of the estimates
  contained in coef, in the same order.}
}
\author{
Viviane Philipps, Marion Medeville, Anais Rouanet, Helene Jacqmin-Gadda
}

\seealso{
\code{\link{weightsMMD}}
}
\examples{
w_simdata <- weightsIMD(data=simdata,Y="Y",X1="X",X2=NULL,subject="id",
death="death",time="time",impute=20,name="w_imd")$data
}
