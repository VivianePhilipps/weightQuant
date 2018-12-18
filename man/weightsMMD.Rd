\name{weightsMMD}
\alias{weightsMMD}

\title{
Estimation of observation-specific weights with monotone missing
data
}
\description{
This function provides stabilized weights for incomplete longitudinal
data selected by death. The procedure assumes monotne missing data and a
MAR-S mechanism, that is the probability of being observed depends also on
further death.
Weights are defined as the inverse of the probability of being
observed. These are obtained by pooled logistic regressions.
}
\usage{
weightsMMD(data, Y, X1, X2, subject, death, time, interval.death = 0, name = "weight")
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
the first visit, 2 for the second visit, etc.
}
  \item{interval.death}{
integer vector, intervals (j-k) to consider for the MAR-S hypothesis (see
details). By default, interval.death=0, estimation under the MAR assumption.
}
  \item{name}{
character indicating the name of the weight variable that will be added
to the data
}
}
\details{
In longitudinal studies, follow-up can be truncated by death. Different
missingness mechanism can be assumed. Missing data can be :
 1. MCAR (completely at random) if the missingness probability is
 independent from the outcome and the death time
 2. MAR (missing at random ) if the probability
is independent from the unobserved values of the outcome and from the
death time
 3. MAR-S if the probability
is independent from the unobserved values but is different according to
the death time
 4. MNAR (missing not at random) if the probability may depend on
 unobserved values.

Denoting T_i the death time, R_{ij} the observation indicator for
subject i and occasion j, t the time, Y the outcome and X1 and X2  the
covariates, we propose weights for monotone missing data defined as :

w_{ij} = P(R_{ij} = 1 | T_i > t_{ij}, X1_{ij}, X2_{ij})  / P(R_{ij} = 1 | T_i >
 t_{ij}, X1_{ij}, X2_{ij}, Y_{ij-1}) 

 The numerator corresponds to the conditional probability of being
 observed in the population currently alive under the MCAR assumption.

 The denominator is computed as a telescoping product :
 
 P(R_{ij} = 1 | T_i > t_{ij}, X1_{ij}, X2_{ij}, Y_{ij-1}) =

 prod_{k=2}^j P(R_{ik} = 1 | R_{ik-1} =  1, T_i > t_{ij}, X1_{ij},
 X2_{ij}, Y_{ij-1}) = 

 prod_{k=2}^j lambda_{ijk}

The probability lambda_{ijk} are obtained by logistic regressions.


Under the MAR-S assumption, the regression model is :

logit(lambda_{ijk}) = b_{0k(j-k)} + b_{1(j-k)} X1_{ik} +  b_{2(j-k)}
Y_{i(k-1)} + b_{3(j-k)} X1_{ik} Y_{i(k-1)} +  b_{4(j-k)} X2_{ik} 

For each interval (j-k), one logistic regression is performed.


Under the MAR assumption, one logistic regression is performed :

logit(lambda_{ikk}) = b_{0k} + b_1 X1_{ik} + b_2 X2_{ik} +  b_3 Y_{i(k-1)} + b_4 X1_{ik} Y_{i(k-1)} 
 
}
\value{
  A list containing :
  \item{data}{the data frame with initial data and estimated weights as
    last column}
  \item{coef}{a list containing the estimates of the logistic
    regressions. The first element of coef contains the estimates under the MCAR
  assumption, the further contain the estimates under the MAR or MAR-S
  assumption.}
\item{se}{a list containing the standard erros of the estimates
  contained in coef, in the same order.}
}
\author{
Viviane Philipps, Marion Medeville, Anais Rouanet, Helene Jacqmin-Gadda
}

\seealso{
\code{\link{weightsIMD}}
}
\examples{
w_simdata <- weightsMMD(data=simdata,Y="Ytrunc",X1="X", X2=NULL,
subject="id", death="death", time="time", interval.death = 0)$data
}
