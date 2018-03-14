\name{test.bootwrq}
\alias{test.bootwrq}
\title{
Test of covariate effects between different quantiles
}
\description{
This function provides a test for the covariate effects estimated for
different quantiles.
}
\usage{
test.bootwrq(x, m)
}
\arguments{
  \item{x}{
results from bootstrap estimations obtained with bootwrq function
}
  \item{m}{
a quantile regression model estimated with rq function from quantreg
package. At least 2 quantiles should be specified in rq function.
}
}
\details{
For 2 quantiles tau1 and tau2, the test of the null hypothesis H0 : b_{tau1} =
b_{tau2} is obtained with the following procedure
:
1. estimate the difference diff = b_{tau1} - b_{tau2} on the initial sample
(ie from model m)
2. estimate the difference diff_b = b_{tau1}^b - b_{tau2}^b on each of the B
bootstrap samples
3. compute se_diff, the empirical standard error of these B differences
4. the associated p-value is obtained with the Gaussian assumption (
p-value = 2*P(N(0,1) > abs(diff/se_diff)) )
}
\value{
  A list containing :
  \item{results0}{a matrix with 3 columns containing the results
  (difference of the coefficients, standard erros of the diffrence and
  assocated p-values) without computing the
  weights in each bootstrap sample. Or NULL if the bootstrap results are
  obtained with wcompute=1.}
  \item{results1}{a matrix with 3 columns containing the results
  (difference of the coefficients, standard erros of the diffrence and
  assocated p-values) with re-estimated weights
  on each bootstrap sample. Or NULL if the bootstrap results are
  obtained with wcompute=0.}
}
\author{
Viviane Philipps
}
