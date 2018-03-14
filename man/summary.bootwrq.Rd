\name{summary.bootwrq}
\alias{summary.bootwrq}

\title{
Summary of a quantile regression model
}
\description{
The function provides a summary of quantile regression
estimation. Standard erros and p values are obtained from a bootstrap procedure.
}
\usage{
\method{summary}{bootwrq}(object, \dots)
}

\arguments{
  \item{object}{
    results from bootstrap estimations obtained with bootwrq function
  }
  \item{...}{additional arguments. If a quantile regression model
  estimated with rq function from quantreg package is specified, the
  function uses these estimated coefficients as results. Otherwise, the
  coefficients are obtained as the mean over the B estimated
  coefficients from the bootstrap results.
}
}
\value{
  A list containing :
  \item{results0}{a matrix with 3 columns containing the results
  (coefficients, standard erros and p-values) without computing the
  weights in each bootstrap sample. Or NULL if the bootstrap results are
  obtained with wcompute=1.}
  \item{results1}{a matrix with 3 columns containing the results
  (coefficients, standard erros and p-values) with re-estimated weights
  on each bootstrap sample. Or NULL if the bootstrap results are
  obtained with wcompute=0.}
}
\author{
Viviane Philipps
}
