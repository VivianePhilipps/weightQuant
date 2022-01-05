\name{weightQuant-package}
\alias{weightQuant-package}
\alias{weightQuant}
\docType{package}
\title{
Weights for incomplete longitudinal data and quantile regression
}
\description{
Functions for the estimation of observation-specific weights for
incomplete longitudinal data. A bootstrap method is also provided to
obtain standard erros of weighted quantile regressions.
}
\details{
  \Sexpr[stage=build,results=hide]{lb <- library(help="weightQuant",
    character.only=TRUE)}
  
  \Sexpr[stage=build,results=hide]{descr <-
    packageDescription("weightQuant")}
  
  \tabular{ll}{
    Package:  \tab \Sexpr[stage=build]{descr$Package}\cr
    Type:     \tab \Sexpr[stage=build]{descr$Type}\cr
    Version:  \tab \Sexpr[stage=build]{descr$Version} \cr
    Date:     \tab \Sexpr[stage=build]{descr$Date}\cr
    License:  \tab \Sexpr[stage=build]{descr$License}\cr
  }

  Index:
  
  \Sexpr[stage=build,results=rd]{paste("\\\\preformatted{", paste(if(!is.null(lb$info[[2]])) lb$info[[2]] else "No entries",collapse="\\n"), "}", sep="")}
}
\references{
Jacqmin-Gadda H, Rouanet A, Mba RD, Philipps V, Dartigues J-F. Quantile regression for incomplete longitudinal data with selection  by death. Statistical Methods in Medical Research. 2020;29(9):2697-2716. doi:10.1177/0962280220909986
}
\author{
Viviane Philipps
}
\keyword{ package }
