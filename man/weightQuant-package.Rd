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
\author{
Viviane Philipps
}
\keyword{ package }
