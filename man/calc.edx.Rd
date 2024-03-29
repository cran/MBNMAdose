% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.functions.R
\name{calc.edx}
\alias{calc.edx}
\title{Calculates values for EDx from an Emax model, the dose at which x\% of the maximal response (Emax)
is reached}
\usage{
calc.edx(mbnma, x = 50)
}
\arguments{
\item{mbnma}{An S3 object of class \code{"mbnma"} generated by running
a dose-response MBNMA model}

\item{x}{A numeric value between 0 and 100 for the dose at which x\% of the maximal response (Emax)
should be calculated}
}
\value{
A data frame of posterior EDx summary values for each agent
}
\description{
Calculates values for EDx from an Emax model, the dose at which x\% of the maximal response (Emax)
is reached
}
