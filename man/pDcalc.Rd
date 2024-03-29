% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.functions.R
\name{pDcalc}
\alias{pDcalc}
\title{Calculate plugin pD from a JAGS model with univariate likelihood for studies
with repeated measurements}
\usage{
pDcalc(
  obs1,
  obs2,
  fups = NULL,
  narm,
  NS,
  theta.result,
  resdev.result,
  likelihood = "normal",
  type = "time"
)
}
\arguments{
\item{obs1}{A matrix (study x arm) or array (study x arm x time point) containing
observed data for \code{y} (normal likelihood) or \code{r} (binomial or poisson likelihood)
in each arm of each study. This will be the same array
used as data for the JAGS model.}

\item{obs2}{A matrix (study x arm) or array (study x arm x time point) containing
observed data for \code{se} (normal likelihood), \code{n} (binomial likelihood) or \code{E} (poisson likelihood)
in each arm of each study. This will be the same array
used as data for the JAGS model.}

\item{fups}{A numeric vector of length equal to the number of studies,
containing the number of follow-up mean responses reported in each study. Required for
time-course MBNMA models (if \code{type="time"})}

\item{narm}{A numeric vector of length equal to the number of studies,
containing the number of arms in each study.}

\item{NS}{A single number equal to the number of studies in the dataset.}

\item{theta.result}{A matrix (study x arm) or array (study x arm x time point)
containing the posterior mean predicted means/probabilities/rate in each arm of each
study. This will be estimated by the JAGS model.}

\item{resdev.result}{A matrix (study x arm) or array (study x arm x time point)
containing the posterior mean residual deviance contributions in each arm of each
study. This will be estimated by the JAGS model.}

\item{likelihood}{A character object of any of the following likelihoods:
\itemize{
\item \code{normal}
\item \code{binomial} (does not work with time-course MBNMA models)
\item \code{poisson} (does not work with time-course MBNMA models)
}}

\item{type}{The type of MBNMA model fitted. Can be either \code{"time"} or \code{"dose"}}
}
\value{
A single numeric value for pD calculated via the plugin method.
}
\description{
Uses results from MBNMA JAGS models to calculate pD via the
plugin method \insertCite{spiegelhalter2002}{MBNMAdose}. Can only be used for models with known
standard errors or covariance matrices. Currently only functions with univariate likelihoods. Function
is identical in MBNMAdose and MBNMAtime packages.
}
\details{
Method for calculating pD via the plugin method proposed by
Spiegelhalter \insertCite{spiegelhalter2002}{MBNMAdose}. Standard errors / covariance matrices must be assumed
to be known. To obtain values for \code{theta.result} and \code{resdev.result} these
parameters must be monitored when running the MBNMA model (using \code{parameters.to.save}).

For non-linear time-course MBNMA models residual deviance contributions may be skewed, which
can lead to non-sensical results when calculating pD via the plugin method.
Alternative approaches are to use pV as an approximation or
pD calculated by Kullback-Leibler divergence \insertCite{plummer2008}{MBNMAdose}.
}
\examples{
\donttest{
# Using the triptans data
network <- mbnma.network(triptans)

# Fit a dose-response MBNMA, monitoring "psi" and "resdev"
result <- mbnma.run(network, fun=dloglin(), method="random",
              parameters.to.save=c("psi", "resdev"))


#### Calculate pD for binomial data ####

# Prepare data for pD calculation
r <- result$model$data()$r
n <- result$model$data()$n
narm <- result$model$data()$narm
NS <- result$model$data()$NS

psi <- result$BUGSoutput$median$psi
resdevs <- result$BUGSoutput$median$resdev

# Calculate pD via plugin method
pD <- pDcalc(obs1=r, obs2=n, narm=narm, NS=NS,
          theta.result=psi, resdev.result=resdevs,
          likelihood="binomial", type="dose")
}

}
\references{
\insertAllCited
}
