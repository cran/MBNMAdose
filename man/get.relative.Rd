% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inconsistency.functions.R
\name{get.relative}
\alias{get.relative}
\title{Calculates league table of effects between treatments in MBNMA and/or NMA models}
\usage{
get.relative(
  lower.diag,
  upper.diag = lower.diag,
  treatments = list(),
  lower.direction = "colvrow",
  upper.direction = "rowvcol",
  regress.vals = NULL,
  eform = FALSE,
  lim = "cred"
)
}
\arguments{
\item{lower.diag}{An S3 object either of class \code{"mbnma"} generated by running
a dose-response MBNMA model (using \code{mbnma.run()}), or of
class \code{"nma"} generated by running an NMA model (using \code{nma.run()}). Indicates that
this model should be used to estimate treatment effects in the \strong{lower-left diagonal} of the league table.}

\item{upper.diag}{Same as for \code{lower.diag}, but indicates that this model should be used
to estimate treatment effects in the \strong{upper-right diagonal} of the league table.}

\item{treatments}{A list whose elements each represent different treatments.
Treatment is defined as a combination of agent and dose. Only agents specified in
\code{lower.diag} and \code{upper.diag} can be included. Each element in \code{treatments} is named corresponding to the
agent and contains a numeric vector of doses. Relative effects will be calculated between
all treatments specified in \code{treatments}. If \code{treatments} is left empty then the set of
treatments in the network on which \code{lower.diag} is estimated will be used as the default.}

\item{lower.direction}{Whether treatment effects should be presented as the column versus the row treatment
for each cell in the \strong{lower-left diagonal} of the league table (\code{"colvrow"}) or as the row versus the
column treatment (\code{"rowvcol"}).}

\item{upper.direction}{Same as for \code{upper.direction} but for results shown in the \strong{upper-right diagonal} of
the league table.}

\item{regress.vals}{A named numeric vector of effect modifier values at which relative effects
should be estimated. Named elements must match variable names specified in regression design matrix
(\code{mbnma$model.arg$regress.mat}), corresponding to variables in \code{regress} within
the MBNMA model (or \strong{both} MBNMA models if different models are used for \code{upper.diag}
and \code{lower.diag}. Meta-regression is not currently possible for NMA models, so only
unadjusted results will be shown for these (and a warning given if they are compared to
adjusted MBNMA results).}

\item{eform}{Whether outputted results should be presented in their exponential form (e.g. for
models with log or logit link functions)}

\item{lim}{Specifies calculation of either 95\% credible intervals (\code{lim="cred"}) or 95\% prediction intervals (\code{lim="pred"}).}
}
\value{
An array of \verb{length(treatments) x length(treatments) x nsims}, where \code{nsims}
is the number of iterations monitored in \code{lower.diag}. The array contains the individual
MCMC values for each relative effect calculated between all \code{treatments} on the link scale
specified in the \code{lower.diag} and \code{upper.diag} models.
}
\description{
Calculates league table of effects between treatments in MBNMA and/or NMA models
}
\examples{
\donttest{

# Using the osteoarthritis data
network <- mbnma.network(osteopain)

# Run an MBNMA model
expon <- mbnma.run(network, fun=dexp(), method="random")

# Calculate relative effects for MBNMA between:
# Celebrex 100mg/d, Celebrex 200mg/d, Tramadol 100mg/d
rel.eff <- get.relative(lower.diag=expon,
  treatments=list("Celebrex"=c(100,200), "Tramadol"=100))

# Run an NMA model
nma <- nma.run(network, method="random")

# Compare results between MBNMA and NMA models
rel.eff <- get.relative(lower.diag=expon, upper.diag=nma,
  treatments=list("Celebrex"=c(100,200), "Tramadol"=100),
  upper.direction="colvrow")
}

}
