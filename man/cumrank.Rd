% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.functions.R
\name{cumrank}
\alias{cumrank}
\title{Plot cumulative ranking curves from MBNMA models}
\usage{
cumrank(x, params = NULL, sucra = TRUE, ...)
}
\arguments{
\item{x}{An object of class \code{"mbnma.rank"} generated by \code{rank.mbnma()}}

\item{params}{A character vector of named parameters in the model that vary by either agent
or class (depending on the value assigned to \code{level}). If left as \code{NULL} (the default), then
ranking will be calculated for all available parameters that vary by agent/class.}

\item{sucra}{A logical object to indicate whether Surface Under Cumulative Ranking Curve (SUCRA)
values should be calculated and returned as a data frame. Areas calculated
using trapezoid approach.}

\item{...}{Arguments to be sent to \code{ggplot::geom_line()}}
}
\value{
Line plots showing the cumulative ranking probabilities for each agent/class and
dose-response parameter in \code{x}. The object returned is a list which contains the plot
(an object of \verb{class(c("gg", "ggplot")}) and a data frame of SUCRA values
if \code{sucra = TRUE}.
}
\description{
Plot cumulative ranking curves from MBNMA models
}
\examples{
\donttest{
# Using the triptans data
network <- mbnma.network(triptans)

# Estimate rankings  from an Emax dose-response MBNMA
emax <- mbnma.run(network, fun=demax(), method="random")
ranks <- rank(emax)

# Plot cumulative rankings for both dose-response parameters simultaneously
# Note that SUCRA values are also returned
cumrank(ranks)
}
}
