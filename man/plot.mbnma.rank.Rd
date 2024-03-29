% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbnma.rank-class.R
\name{plot.mbnma.rank}
\alias{plot.mbnma.rank}
\title{Plot histograms of rankings from MBNMA models}
\usage{
\method{plot}{mbnma.rank}(x, params = NULL, treat.labs = NULL, ...)
}
\arguments{
\item{x}{An object of class \code{"mbnma.rank"} generated by \code{rank.mbnma()}}

\item{params}{A character vector of named parameters in the model that vary by either agent
or class (depending on the value assigned to \code{level}). If left as \code{NULL} (the default), then
ranking will be calculated for all available parameters that vary by agent/class.}

\item{treat.labs}{A vector of treatment labels in the same order as treatment codes.
Easiest to use treatment labels stored by \code{mbnma.network()}}

\item{...}{Arguments to be sent to \code{ggplot::geom_bar()}}
}
\value{
A series of histograms that show rankings for each treatment/agent/prediction, with a
separate panel for each parameter.
The object returned is a list containing a separate element for each parameter in \code{params}
which is an object of \code{class(c("gg", "ggplot"))}.
}
\description{
Plot histograms of rankings from MBNMA models
}
\examples{
\donttest{
# Using the triptans data
network <- mbnma.network(triptans)

# Estimate rankings  from an Emax dose-response MBNMA
emax <- mbnma.run(network, fun=demax(), method="random")
ranks <- rank(emax)

# Plot rankings for both dose-response parameters (in two separate plots)
plot(ranks)

# Plot rankings just for ED50
plot(ranks, params="ed50")

# Plot rankings from prediction
doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
            exact.doses=doses)
rank <- rank(pred)
plot(rank)


# Trying to plot a parameter that has not been ranked will return an error
#### ERROR ####
# plot(ranks, params="not.a.parameter")
}
}
