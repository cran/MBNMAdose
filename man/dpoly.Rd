% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dose.functions.R
\name{dpoly}
\alias{dpoly}
\title{Polynomial dose-response function}
\usage{
dpoly(
  degree = 1,
  beta.1 = "rel",
  beta.2 = "rel",
  beta.3 = "rel",
  beta.4 = "rel"
)
}
\arguments{
\item{degree}{The degree of the polynomial - e.g. \code{degree=1} for linear, \code{degree=2} for quadratic, \code{degree=3} for cubic.}

\item{beta.1}{Pooling for the 1st polynomial coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}

\item{beta.2}{Pooling for the 2nd polynomial coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}

\item{beta.3}{Pooling for the 3rd polynomial coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}

\item{beta.4}{Pooling for the 4th polynomial coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}
}
\value{
An object of \code{class("dosefun")}
}
\description{
Polynomial dose-response function
}
\details{
\itemize{
\item \eqn{\beta_1} represents the 1st coefficient.
\item \eqn{\beta_2} represents the 2nd coefficient.
\item \eqn{\beta_3} represents the 3rd coefficient.
\item \eqn{\beta_4} represents the 4th coefficient.
}

Linear model:
\deqn{\beta_1{x}}

Quadratic model:
\deqn{\beta_1{x} + \beta_2{x^2}}

Cubic model:
\deqn{\beta_1{x} + \beta_2{x^2} + \beta_3{x^3}}

Quartic model:
\deqn{\beta_1{x} + \beta_2{x^2} + \beta_3{x^3} + \beta_4{x^4}}
}
\section{Dose-response parameters}{
\tabular{ll}{
   \strong{Argument} \tab \strong{Model specification} \cr
   \code{"rel"} \tab Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. \cr
   \code{"common"} \tab Implies that all agents share the same common effect for this dose-response parameter. \cr
   \code{"random"} \tab Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. \cr
   \code{numeric()} \tab Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. \cr
}


When relative effects are modelled on more than one dose-response parameter,
correlation between them is automatically estimated using a vague inverse-Wishart prior.
This prior can be made slightly more informative by specifying the scale matrix \code{omega}
and by changing the degrees of freedom of the inverse-Wishart prior
using the \code{priors} argument in \code{mbnma.run()}.
}

\examples{
# Linear model with random effects
dpoly(beta.1="rel")

# Quadratic model dose-response function
# with an exchangeable (random) absolute parameter estimated for the 2nd coefficient
dpoly(beta.1="rel", beta.2="random")

}
\references{
\insertAllCited
}
