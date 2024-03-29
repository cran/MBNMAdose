% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare.functions.R
\name{mbnma.validate.data}
\alias{mbnma.validate.data}
\title{Validates that a dataset fulfills requirements for MBNMA}
\usage{
mbnma.validate.data(data.ab, single.arm = FALSE)
}
\arguments{
\item{data.ab}{A data frame of arm-level data in "long" format containing the columns:
\itemize{
\item \code{studyID} Study identifiers
\item \code{dose} Numeric data indicating the dose (must take positive values)
\item \code{agent} Agent identifiers (can be numeric, factor or character)
\item \code{y} Numeric data indicating the aggregate response for a continuous outcome. Required for
continuous data.
\item \code{se} Numeric data indicating the standard error for a given observation. Required for
continuous data.
\item \code{r} Numeric data indicating the number of responders within a study arm. Required for
binomial or poisson data.
\item \code{n} Numeric data indicating the total number of participants within a study arm. Required for
binomial data or when modelling Standardised Mean Differences
\item \code{E} Numeric data indicating the total exposure time for participants within a study arm. Required
for poisson data.
\item \code{class} An optional column indicating a particular class code. Agents with the same identifier
must also have the same class code.
\item \code{standsd} An optional column of numeric data indicating reference SDs used to standardise
treatment effects when modelling using Standardised Mean Differences (SMD).
}}

\item{single.arm}{A boolean object to indicate whether to allow single arm studies in the dataset (\code{TRUE})
or not (\code{FALSE})}
}
\value{
An error if checks are not passed. Runs silently if checks are passed
}
\description{
Validates that a dataset fulfills requirements for MBNMA
}
\details{
Checks done within the validation:
\itemize{
\item Checks data.ab has required column names
\item Checks there are no NAs
\item Checks that all SEs are >0 (if variables are included in dataset)
\item Checks that all doses are >=0
\item Checks that all r and n are positive (if variables are included in dataset)
\item Checks that all y, se, r, n and E are numeric
\item Checks that class codes are consistent within each agent
\item Checks that agent/class names do not contain restricted characters
\item Checks that studies have at least two arms (if \code{single.arm = FALSE})
\item Checks that each study includes at least two treatments
\item Checks that agent names do not include underscores
\item Checks that standsd values are consistent within a study
}
}
