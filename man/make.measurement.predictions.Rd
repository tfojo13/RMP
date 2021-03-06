% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{make.measurement.predictions}
\alias{make.measurement.predictions}
\title{Make predictions of Repeated Measurements}
\usage{
make.measurement.predictions(dm, observations, observed.times, predict.times,
  covariates = NULL, ci.coverage = 0.95, iq.coverage = 0.5, sum = F)
}
\arguments{
\item{observations}{an array indexed [i,time,test] - where i denotes an individual}

\item{observed.times}{a numeric vector such that observed.times[t] is the time at which observations[,t,] were taken}

\item{covariates}{A row for each individual with the covariates used to fit the model}
}
\value{
A 4d array indexed[individual, time, test, statistic], where statistic is one of: mean, ci.lower, ci.upper, iq.lower, iq.upper
}
