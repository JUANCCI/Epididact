% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/epididact.R
\name{plot2logis}
\alias{plot2logis}
\title{plot2logis}
\usage{
plot2logis(y0_1, y0_2, r_1, r_2, maxt)
}
\arguments{
\item{y0_1}{Epidemic 1: initial inoculum}

\item{y0_2}{Epidemic 2: initial inoculum}

\item{r_1}{Epidemic 1: rate of progress}

\item{r_2}{Epidemic 2: rate of progress}

\item{maxt}{epidemic duration (days)}
}
\value{
two logistic curves
}
\description{
Comparison of two epidemics (logistic models) in the same plot
}
\examples{
plot2logis(0.001, 0.001, 0.4,0.2, 60)
}

