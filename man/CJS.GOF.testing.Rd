% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CJS_goodness.of.fit.R
\name{CJS.GOF.testing}
\alias{CJS.GOF.testing}
\title{Goodness of fit testing of CMR data for a CJS model}
\usage{
CJS.GOF.testing(data.processed, verbose = FALSE)
}
\arguments{
\item{data.processed}{(required) a list with several variables
(e.g., number of capture occasions, time intervals) returned by
\code{\link[RMark]{process.data}()} based on live-only capture histories and
model = "CJS".}

\item{verbose}{(optional, default = FALSE) if TRUE, prints test details}
}
\value{
Returns a dataframe with goodness of fit testing results from overall
CJS, and tests 2CT, 3SR, 3SM, 2CL
}
\description{
Run goodness of fit tests from R2ucare on CJS processed data
}
\details{
This function produces a dataframe of summary stats from each R2ucare
test with the option to print details of each test to the console.
Note: future versions will allow the user to specify group covariates, but
currently this feature is unavailable.
}
\author{
Sarah Gutowsky
}
