% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssmc_data_description.R
\docType{data}
\name{site.dat}
\alias{site.dat}
\title{Example population data}
\format{
A 25 x 8 data frame
}
\usage{
data(site.dat)
}
\description{
A data frame containing example data to use with the \code{\link[ssmc]{site.analysis}} and \code{\link[ssmc]{best.locale}} functions
}
\details{
\itemize{
   \item site - Numeric site ID
   \item x -  x-coordinate
   \item y -  y-coordinate
   \item A - area of each site
   \item ab.mean Mean population abundance of each site
   \item ab.sd Standard deviation of mean population abundance at each site
   \item svl.mean Mean size (snout-vent length = SVL) for each population. Provided in standard units such that the mean of all observation equals 0 with a standard deviation of 1 (i.e. scaled and centered)
   \item svl.sd Standard deviation of standardized mean size variable
   }
}
\keyword{datasets}
