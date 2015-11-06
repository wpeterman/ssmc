# About file
#' @docType package
#' @name ssmc-package
#' @title  About this R package
#' @description  This package contains two functions: (1) Assess the average source-sink status of existing populations as well as their rank importance to the metapopulation; (2) Determine the optimal location for habitat creation that will provide the greatest benefit to the metapopulation.
#' @details  \tabular{ll}{
#'    Package: \tab ssmc\cr
#'    Type: \tab Package\cr
#'    License: \tab >=GPL-2\cr
#'  }
#'  This package provides functions to conduct Monte Carlo simulations of demographic connectivity models. These functions can assess the contributions and importance of existing populations and can also assess the potential contributions of populations following habitat creation or restoration.
#'
#' @import dplyr
#' @importFrom msm rtnorm
#' @importFrom plyr create_progress_bar ldply llply progress_text
#' @references Please cite:
#' Manuscript is being prepared for publication
#'
#'
#' @author Bill Peterman \email{Peterman.73@@osu.edu}
#'
NULL

#' Example population data
#'
#' A data frame containing example to use with the \code{\link[ssmc]{site.analysis}} and \code{\link[ssmc]{best.locale}} functions
#'
#' \itemize{
#'    \item site - Numeric site ID
#'    \item x -  x-coordinate
#'    \item y -  y-coordinate
#'    \item A - area of each site
#'    \item ab.mean Mean population abundance of each site
#'    \item ad.sd Standard deviation of mean population abundance at each site
#'    \item svl.mean Mean size (snout-vent length = SVL) for each population. Provided in standard units such that the mean of all observation equals 0 with a standard deviation of 1 (i.e. scaled and centered)
#'    \item svl.sd Standard deviation of standardized mean size variable
#'    }
#'
#' @docType data
#' @name site.dat
#' @format A 25 x 8 data frame
#' @usage data(site.dat)
#' @keywords datasets
#'
NULL

#' Example potential site data
#'
#' A data frame containg the necessary components to use the \code{\link[ssmc]{best.locale}} function
#'
#'
#' \itemize{
#'    \item site - Numeric site ID
#'    \item x -  x-coordinate
#'    \item y -  y-coordinate
#'    }
#'
#' @docType data
#' @name potential.dat
#' @format A 25 x 3 data frame
#' @usage data(potential.dat)
#' @description  Sample file to be used with \code{\link[ssmc]{best.locale}}
#' @keywords datasets
NULL
