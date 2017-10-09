## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install.package, eval=FALSE-----------------------------------------
#  # Install 'devtools' package, if needed
#  if(!("devtools" %in% list.files(.libPaths()))) {
#      install.packages("devtools",
#                       repo = "http://cran.rstudio.com",
#                       dep = TRUE)
#  }
#  
#  library(devtools) # Loads devtools
#  
#  install_github("wpeterman/ssmc") # Download package

## ----results='hide',message=FALSE, warning=FALSE-------------------------
library(ssmc)
rm(list = ls())

## ------------------------------------------------------------------------
str(site.dat)

## ----results='hide',message=FALSE, warning=FALSE-------------------------
    site_results <- site.analysis(sites = site.dat[,1:4],
                                  pop.abun = site.dat[,5:6],
                                  met.size = site.dat[,7:8],
                                  prop.philo = 0.95,
                                  sd.philo = 0.05,
                                  lower.upper_philo=c(lower=0, upper=1),
                                  lower.upper_survive = c(lower=0, upper=1),
                                  dispersal = 25,
                                  sd.dispersal = 10,
                                  lower.upper_dispersal = c(lower=10, upper=Inf),
                                  eps = 1,
                                  mu = 2,
                                  eta = 0.5,
                                  iterations = 10,
                                  seed = 123)

## ------------------------------------------------------------------------
site_summary <- ssmc_summary(ssmc_results = site_results)
str(site_summary)

## ----results='hide',message=FALSE, warning=FALSE-------------------------
best_results <- best.locale(sites = site.dat[,1:4],
                            potential.sites = potential.dat,
                            pop.abun = site.dat[,5:6],
                            met.size = site.dat[,7:8],
                            prop.philo = 0.95,
                            sd.philo = 0.05,
                            lower.upper_philo=c(lower=0, upper=1),
                            lower.upper_survive = c(lower=0, upper=1),
                            dispersal = 25,
                            sd.dispersal = 10,
                            lower.upper_dispersal = c(lower=10, upper=Inf),
                            eps = 1,
                            mu = 2,
                            eta = 0.5,
                            iterations = 10,
                            seed = 123)

## ------------------------------------------------------------------------
best_summary <- ssmc_summary(ssmc_results = best_results)
str(best_summary)

