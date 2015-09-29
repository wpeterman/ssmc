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
site_results <- site.analysis(sites = site.dat[,1:3],
                              pop.abun = site.dat[,4:5],
                              met.size = site.dat[,6:7],
                              prop.philo = 0.95,   
                              sd.philo = 0.05,
                              lower.upper_philo=c(lower=0, upper=1),
                              lower.upper_survive = c(lower=0, upper=1),
                              dispersal = 25,
                              sd.dispersal = 10,
                              lower.upper_dispersal = c(lower=10, upper=Inf),
                              iterations = 5,
                              seed = 123)

## ------------------------------------------------------------------------
site_summary <- ssmc_summary(ssmc_results = site_results,
                             sort = "avg_rank")
str(site_summary)

## ----results='hide',message=FALSE, warning=FALSE-------------------------
best_results <- best.locale(sites = site.dat[,1:3],
                            potential.sites = potential.dat,
                            pop.abun = site.dat[,4:5],
                            met.size = site.dat[,6:7],
                            prop.philo = 0.95,
                            sd.philo = 0.05,
                            lower.upper_philo=c(lower=0, upper=1),
                            lower.upper_survive = c(lower=0, upper=1),
                            dispersal = 25,
                            sd.dispersal = 10,
                            lower.upper_dispersal = c(lower=10, upper=Inf),
                            iterations = 5,
                            seed = 123)

## ------------------------------------------------------------------------
best_summary <- ssmc_summary(ssmc_results = best_results,
                             sort = "avg_rank")
str(best_summary)

