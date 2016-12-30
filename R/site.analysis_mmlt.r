#'  A function to assess average source-sink status of populations as well as rank importance of populations in contributing to the metapopulation mean lifetime
#'
#' @param sites A four-column matrix or data frame providing the name of the site (column 1) and the xy coordinates of current sites (columns 2 & 3), and area of the habitat patch (column 4)
#' @param pop.abun A two-column data frame with mean population size in column 1 and standard deviation in column 2
#' @param met.size (Default = NULL) If specified, must be a two-column data frame with the mean (column 1) and standard deviation (column 2) of the size of metamoprhs or late-stage larvae. See Details for more information concerning this parameter.
#' @param prop.philo (Default = 0.85) Mean proportion of population that are philopatric to their natal population
#' @param sd.philo (Default = 0.05) Standard deviation of proportion of population that are philopatric to their natal population
#' @param lower.upper_philo (Default = c(lower=0, upper=1)) Threshold lower and upper values for the proportion of philopatric individuals. Must be provided as a two-element vector with lower value first. See Details for use
#' @param prop.survive (Default = 0.2) Mean survival of to adulthood. Ignored if \code{met.size} is specified. See Details for more information.
#' @param sd.survive (Default = 0.075) Standard deviation of survival of to adulthood.
#' @param lower.upper_survive (Default = c(lower=0, upper=1))Threshold lower and upper values for the proportion of individuals surviving to adulthood. Must be provided as a two-element vector with lower value first. See Details
#' @param dispersal Mean dispersal distance. Should be specified in units meaningful to the coordinate system describing the spatial location of populations. See Details for how dispersal is estimated.
#' @param sd.dispersal Standard deviation of dispersal distance.
#' @param lower.upper_dispersal (Default = c(lower=1, upper=Inf))Threshold lower and upper values for average dispersal distance. Must be provided as a two-element vector with minimum value first. See Details
#' @param eps Coefficient ([0,1]) relating to minimum patch size (Default = 1). See Details
#' @param mu Number of immigrants needed for successful colonization (Default = 2). See Details
#' @param eta Scaling parameter (Default = 0.5). See Details
#' @param iterations Number of Monte Carlo iterations to run
#' @param seed Optional to set the seed for repeatability  among model runs (Default = NULL)
#'
#' @usage site.analysis(sites,
#' pop.abun,
#' met.size=NULL,
#' prop.philo=.85,
#' sd.philo = 0.05,
#' lower.upper_philo=c(lower=0, upper=1),
#' prop.survive=0.2,
#' sd.survive = 0.075,
#' lower.upper_survive = c(lower=0, upper=1),
#' dispersal,
#' sd.dispersal,
#' lower.upper_dispersal = c(lower=1, upper=Inf),
#' eps = 1,
#' mu = 2,
#' eta = 0.5,
#' iterations,
#' seed=NULL)
#'
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @return
#' This function returns a list with two components: (1) \code{$summary.df} is a data frame with the results averaged over all Monte Carlo iterations, and (2) \code{$results.list} is a list of length equal to the number of Monte Carlo iterations and contains the raw results from the analysis.
#'
#' Elements of the summary data frame include: \cr
#' (1) \code{immig}: The average number of immigrants entering each population; \cr
#' (2) \code{philo}: The average number of philopatric individuals; \cr
#' (3) \code{emig}: The average number of emigrants dispersing from a population; \cr
#' (4) \code{pct_immig}: Average percentage of a population that is comprised of immigrants; \cr
#' (5) \code{pct_sink}: The frequency that a population acted as a sink (i.e. immigrants > philopatric) across all simulations; \cr
#' (6) \code{pct_src}: The frequency that a population acted as a source (i.e. immigrant <= philopatric) across all simulations; \cr
#' (7) \code{delt_mmlt}: Average change in the log of the metapopulation mean lifetime; \cr
#' (8) \code{rank}: The rank order importance of each location based on the change in mmlt
#'
#' @export
#'
#' @details
#' If \code{met.size} is specified, the probability of surviving to adulthood is determined using the equation:\cr
#'
#' logit(p.survive) =  -1.366 + 0.87 * size\cr
#'
#' This equation comes from Altwegg & Reyer (2003). Mean and standard deviation values for \code{met.size} must be reported in standard units such that the mean and standard deviation of observations equal zero and one, respectively (i.e. scale and center observations). In the absence of metamorph, survival probability and variation can be specified using \code{prop.survive}.
#'
#' This model assumes uncertainty or variability in:\cr
#'  (1) population size; \cr
#'  (2) size of metamorphs OR proportion surviving; \cr
#'  (3) proportion of population that is philopatric; \cr
#'  (4) mean dispersal distance.
#'
#'  Uncertainty in these parameters is incorporated through repeated draws from normal distibutions with a mean and standard deviation as specified. Because some values are unrealistic (e.g., survival > 1), a truncated normal distribution is used, which requires the specification of lower and upper values. If there are no limits on the lower or upper values, then \code{-Inf} or \code{Inf} should be specified. Lower and upper values must be provided as a two-element vector (e.g., c(0,1)) for \code{lower.upper_philo}, \code{lower.upper_survive}, and \code{lower.upper_dispersal}
#'
#' Probability of dispersal between two populations is determined using an incidence function wherein the probability of connectivity is a negative exponential relationship with 1/mean dispersal controlling the rate of decay.
#'
#' @seealso
#' \code{\link[ssmc]{ssmc_summary}}
#'
#' @references Altwegg, R., and H.-U. Reyer. 2003. Patterns of natural selection on size at metamorphosis in water frogs. Evolution 57:872-882.
#'
#' @examples
#'    # Assess existing populations
#'    site_results <- site.analysis(sites = site.dat[,1:4],
#'    pop.abun = site.dat[,5:6],
#'    met.size = site.dat[,7:8],
#'    prop.philo = 0.95,
#'    sd.philo = 0.05,
#'    lower.upper_philo=c(lower=0, upper=1),
#'    # prop.survive=0.2,      ## Not needed b/c met.size is specified
#'    # sd.survive = 0.075,    ## Not needed b/c met.size is specified
#'    lower.upper_survive = c(lower=0, upper=1),
#'    dispersal = 25,
#'    sd.dispersal = 10,
#'    lower.upper_dispersal = c(lower=10, upper=Inf),
#'    eps = 1,
#'    mu = 2,
#'    eta = 0.5,
#'    iterations = 10,
#'    seed = 123)


site.analysis <- function(sites, # Pond names in column 1, coords in col 2&3, area in column 4
                          pop.abun, # Mean and sD for each pond
                          met.size=NULL, # Standardized mean and sd for each population
                          prop.philo=0.85,
                          sd.philo = 0.05,
                          lower.upper_philo=c(lower=0, upper=1),
                          prop.survive=0.2,
                          sd.survive = 0.075,
                          lower.upper_survive = c(lower=0, upper=1),
                          dispersal,
                          sd.dispersal,
                          lower.upper_dispersal = c(lower=1, upper=Inf),
                          eps = 1,
                          mu = 2,
                          eta = 0.5,
                          iterations,
                          seed=NULL){
  if(!is.null(seed)) {
    set.seed(seed)
  }

  colnames(sites) <-  c('site','x','y','area')
  colnames(pop.abun) <- c('mean.ab','sd')

  if(!is.null(met.size)){
    colnames(met.size) <- c('svl.mean','svl.sd')
  }

  # Process lower.upper values
  if(length(lower.upper_philo)!=2 |
     length(lower.upper_dispersal)!=2) {
    stop("Must specify minimum and maximum values for philopatry and dispersal. See Details for help")
  }

  if(is.null(met.size) & length(lower.upper_survive)!=2) {
    stop("Must specify minimum and maximum values for survival. See Details for help")
  }

  n <- nrow(sites)
  results.list <- vector("list",iterations)
  #   results.df <- vector("list",iterations)
  site.coords <- sites[,2:3]
  colnames(site.coords) <- c("x","y")

  progress_bar <- plyr::create_progress_bar("text")
  progress_bar$init(iterations)

  for(i in 1:iterations){
    sample.seed <- sample(1:10e6, 1) + i

    # Get values for philopatry and mean dispersal from distribution
    p.philo <- rtnorm(1,prop.philo,sd.philo,lower = lower.upper_philo[[1]], upper = lower.upper_philo[[2]])
    mean.disperse <- rtnorm(1,dispersal,sd.dispersal,lower = lower.upper_dispersal[[1]], lower.upper_dispersal[[2]])

    # Survival--> either based on provided values or result of function
    if(!is.null(met.size)){
      size.dat <- met.size %>% rowwise() %>% mutate(met.size = rnorm(1,svl.mean,svl.sd))  %>%  select(met.size) %>% data.frame(.)

      p.survive <- plogis(-1.366 + (0.87*size.dat[,1])) %>% data.frame(.)
      p.survive <- (rbind(p.survive, mean(p.survive[p.survive != plogis(-1.366)])))
      colnames(p.survive) <- c('p.survive')

    } else {
      p.survive <- rtnorm(1,prop.survive,sd.survive,lower = lower.upper_survive[[1]], upper = lower.upper_survive[[2]])
    }

    # Determine population size from distribution
    pop.size <- pop.abun %>%
      rowwise() %>%
      mutate(pop.size = rnorm(1,mean.ab,sd))  %>%
      mutate(pop.size=ifelse(pop.size<0,0,pop.size)) %>%
      select(pop.size) %>%
      round(.) %>%
      data.frame(.)

    # Determine philopatric and emigrant populations after accounting for survival
    if(!is.null(met.size)){
      philo.emig <- pop.size %>%
        mutate(survive = rbinom(n(), pop.size$pop.size, p.survive$p.survive),
               philo=rbinom(n(), survive, p.philo),
               emigrant = survive-philo)
    } else {
      philo.emig <- pop.size %>%
        mutate(survive = rbinom(n(), pop.size$pop.size, p.survive),
               philo=rbinom(n(), survive, p.philo),
               emigrant = survive-philo)
    }

    # Calculate distance matrix, get connectivity
    dist.mat <- as.matrix(dist(site.coords))
    connect.mat <- exp((-1/mean.disperse)*dist.mat)
    diag(connect.mat) <- 0
    connect.mat2 <- t(apply(connect.mat,1,function(x){x/sum(x)}))

    e.mat <- matrix(NaN,nrow(connect.mat2),nrow(connect.mat2))

    set.seed(sample.seed)
    for (r in 1:nrow(connect.mat2)){
      e.mat[r,] <- t(rmultinom(1,philo.emig[r,4],connect.mat2[r,]))
    }


    # Philopatric total
    e.mat_bin <- ifelse(e.mat==0, 0, 1)
    results <- data.frame(immig = colSums(e.mat),
                          philo = philo.emig[,3],
                          emig = philo.emig[,4],
                          in_deg = colSums(e.mat_bin),
                          out_deg = rowSums(e.mat_bin))

    results <- results %>%
      mutate(pct_immig = ifelse(philo>0 & immig!=0, immig/(immig+philo),
                                ifelse(philo>=0 & immig==0,0,1)))


    results[is.na(results)] <- 0
    results$sink <- ifelse(results$immig>=results$philo,1,0)


    #~#~#~#~#~#~#~#~#~#

    ## METAPOPULATION Mean Lifetime ##
    mmlt <- MMLT(N = dim(e.mat)[1],
                 S_i = results$immig,
                 S_o = results$emig,
                 eps = eps,
                 A = sites$area,
                 mu = mu,
                 eta = eta)

    # lam_c <- philo.emig$survive
    # lam_M <- sum(lam_c * lam_c/sum(lam_c))

    ## Assess importance of each pond by omitting it
    for(j in 1:dim(connect.mat)[1]){

      connect.mat3 <- t(apply(connect.mat[-j,-j],1,function(x){x/sum(x)}))
      e.mat2 <- matrix(NaN,nrow(connect.mat3),nrow(connect.mat3))

      set.seed(sample.seed)
      for (r in 1:nrow(connect.mat3)){
        p.e <- philo.emig[-j,]
        S_i = results$immig[-j]
        S_o = results$emig[-j]
        A = sites$area[-j]
        e.mat2[r,] <- t(rmultinom(1,p.e[r,4],connect.mat3[r,]))
      }

      mmlt2 <- MMLT(N = dim(e.mat2)[1],
                    S_i = S_i,
                    S_o = S_o,
                    eps = eps,
                    A = A,
                    mu = mu,
                    eta = eta)

      # lam_c2 <- p.e$survive
      # lam_M2 <- sum(lam_c2 * lam_c2/sum(lam_c2))

      # MMLT Change
      results$delt_mmlt[j] <- (mmlt2 - mmlt) / mmlt
      # results$delt_mmlt[j] <- (log(mmlt2 + 1, 10) - log(mmlt + 1, 10)) / log(mmlt + 1, 10)
      # results$delt_lam[j] <- (lam_M2 - lam_M)/lam_M


      # if(lambda.M==lambda.M2){
      #   results$lambda_change[j] <- 0
      # } else if(lambda.M==0 & lambda.M2>0){
      #   results$lambda_change[j] <- log(100)
      # } else if(lambda.M>0 & lambda.M2==0){
      #   results$lambda_change[j] <- log(100)*-1
      # } else {
      #   results$lambda_change[j] <- log((lambda.M2)/lambda.M)
      # }

    } # End leave-one-out loop


    results <- cbind(data.frame(iter=i),results)

    results.list[[i]] <- results

    progress_bar$step()
  } # End iterations loop

  immig <- llply(results.list,function(x) x[,2])
  philo <- llply(results.list,function(x) x[,3])
  emig <- llply(results.list,function(x) x[,4])
  in_deg <- llply(results.list,function(x) x[,5])
  out_deg <- llply(results.list,function(x) x[,6])
  pct_immig <- llply(results.list,function(x) x[,7])
  sink <- llply(results.list,function(x) x[,8])
  delta_mmlt <- llply(results.list,function(x) x[,9])

  immig <- ldply(immig,rbind) %>%
    apply(.,2,mean,na.rm=T)
  philo <- ldply(philo,rbind) %>%
    apply(.,2,mean,na.rm=T)
  emig <- ldply(emig,rbind) %>%
    apply(.,2,mean,na.rm=T)
  in_deg <- ldply(in_deg,rbind) %>%
    apply(.,2,mean,na.rm=T)
  out_deg <- ldply(out_deg,rbind) %>%
    apply(.,2,mean,na.rm=T)
  pct_immig <- ldply(pct_immig,rbind) %>%
    apply(.,2,mean,na.rm=T)
  pct_sink <- ldply(sink,rbind) %>%
    apply(.,2,mean,na.rm=T)
  delt_mmlt <- ldply(delta_mmlt,rbind) %>%
    apply(.,2,mean,na.rm=T)

  rank_mmlt <- rank(delt_mmlt, ties.method = 'min')

  out <- data.frame(site = sites[,1],
                    x = site.coords[,1],
                    y = site.coords[,2],
                    immig = immig,
                    philo = philo,
                    emig = emig,
                    in_deg = in_deg,
                    out_deg = out_deg,
                    pct_immig=pct_immig,
                    pct_sink=pct_sink,
                    pct_src=1-pct_sink,
                    delt_mmlt=delt_mmlt,
                    rank=rank_mmlt
  )

  ss.results <- list(summary.df=out,
                     results.list=results.list)
  return(ss.results)
  gc()  ## Flush memory

} # End function

MMLT <- function(N, S_i, S_o, eps, A, mu, eta, method = "approx") {
  # Equations from Kininmonth et al. 2010
  # Extinction rate
  Vi <- eps*(A^-eta) # eq 1
  C_out <- S_o/mu # eq 3
  C_in <- S_i/mu # eq 4
  u_out <- C_out/Vi # eq 6
  u_in <- C_in/Vi # eq 7

  # Harmonic mean
  Ui <- (.5*(u_in)^-2 + .5*(u_out)^-2)^-.5 # (Drechsler 2009 eq 8)
  Ui[Ui<sqrt(2)] <- sqrt(2) # Take Ui or sqrt(2), whichever is larger, eq 9, part 1

  q <- prod(Ui^(1/N)) # eq 9 --> Colinization-extinction ratio
  v <- prod(Vi^(1/N)) # eq 10 (Geometric local extinction rates)


  ## Default is to use approximation
  if(method == "approx") {
    MMLT <- N*((1/q) + log(q) - 1) # Eq 12
    MMLT <- exp(MMLT) / v
    MMLT <- log10(MMLT + 1)
    MMLT <- ifelse(is.infinite(MMLT),NA, MMLT)
  } else {
    if(q>2.5 & N>150){
      MMLT <- N*((1/q) + log(q) - 1) # Eq 12
      MMLT <- exp(MMLT) / v
      MMLT <- log10(MMLT + 1)
      MMLT <- ifelse(is.infinite(MMLT),NA, MMLT)

      return(MMLT)

    } else {   # Requiring factorial limits size of network that can be analyzed (eq. 11)
      ## Eq 11
      MMLT <- (1/v)*sum(outer(1:N, 1:N, mmlt_fun, q, N))
      MMLT <- log10(MMLT + 1)
      return(MMLT)
    } ## End inner ifelse
  } ## End outer ifelse
} ## End function

#!#!#!#!#!#!#!
# Compliments of Grant Connette...nearly 100x faster
# Grant <- function(){
#   MMT_fun <- function(i,k){(1/k * (factorial(N-i)/factorial(N-k)) * (1/(N-1)^(k-i)) * q^(k-i))}
#   MMT <- (1/v)*sum(outer(1:N, 1:N, MMT_fun))
#   log(MMT)
# }
#
# q <- 4.69
# v <- 8.4e-15
# N <- 140

mmlt_fun <- function(i,k,q,N){(1/k * (factorial(N-i)/factorial(N-k)) * (1/(N-1)^(k-i)) * q^(k-i))}
  # MMLT2 <- (1/v)*sum(outer(1:N, 1:N, mmlt_fun, q, N))

# Dan <- function() {
# mmlt_fun <- function(x) {
#   mmlt <- (1/x["k"] * (factorial(x["N"]-x["i"])/factorial(x["N"]-x["k"])) * (1/(x["N"]-1)^(x["k"]-x["i"])) * x["q"]^(x["k"]-x["i"]))
# }
#
#   vec <- apply(df, MARGIN = 1, FUN = mmlt_fun)
#   MMLT <- (1/v) * sum(vec)
#   log(MMLT)
#
# }

# Dan's function
# mmlt_fun2 <- function(x) {
#   mmlt <- (1/x["k"] * (gamma(as(x["N"]-x["i"],'mpfr'))/gamma(as(x["N"]-x["k"],'mpfr'))) * (1/(x["N"]-1)^(x["k"]-x["i"])) * x["q"]^(x["k"]-x["i"]))
#   return(mmlt)
# }

