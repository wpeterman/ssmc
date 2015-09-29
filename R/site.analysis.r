#'  A function to assess average source-sink status of populations as well as rank importance of populations in contributing to the metapopulation
#'
#'  @param sites A three-column matrix or data frame providing the name of the site (column 1) and the xy coordinates of current sites (columns 2 & 3)
#'  @param pop.abun A two-column data frame with mean population size in column 1 and standard deviation in column 2
#'  @param met.size (Default = NULL) If specified, must be a two-column data frame with the mean (column 1) and standard deviation (column 2) of the size of metamoprhs or late-stage larvae. See Details for more information concerning this parameter.
#'  @param pop.philo (Default = 0.85) Mean proportion of population that are philopatric to their natal population
#'  @param sd.philo (Default = 0.05) Standard deviation of proportion of population that are philopatric to their natal population
#'  @param lower.upper_philo (Default = c(lower=0, upper=1)) Threshold lower and upper values for the proportion of philopatric individuals. Must be provided as a two-element vector with lower value first. See Details for use
#'  @param prop.survive (Default = 0.2) Mean survival of to adulthood. Ignored if \code{met.size} is specified. See Details for more information.
#'  @param sd.survive (Default = 0.075) Standard deviation of survival of to adulthood.
#'  @param lower.upper_survive (Default = c(lower=0, upper=1))Threshold lower and upper values for the proportion of individuals surviving to adulthood. Must be provided as a two-element vector with lower value first. See Details
#'  @param dispersal Mean dispersal distance. Should be specified in units meaningful to the coordinate system describing the spatial location of populations. See Details for how dispersal is estimated.
#'  @param sd.dispersal Standard deviation of dispersal distance.
#'  @param lower.upper_dispersal (Default = c(lower=1, upper=Inf))Threshold lower and upper values for average dispersal distance. Must be provided as a two-element vector with minimum value first. See Details
#'  @param iterations Number of Monte Carlo iterations to run
#'  @param seed Optional to set the seed for repeatability  among model runs (Default = NULL)
#'
#'  @usage site.analysis(sites,
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
#' (7) \code{cont_lamb}: The contribution of each location to the metapopulation; \cr
#' (8) \code{rank_cont}: The rank order importance of each location to the metapopulation; \cr
#' (9) \code{delta_lam}: The change in the value of the the leading eigen value of the landscape matrix when each population is removed, which is interpreted as the metapopulation capacity. This is calculated as the log ratio of values [delta_lam = log(lambda.new / lambda.old)]; \cr
#' (10) \code{rank_lamb}: The rank order importance of each location to metapopulation capacity; \cr
#' (11) \code{avg_rank}: The rank importance of each potential location, determined by averaging \code{rank_cont} and \code{rank_lamb}.
#'
#' @export
#'
#' @details
#' If \code{met.size} is specified, the probability of surviving to adulthood is determined using the equation:\cr
#' logit(p.survive) =  -1.366 + 0.87 * size\cr
#' This equation comes from Altwegg & Reyer (2003). Mean and standard deviation values for \code{met.size} must be reported in standard units such that the mean and standard deviation of observations equal zero and one, respectively (i.e. scale and center observations). In the absence of metamorph, survival probability and variation can be specified using \code{prop.survive}.
#'
#' This model assumes uncertainty or variability in (1) population size; \cr \cr (2) size of metamorphs OR proportion surviving; \cr \cr (3) proportion of population that is philopatric; \cr \cr (4) mean dispersal distance. Uncertainty in these parameters is incorporated through repeated draws from normal distibutions with a mean and standard deviation as specified. Because some values are unrealistic (e.g., survival > 1), a truncated normal distribution is used, which requires the specification of lower and upper values. If there are no limits on the lower or upper values, then \code{-Inf} or \code{Inf} should be specified. Lower and upper values must be provided as a two-element vector (e.g., c(0,1)) for \code{lower.upper_philo}, \code{lower.upper_survive}, and \code{lower.upper_dispersal}
#'
#' Probability of dispersal between two populations is determined using an incidence function wherein the probability of connectivity is a negative exponential relationship with 1/mean dispersal controlling the rate of decay.
#'
#' @seealso
#' \code{\link[ssmc]{ssmc_summary}}
#'
#'    @references Altwegg, R., and H.-U. Reyer. 2003. Patterns of natural selection on size at metamorphosis in water frogs. Evolution 57:872-882.
#'
#' @examples
#'    # Assess existing populations
#'    site_results <- site.analysis(sites = site.dat[,1:3],
#'    pop.abun = site.dat[,4:5],
#'    met.size = site.dat[,6:7],
#'    prop.philo = 0.95,
#'    sd.philo = 0.05,
#'    lower.upper_philo=c(lower=0, upper=1),
#'    # prop.survive=0.2,      ## Not needed b/c met.size is specified
#'    # sd.survive = 0.075,    ## Not needed b/c met.size is specified
#'    lower.upper_survive = c(lower=0, upper=1),
#'    dispersal = 25,
#'    sd.dispersal = 10,
#'    lower.upper_dispersal = c(lower=10, upper=Inf),
#'    iterations = 10,
#'    seed = 123)


site.analysis <- function(sites, # Pond names in column 1, coords in col 2&3
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
                          iterations,
                          seed=NULL){
  if(!is.null(seed)) {
    set.seed(seed)
  }

  colnames(sites) <-  c('site','x','y')
  colnames(pop.abun) <- c('mean.ab','sd')
  colnames(met.size) <- c('svl.mean','svl.sd')

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
    sample.seed <- sample(1:10e6, 1)

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
    results <- data.frame(immig = colSums(e.mat),
                          philo = philo.emig[,3],
                          emig = philo.emig[,4])

    results <- results %>%
      mutate(pct_immig = ifelse(philo>0 & immig!=0, immig/(immig+philo),
                                ifelse(philo>=0 & immig==0,0,1)))


    results[is.na(results)] <- 0
    results$sink <- ifelse(results$immig>=results$philo,1,0)


    #~#~#~#~#~#~#~#~#~#
    ## MODULARITY ##
    # g <- graph.adjacency(e.mat,mode = "directed",weighted = T,diag = F)
    # mod <- modularity(cluster_edge_betweenness(g,modularity = T))

    ## METAPOPULATION ##

    # Get eigen values of emigrant matrix


    # diag(e.mat) <- philo.emig$philo # Set diagonal to philopatric

    tmp <- eigen(e.mat,symmetric = F)
    lambda.M <- Re(tmp$value[1])
    results$lambda <- Re(tmp$vector[, 1]^2)


    ## Assess importance of each pond by omitting it
    for(j in 1:dim(connect.mat)[1]){

      connect.mat3 <- t(apply(connect.mat[-j,-j],1,function(x){x/sum(x)}))
      e.mat2 <- matrix(NaN,nrow(connect.mat3),nrow(connect.mat3))

      set.seed(sample.seed)
      for (r in 1:nrow(connect.mat3)){
        p.e <- philo.emig[-j,]
        e.mat2[r,] <- t(rmultinom(1,p.e[r,4],connect.mat3[r,]))
      }

      # diag(e.mat2) <- p.e$philo # Set diagonal to philopatric
      tmp <- eigen(e.mat2,symmetric = F)
      lambda.M2 <- Re(tmp$value[1])

      # Metapopulation capacity change
      oldw <- getOption("warn")
      options(warn = -1)
      if(lambda.M2!=0 & lambda.M!=0) {
        results$lambda_change[j] <- log((lambda.M2)/lambda.M)
      } else if(lambda.M2 > lambda.M){
        results$lambda_change[j] <- log(100)
      } else {
        results$lambda_change[j] <- log(100)*-1
      }
      options(warn = oldw)

      # if(lambda.M==lambda.M2){
      #   results$lambda_change[j] <- 0
      # } else if(lambda.M==0 & lambda.M2>0){
      #   results$lambda_change[j] <- log(100)
      # } else if(lambda.M>0 & lambda.M2==0){
      #   results$lambda_change[j] <- log(100)*-1
      # } else {
      #   results$lambda_change[j] <- log((lambda.M2)/lambda.M)
      # }

      ## MODULARITY ##
      # g2 <- graph.adjacency(e.mat2,mode = "directed",weighted = T,diag = F)
      # mod2 <- modularity(cluster_edge_betweenness(g2,modularity = T))

      # Change in modularity
      # Metapopulation capacity change
      # if(mod==mod2){
      #   results$mod_change[j] <- 0
      # } else if(mod==0 & mod2>0){
      #   results$mod_change[j] <- log(100)
      # } else if(mod==0 & mod2<0){
      #   results$mod_change[j] <- log(100)*-1
      # } else {
      #   results$mod_change[j] <- log((mod2)/mod)
      # }

    } # End leave-one-out loop


    results <- cbind(data.frame(iter=i),results)

    results.list[[i]] <- results

    progress_bar$step()
  } # End iterations loop

  immig <- llply(results.list,function(x) x[,2])
  philo <- llply(results.list,function(x) x[,3])
  emig <- llply(results.list,function(x) x[,4])
  pct_immig <- llply(results.list,function(x) x[,5])
  sink <- llply(results.list,function(x) x[,6])
  lambda_contrib <- llply(results.list,function(x) x[,7])
  lambda_change <- llply(results.list, function(x) x[,8])
  # mod_change <- llply(results.list, function(x) x[,9])

  immig <- ldply(immig,rbind) %>%
    apply(.,2,mean,na.rm=T)
  philo <- ldply(philo,rbind) %>%
    apply(.,2,mean,na.rm=T)
  emig <- ldply(emig,rbind) %>%
    apply(.,2,mean,na.rm=T)
  pct_immig <- ldply(pct_immig,rbind) %>%
    apply(.,2,mean,na.rm=T)
  pct_sink <- ldply(sink,rbind) %>%
    apply(.,2,mean,na.rm=T)
  lambda_cont <- ldply(lambda_contrib,rbind) %>%
    apply(.,2,mean,na.rm=T)
  # lambda_cont.sd <- ldply(lambda_contrib,rbind) %>%
  # apply(.,2,sd,na.rm=T)
  pc_lambda <- ldply(lambda_change,rbind) %>%
    apply(.,2,mean,na.rm=T)
  # delta_mod <- ldply(mod_change,rbind) %>%
  # apply(.,2,mean,na.rm=T)
  # sd_lambda <- ldply(lambda_change,rbind) %>%
  # apply(.,2,sd,na.rm=T)
  rank_lambda_cont <- rank(-lambda_cont)
  rank_lambda_change <- rank(pc_lambda)
  rank.df <- data.frame(rank_lambda_change, rank_lambda_cont)
  avg_rank <- rank(apply(rank.df,1,mean), ties.method = 'min')

  out <- data.frame(site = sites[,1],
                    x = site.coords[,1],
                    y = site.coords[,2],
                    immig = immig,
                    philo = philo,
                    emig = emig,
                    pct_immig=pct_immig,
                    pct_sink=pct_sink,
                    pct_src=1-pct_sink,
                    # delta_mod=delta_mod,
                    cont_lam=lambda_cont,
                    # cont_sd=lambda_cont.sd,
                    rank_cont=rank_lambda_cont,
                    delta_lam=pc_lambda,
                    # sd_lambda = sd_lambda,
                    rank_lamb=rank_lambda_change,
                    avg_rank=avg_rank
  )

  ss.results <- list(summary.df=out,
                     results.list=results.list)
  return(ss.results)

} # End function
