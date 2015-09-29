#'  A function to assess and rank possible locations for their contribution to the metapopulation
#' @inheritParams site.analysis
#'  @param potential.sites A three column matrix providing the name of the site (column 1) and the xy coordinates of current sites (columns 2 & 3)

#'
#'  @usage best.locale(sites,
#' potential.sites,
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
#' @return This function returns a list with two components: (1) \code{$summary.df} is a data frame with the results averaged over all Monte Carlo iterations, and (2) \code{$results.list} is a list of length equal to the number of Monte Carlo iterations and contains the raw results from the analysis.
#'
#' Elements of the summary data frame include (1) \code{cont_lamb}: The contribution of each potential location to the metapopulation; (2) \code{rank_cont}: The rank order importance of each location to the metapopulation; (3) \code{lambda}: The the leading eigen value of the dispersal matrix, which is interpreted as the metapopulation capacity; (4) \code{rank_lamb}: The rank order importance of each location to metapopulation capacity; and (5) \code{avg_rank}: The rank importance of each potential location, determined by averaging \code{rank_cont} and \code{rank_lamb}.
#'
#' @export
#'
#' @details
#' If \code{met.size} is specified, the probability of surviving to adulthood is determined using the equation:\cr
#' logit(p.survive) =  -1.366 + 0.87 * size\cr
#' This equation comes from Altwegg & Reyer (2003). Mean and standard deviation values for \code{met.size} must be reported in standard units such that the mean and standard deviation of observations equal zero and one, respectively (i.e. scale and center observations). In the absence of metamorph, survival probability and variation can be specified using \code{prop.survive}.
#'
#' This model assumes uncertainty or variability in (1) population size; (2) size of metamorphs OR proportion surviving; (3) proportion of population that is philopatric; (4) mean dispersal distance. Uncertainty in these parameters is incorporated through repeated draws from normal distibutions with a mean and standard deviation as specified. Because some values are unrealistic (e.g., survival > 1), a truncated normal distribution is used, which requires the specification of lower and upper values. If there are no limits on the lower or upper values, then \code{-Inf} or \code{Inf} should be specified. Lower and upper values must be provided as a two-element vector (e.g., c(0,1)) for \code{lower.upper_philo}, \code{lower.upper_survive}, and \code{lower.upper_dispersal}
#'
#' Probability of dispersal between two populations is determined using an incidence function wherein the probability of connectivity is a negative exponential relationship with 1/mean dispersal controlling the rate of decay.
#'
#'
#' @seealso
#' \code{\link[ssmc]{ssmc_summary}}
#'
#' @examples
#'    # Assess potential new locations
#'    best_results <- best.locale(sites = site.dat[,1:3],
#'    potential.sites = potential.dat,
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
#'
#'    @references Altwegg, R., and H.-U. Reyer. 2003. Patterns of natural selection on size at metamorphosis in water frogs. Evolution 57:872-882.

best.locale <- function(sites,
                        potential.sites,
                        pop.abun,
                        met.size=NULL,
                        prop.philo=.85,
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

  if(!is.null(seed)){
    set.seed(seed)
  }

  n <- nrow(potential.sites)
  results <- vector("list",n)
  results.df <- vector("list",iterations)
  site.coords <- sites[,2:3]
  pot.sites <- potential.sites
  potential.sites <- potential.sites[,2:3]

  colnames(site.coords) <- colnames(potential.sites) <- c("x","y")
  colnames(pop.abun) <- c('mean.ab','sd')
  colnames(met.size) <- c('svl.mean','svl.sd')


  progress_bar <- plyr::create_progress_bar("text")
  progress_bar$init(n*iterations)

  for (j in 1:iterations) {
    sample.seed <- sample(1:10e6, 1)

    p.philo <- rtnorm(1,prop.philo,sd.philo,lower = lower.upper_philo[[1]], upper = lower.upper_philo[[2]])
    mean.disperse <- rtnorm(1,dispersal,sd.dispersal,lower = lower.upper_dispersal[[1]], lower.upper_dispersal[[2]])

    # Make potential ponds of average population size and metamorph size
    pop.size <- pop.abun %>% rowwise() %>%
      mutate(pop.size = rnorm(1,mean.ab,sd))  %>%
      mutate(pop.size=ifelse(pop.size < 0,0,pop.size)) %>%
      select(pop.size) %>%
      data.frame(.)

    # Set new pond to mean of all occupied ponds
    pop.size <- round(rbind(pop.size, mean(pop.size[pop.size > 0])))

    # Survival--> either based on provided values or result of function
    if(!is.null(met.size)){
      size.dat <- met.size %>%
        rowwise() %>%
        mutate(met.size = rnorm(1,svl.mean,svl.sd))  %>%
        select(met.size) %>%
        data.frame(.)

      p.survive <- plogis(-1.366 + (0.87*size.dat[,1])) %>%
        data.frame(.)

      # Make size-based survival at new pond the average of occupied ponds
      p.survive <- (rbind(p.survive, mean(p.survive[p.survive != plogis(-1.366)])))
      colnames(p.survive) <- c('p.survive')

    } else {

      p.survive <- rtnorm(1,prop.survive,sd.survive,lower = lower.upper_survive[[1]], upper = lower.upper_survive[[2]])

    }

    # Determine philopatric and emigrant populations after accounting for survival
    if(!is.null(met.size)){
      philo.emig <- pop.size %>%
        mutate(survive = rbinom(n(),pop.size$pop.size,p.survive$p.survive),
               philo=rbinom(n(),survive,p.philo),
               emigrant = survive-philo)
    } else {
      philo.emig <- pop.size %>%
        mutate(survive = rbinom(n(),pop.size$pop.size,p.survive),
               philo=rbinom(n(),survive,p.philo),
               emigrant = survive-philo)
    }

    for(i in 1:n){
      # Combine random potential ponds with existing ponds
      site.dat <- rbind(site.coords,potential.sites[i,])

      # Calculate distance matrix with potential pond added, get connectivity
      dist.mat <- as.matrix(dist(site.dat))
      connect.mat <- exp((-1/mean.disperse)*dist.mat)

      # diag(e.mat) <- philo.emig$philo # Set diagonal to philopatric

      connect.mat2 <- t(apply(connect.mat,1,function(x){x/sum(x)}))

      e.mat <- matrix(NaN,nrow(connect.mat2),nrow(connect.mat2))

      set.seed(sample.seed)
      for (r in 1:nrow(connect.mat2)){
        e.mat[r,] <- t(rmultinom(1,philo.emig[r,4],connect.mat2[r,]))
      }

      #~#~#~#~#~#~#~#~#~#
      ## MODULARITY ##
      # g <- graph.adjacency(e.mat,mode = "directed",weighted = T,diag = F)
      # mod <- modularity(cluster_edge_betweenness(g,modularity = T))

      ## METAPOPULATION ##
      # Calculate metapopulaton capacity

      # Get eigen values of emigrant matrix
      tmp <- eigen(e.mat)
      lambda.M <- Re(tmp$value[1])
      lambda.vec <- Re(tmp$vector[, 1])^2 # Contribution of each pond to metapopulation

      metrics <- data.frame(row.names = NULL,
                            iter = j,
                            id = i,
                            x = potential.sites[i,1],
                            y = potential.sites[i,2],
                            met.cap = lambda.M, # Metapopulation capacity
                            met.cont = lambda.vec[nrow(e.mat)]#, # Metapop contribution
                            # mod = mod

      )

      #~#~#~#~#~#~#~#
      # Omit added pond and calculate percent change
      # connect.mat3 <- t(apply(connect.mat[-nrow(connect.mat2),-nrow(connect.mat2)],1,function(x){x/sum(x)}))
      # e.mat <- matrix(NaN,nrow(connect.mat3),nrow(connect.mat3))
      #
      # set.seed(sample.seed)
      # for (r in 1:nrow(connect.mat3)){
      #   e.mat[r,] <- t(rmultinom(1,philo.emig[r,4],connect.mat3[r,]))
      # }
      #
      # tmp <- eigen(e.mat,symmetric = F)
      # lambda.M2 <- Re(tmp$value[1])
      #
      # # Metapopulation capacity change with new population added
      # if(lambda.M==lambda.M2){
      #   metrics$lambda_change <- 0
      # } else if(lambda.M==0 & lambda.M2>0){
      #   metrics$lambda_change <- log(100)*-1
      # } else if(lambda.M==0 & lambda.M2<0){
      #   metrics$lambda_change <- log(100)
      # } else {
      #   metrics$lambda_change <- log((lambda.M)/lambda.M2)
      # }
      #~#~#~#~#~#~#~#

      results[[i]] <- metrics

      progress_bar$step()

    } # Close potential site loop

    results.df[[j]] <- ldply(results,"identity")

  } # End iterations loop

  lambda <- llply(results.df,function(x) x[,5])
  lambda_contrib <- llply(results.df,function(x) x[,6])
  # mod <- llply(results.df,function(x) x[,7])

  lambda_avg <- ldply(lambda,rbind) %>% apply(.,2,mean,na.rm=T)
  lambda_cont <- ldply(lambda_contrib,rbind) %>% apply(.,2,mean,na.rm=T)
  # mod_avg <- ldply(mod,rbind) %>% apply(.,2,mean,na.rm=T)
  # lambda_cont.sd <- ldply(lambda_contrib,rbind) %>% apply(.,2,sd,na.rm=T)
  rank_lambda_avg <-  rank(-lambda_avg)
  rank_lambda_cont <- rank(-lambda_cont)
  # rank_mod <- rank(mod_avg)
  rank.df <- data.frame(rank_lambda_cont, rank_lambda_avg)
  avg_rank <- rank(apply(rank.df,1,mean), ties.method = 'min')

  summary.df <- data.frame(site = pot.sites[,1],
                           x =  potential.sites[,1],
                           y =  potential.sites[,2],
                           # mod_avg = mod,
                           # rank_mod = rank_mod,
                           cont_lam=lambda_cont,
                           rank_cont=rank_lambda_cont,
                           lambda=lambda_avg,
                           rank_lamb=rank_lambda_avg,
                           avg_rank
  )

  out <- list(summary.df = summary.df,
              results.list = results.df)
  return(out)
} # End function
