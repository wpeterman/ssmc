#'  A function to assess and rank possible locations for their contribution to the metapopulation
#' @inheritParams site.analysis
#' @param potential.sites A three column matrix providing the name of the site (column 1) and the xy coordinates of current sites (columns 2 & 3)
#' @param restore Logical. If TRUE, then \code{potential.sites} does not need to be specified. Instead, existing sites provided in \code{sites} that are currently unoccupied (abundance = 0) will be iteratively assessed at the mean of all occupied sites.

#' @usage best.locale(sites,
#' restore = FALSE,
#' potential.sites = NULL,
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
#' @return This function returns a list with two components: (1) \code{$summary.df} is a data frame with the results averaged over all Monte Carlo iterations, and (2) \code{$results.list} is a list of length equal to the number of Monte Carlo iterations and contains the raw results from the analysis.
#'
#' Elements of the summary data frame include: \cr
#' (1) \code{mmlt}: The log10 of the Metapopulation Mean Lifetime;\cr
#' (2) \code{rank}: The rank order importance based on mmlt;\cr
#' (3) \code{frq_col}: The frequency that 2 or more individuals immigrated into a population. This metric is used to adjust mmlt for likelihood of colonization;\cr
#' (4) \code{mmlt_col}: The product of mmlt and frq_col;\cr
#' (5) \code{rank_adj}: The rank importance of each potential location based on the product of mmlt and the frequency of colonization (frq_col)
#'
#' @export
#'
#' @details
#' If \code{met.size} is specified, the probability of surviving to adulthood is determined using the equation:\cr
#'
#' logit(p.survive) =  -1.366 + 0.87 * size\cr
#'
#' This equation comes from Altwegg & Reyer (2003). Mean and standard deviation values for \code{met.size} must be reported in standard units such that the mean and standard deviation of observations equal zero and one, respectively (i.e. scale and center observations). In the absence of metamorph data, survival probability and variation can be specified using \code{prop.survive}.
#'
#' This model assumes uncertainty or variability in: \cr
#' (1) population size; \cr
#' (2) size of metamorphs OR proportion surviving; \cr
#' (3) proportion of population that is philopatric; \cr
#' (4) mean dispersal distance. \cr
#'
#' Uncertainty in these parameters is incorporated through repeated draws from normal distibutions with a mean and standard deviation as specified. Because some values are unrealistic (e.g., survival > 1), a truncated normal distribution is used, which requires the specification of lower and upper values. If there are no limits on the lower or upper values, then \code{-Inf} or \code{Inf} should be specified. Lower and upper values must be provided as a two-element vector (e.g., c(0,1)) for \code{lower.upper_philo}, \code{lower.upper_survive}, and \code{lower.upper_dispersal}
#'
#' Probability of dispersal between two populations is determined using an incidence function wherein the probability of connectivity is a negative exponential relationship with 1/mean dispersal controlling the rate of decay.
#'
#'
#' @seealso
#' \code{\link[ssmc]{ssmc_summary}}
#'
#' @examples
#'    ## Assess potential new locations
#'    best_results <- best.locale(sites = site.dat[,1:4],
#'    potential.sites = potential.dat,
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
#'
#'    ## Assess existing locations for restoration
#'    best_results <- best.locale(sites = site.dat[,1:4],
#'    restore = TRUE,
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
#'
#' @references Altwegg, R., and H.-U. Reyer. 2003. Patterns of natural selection on size at metamorphosis in water frogs. Evolution 57:872-882.

best.locale <- function(sites,
                        restore = FALSE,
                        potential.sites = NULL,
                        n.sites = 10,
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
                        eps = 1,
                        mu = 2,
                        eta = 0.5,
                        iterations,
                        seed=NULL){

  if(is.null(potential.sites) & isFALSE(restore))
    stop("Must specify either 'potential.sites' OR set 'restore = TRUE' when using this function")

  if(!is.null(seed)){
    set.seed(seed)
  }

  #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
  #### New Site Option ####

  if(isFALSE(restore)) {
    n <- nrow(potential.sites)
    results <- vector("list",n)
    results.df <- vector("list",iterations)
    site.coords <- sites[,2:3]
    pot.sites <- potential.sites
    potential.sites <- potential.sites[,2:3]

    area <- c(sites[,4],median(sites[,4]))
    area[area<=0] <- min(area[area>0])
    colnames(site.coords) <- colnames(potential.sites) <- c("x","y")
    colnames(pop.abun) <- c('mean.ab','sd')

    if(!is.null(met.size)){
      colnames(met.size) <- c('svl.mean','svl.sd')
    }

    progress_bar <- plyr::create_progress_bar("text")
    progress_bar$init(n*iterations)

    for (j in 1:iterations) {
      sample.seed <- sample(1:10e6, 1) + j

      p.philo <- rtnorm(1,prop.philo,sd.philo,lower = lower.upper_philo[[1]], upper = lower.upper_philo[[2]])
      mean.disperse <- rtnorm(1,dispersal,sd.dispersal,lower = lower.upper_dispersal[[1]], lower.upper_dispersal[[2]])

      # Make potential ponds of average population size and metamorph size
      pop.size <- pop.abun %>% rowwise() %>%
        dplyr::mutate(pop.size = rnorm(1,mean.ab,sd))  %>%
        dplyr::mutate(pop.size=ifelse(pop.size < 0,0,pop.size)) %>%
        select(pop.size) %>%
        data.frame(.)

      # Set new pond to mean of all occupied ponds
      pop.size <- round(rbind(pop.size, median(pop.size[pop.size > 0])))


      # Survival--> either based on provided values or result of function
      if(!is.null(met.size)){
        size.dat <- met.size %>%
          rowwise() %>%
          dplyr::mutate(met.size = rnorm(1,svl.mean,svl.sd))  %>%
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
        pop.vec <- pop.size$pop.size
        surv.vec <- p.survive$p.survive

        philo.emig <- pop.size %>%
          dplyr::mutate(survive = rbinom(n(),pop.vec,surv.vec),
                 philo=rbinom(n(),survive,p.philo),
                 emigrant = survive-philo,
                 dead = pop.size - survive)
      } else {
        pop.vec <- pop.size$pop.size

        philo.emig <- pop.size %>%
          dplyr::mutate(survive = rbinom(n(),pop.vec,p.survive),
                 philo=rbinom(n(),survive,p.philo),
                 emigrant = survive-philo,
                 dead = pop.size - survive)
      }

      for(i in 1:n){
        # Combine random potential ponds with existing ponds
        site.dat <- rbind(site.coords,potential.sites[i,])

        # Calculate distance matrix with potential pond added, get connectivity
        dist.mat <- as.matrix(dist(site.dat))
        connect.mat <- exp((-1/mean.disperse)*dist.mat)
        diag(connect.mat) <- 0

        # diag(e.mat) <- philo.emig$philo # Set diagonal to philopatric

        connect.mat2 <- t(apply(connect.mat,1,function(x){x/sum(x)}))

        e.mat <- matrix(NaN,nrow(connect.mat2),nrow(connect.mat2))

        set.seed(sample.seed)
        for (r in 1:nrow(connect.mat2)) {
          e.mat[r,] <- t(rmultinom(1,philo.emig[r,4],connect.mat2[r,]))
        }

        e.mat_orig <- e.mat[-dim(e.mat),-dim(e.mat)]

        colonize <- ifelse(sum(e.mat[-nrow(connect.mat2),nrow(connect.mat2)]) > 2, 1, 0)

        #~#~#~#~#~#~#~#~#~#
        ## METAPOPULATION Mean Lifetime ##
        mmlt <- MMLT(N = dim(e.mat)[1],
                     S_i = colSums(e.mat),
                     S_o = rowSums(e.mat),
                     eps = eps,
                     A = area,
                     mu = mu,
                     eta = eta)

        mmlt_orig <- MMLT(N = dim(e.mat_orig)[1],
                          S_i = colSums(e.mat_orig),
                          S_o = rowSums(e.mat_orig),
                          eps = eps,
                          A = area[-dim(e.mat)],
                          mu = mu,
                          eta = eta)

        mmlt.improve <- ((mmlt - mmlt_orig) / mmlt_orig) * 100

        metrics <- data.frame(row.names = NULL,
                              iter = j,
                              id = i,
                              x = potential.sites[i,1],
                              y = potential.sites[i,2],
                              mmlt = mmlt, # Mean Metapopulation Lifetime
                              mmlt.improve = mmlt.improve,
                              colonize = colonize
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

    mmlt.list <- llply(results.df, function(x) x[,5])
    mmlt <- ldply(mmlt.list, rbind) %>% apply(., 2, mean, na.rm = T)
    rank_mmlt <- rank(-mmlt)

    mmlt_improve.list <- llply(results.df, function(x) x[,6])
    mmlt_improve <- ldply(mmlt_improve.list, rbind) %>% apply(., 2, mean, na.rm = T)

    col.list <- llply(results.df, function(x) x[,7])
    frq_col <- ldply(col.list, rbind) %>% apply(., 2, mean, na.rm = T)
    mmlt_col <- frq_col * mmlt
    rank_adj <- rank(-mmlt_col)
    # lambda <- llply(results.df,function(x) x[,5])
    # lambda_contrib <- llply(results.df,function(x) x[,6])
    # mod <- llply(results.df,function(x) x[,7])

    # lambda_avg <- ldply(lambda,rbind) %>% apply(.,2,mean,na.rm=T)
    # lambda_cont <- ldply(lambda_contrib,rbind) %>% apply(.,2,mean,na.rm=T)
    # mod_avg <- ldply(mod,rbind) %>% apply(.,2,mean,na.rm=T)
    # lambda_cont.sd <- ldply(lambda_contrib,rbind) %>% apply(.,2,sd,na.rm=T)
    # rank_lambda_avg <-  rank(-lambda_avg)
    # rank_lambda_cont <- rank(-lambda_cont)
    # rank_mod <- rank(mod_avg)
    # rank.df <- data.frame(rank_lambda_cont, rank_lambda_avg)
    # avg_rank <- rank(apply(rank.df,1,mean), ties.method = 'min')

    summary.df <- data.frame(site = pot.sites[,1],
                             x =  potential.sites[,1],
                             y =  potential.sites[,2],
                             mmlt = mmlt,
                             rank = rank_mmlt,
                             frq_col = frq_col,
                             mmlt_col = mmlt_col,
                             rank_adj = rank_adj,
                             mmlt_pct.improve = mmlt_improve

                             # mod_avg = mod,
                             # rank_mod = rank_mod,
                             # cont_lam=lambda_cont,
                             # rank_cont=rank_lambda_cont,
                             # lambda=lambda_avg,
                             # rank_lamb=rank_lambda_avg,
                             # avg_rank
    )

    out <- list(summary.df = summary.df,
                results.list = results.df)

    return(out)

    gc() ## Flush memory


    #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
    #### Restore Option ####


    ## If RESTORE is specified
  } else {
    area <- c(sites[,4])
    area[area<=0] <- min(area[area>0])

    colnames(pop.abun) <- c('mean.ab','sd')

    restore.sites <- which(pop.abun$mean.ab==0)
    n <- length(restore.sites)
    results <- vector("list",n)
    results.df <- vector("list",iterations)
    site.coords <- sites[,2:3]
    # pot.sites <- potential.sites
    # potential.sites <- potential.sites[,2:4]

    colnames(site.coords) <- c("x","y")
    # colnames(potential.sites) <- c("x","y","area")

    if(!is.null(met.size)){
      colnames(met.size) <- c('svl.mean','svl.sd')
    }

    progress_bar <- plyr::create_progress_bar("text")
    progress_bar$init(n*iterations)

    for (j in 1:iterations) {
      sample.seed <- sample(1:10e6, 1) + j

      p.philo <- rtnorm(1,prop.philo,sd.philo,lower = lower.upper_philo[[1]], upper = lower.upper_philo[[2]])
      mean.disperse <- rtnorm(1,dispersal,sd.dispersal,lower = lower.upper_dispersal[[1]], lower.upper_dispersal[[2]])

      ## Determine population size from distribution
      pop.size <- pop.abun %>% rowwise() %>%
        dplyr::mutate(pop.size = rnorm(1,mean.ab,sd))  %>%
        dplyr::mutate(pop.size=ifelse(pop.size < 0,0,pop.size)) %>%
        select(pop.size) %>%
        data.frame(.)

      # Set new pond to mean of all occupied ponds
      # mean.pop <- mean(pop.size[pop.size > 0])

      # Set new pond to median of all occupied ponds
      mean.pop <- median(pop.size[pop.size > 0])
      # pop.size <- round(rbind(pop.size, mean(pop.size[pop.size > 0])))


      # Survival--> either based on provided values or result of function
      if(!is.null(met.size)){
        size.dat <- met.size %>%
          rowwise() %>%
          dplyr::mutate(met.size = rnorm(1,svl.mean,svl.sd))  %>%
          select(met.size) %>%
          data.frame(.)

        p.survive <- plogis(-1.366 + (0.87*size.dat[,1])) %>%
          data.frame(.)

        # Make size-based survival at new pond the average of occupied ponds
        mean.survive <- mean(p.survive[p.survive != plogis(-1.366)])
        # p.survive <- (rbind(p.survive, mean(p.survive[p.survive != plogis(-1.366)])))
        colnames(p.survive) <- c('p.survive')

      } else {

        p.survive <- rtnorm(1,prop.survive,sd.survive,lower = lower.upper_survive[[1]], upper = lower.upper_survive[[2]])

        # p.survive <- as.data.frame(p.survive)
        # names(p.survive) <- 'p.survive'
        # mean.survive <- p.survive
      }

      # Calculate distance matrix, get connectivity
      dist.mat <- as.matrix(dist(site.coords))
      connect.mat <- exp((-1/mean.disperse)*dist.mat)
      diag(connect.mat) <- 0

      connect.mat2 <- t(apply(connect.mat,1,function(x){x/sum(x)}))

      e.mat <- matrix(NaN,nrow(connect.mat2),nrow(connect.mat2))

      ## Assess each potential restoration site
      for(i in 1:n){
        set.seed(sample.seed)
        pop.size$pop.size[restore.sites[i]] <- mean.pop

        ## Determine philopatric and emigrant populations after accounting for survival
        if(!is.null(met.size)){
          p.survive$p.survive[restore.sites[i]] <- mean.survive

          pop.vec <- pop.size$pop.size
          surv.vec <- p.survive$p.survive

          philo.emig <- round(pop.vec) %>% data.frame() %>%
            dplyr::mutate(survive = rbinom(n(),round(pop.vec),surv.vec),
                   philo=rbinom(n(),survive,p.philo),
                   emigrant = survive-philo,
                   dead = pop.size - survive)
        } else {
          pop.vec <- pop.size$pop.size

                    philo.emig <- round(pop.vec) %>%  data.frame() %>%
            dplyr::mutate(survive = rbinom(n(),round(pop.vec),p.survive),
                   philo=rbinom(n(),survive,p.philo),
                   emigrant = survive-philo,
                   dead = pop.size - survive)
        }



        for (r in 1:nrow(connect.mat2)) {
          e.mat[r,] <- t(rmultinom(1,philo.emig[r,4],connect.mat2[r,]))
        }

        e.mat_orig <- e.mat
        e.mat_orig[restore.sites[i],] <- 0

        colonize <- ifelse(colSums(e.mat)[restore.sites[i]] > 2, 1, 0)

        #~#~#~#~#~#~#~#~#~#
        ## METAPOPULATION Mean Lifetime ##
        mmlt <- MMLT(N = dim(e.mat)[1],
                     S_i = colSums(e.mat),
                     S_o = philo.emig$emigrant,
                     eps = eps,
                     A = area,
                     mu = mu,
                     eta = eta)

        mmlt_orig <- MMLT(N = dim(e.mat_orig)[1],
                          S_i = colSums(e.mat_orig),
                          S_o = rowSums(e.mat_orig),
                          eps = eps,
                          A = area,
                          mu = mu,
                          eta = eta)

        mmlt.improve <- ((mmlt - mmlt_orig) / mmlt_orig) * 100

        metrics <- data.frame(row.names = NULL,
                              iter = j,
                              id = sites[restore.sites[i],1],
                              x = sites[restore.sites[i],2],
                              y = sites[restore.sites[i],3],
                              mmlt = mmlt, # Mean Metapopulation Lifetime
                              mmlt.improve = mmlt.improve,
                              colonize = colonize
        )



        results[[i]] <- metrics

        progress_bar$step()

      } # Close potential site loop

      results.df[[j]] <- ldply(results,"identity")

    } # End iterations loop

    mmlt.list <- llply(results.df, function(x) x[,5])
    mmlt <- ldply(mmlt.list, rbind) %>% apply(., 2, mean, na.rm = T)
    rank_mmlt <- rank(-mmlt)

    mmlt_improve.list <- llply(results.df, function(x) x[,6])
    mmlt_improve <- ldply(mmlt_improve.list, rbind) %>% apply(., 2, mean, na.rm = T)

    col.list <- llply(results.df, function(x) x[,7])
    frq_col <- ldply(col.list, rbind) %>% apply(., 2, mean, na.rm = T)
    mmlt_col <- frq_col * mmlt
    rank_adj <- rank(-mmlt_col)


    summary.df <- data.frame(site = sites[restore.sites,1],
                             x =  sites[restore.sites,2],
                             y =  sites[restore.sites,3],
                             mmlt = mmlt,
                             rank = rank_mmlt,
                             frq_col = frq_col,
                             mmlt_col = mmlt_col,
                             rank_adj = rank_adj,
                             mmlt_pct.improve = mmlt_improve
    )

    out <- list(summary.df = summary.df,
                results.list = results.df)
    return(out)
    gc()  ## Flush memory
  }
} # End function
