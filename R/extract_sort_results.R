#'  A function to extract the summary data from \code{\link[ssmc]{site.analysis}} or \code{\link[ssmc]{best.locale}} functions
#'
#'  @param ssmc_results Output object from \code{\link[ssmc]{site.analysis}} or \code{\link[ssmc]{best.locale}}
#'
#'  @usage ssmc_summary(ssmc_results)
#'
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @return This function returns a data frame of the summarized results from either the \code{\link[ssmc]{site.analysis}} or \code{\link[ssmc]{best.locale}} sorted by the specified column.
#'
#' @export
#'

#'
#'    @usage ssmc_summary(ssmc_results)
#'
#'    @examples
#'    # Assess potential new locations
#'    best_results <- best.locale(sites = site.dat[,1:3],
#'    potential.sites = potential.dat,
#'    pop.abun = site.dat[,4:5],
#'    met.size = site.dat[,6:7],
#'    prop.philo = 0.95,
#'    sd.philo = 0.05,
#'    min.max_philo=c(min=0, max=1),
#'    # prop.survive=0.2,      ## Not needed b/c met.size is specified
#'    # sd.survive = 0.075,    ## Not needed b/c met.size is specified
#'    min.max_survive = c(min=0, max=1),
#'    dispersal = 25,
#'    sd.dispersal = 10,
#'    min.max_dispersal = c(min=10, max=Inf),
#'    iterations = 10,
#'    seed = 123)
#'
#'    best_locales <- ssmc_summary(best_results)

ssmc_summary <- function(ssmc_results, ...) {
  dots <- list(...)                   #1
  if(length(dots) == 0) {
    sort.by <- 'rank'
  } else {
    sort.by <- dots[[1]]
  }

  out <- arrange_(ssmc_results$summary.df, sort.by)
  return(out)
}
