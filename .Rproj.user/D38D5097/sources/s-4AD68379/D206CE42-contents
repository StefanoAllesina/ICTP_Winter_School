# find params for which we have either limit cycles or chaos
source("general_code.R")

find_next_interesting_param <- function(i){
  success <- FALSE
  while(!success){
    i <- i + 1
    print(i)
    set.seed(i)
    n <- 5
    A <- matrix(rnorm(n * n), n, n)
    r <- rnorm(n)
    diag(A)[r < 0] <- 0
    diag(A) <- -abs(diag(A))
    tmp <- GLV_dynamics(A = A, r = r, maxtime = 500)
    if (tmp$pruned_feasible == TRUE){
      if (tmp$pruned_stable == FALSE){
        success <- TRUE
        show(plot_dynamics(tmp))
      }
    }
  }
  return(list(r = tmp$pruned_r, A = tmp$pruned_A))
}

limit_cycle_4spp <- find_next_interesting_param(59)
save(limit_cycle_4spp, file = "four_spp_cycle.RData")

limit_cycke_3spp <- find_next_interesting_param(150)
tmp <- find_next_interesting_param(620)
