# import libraries
library(igraph) # graph visualization
library(tidyverse) # data manipulation, plotting
library(deSolve) # integrate differential equations
# constants
THRESH <- 10^(-8) # consider a species extinct if it falls below

# Generalized Lotka-Volterra model
GLV <- function(t, z, params){
  # t is the time, z is the density of all species
  # params contains the interaction matrix A and the growth rates r
  r <- params$r
  A <- params$A
  # remove extinct species
  z[z < THRESH] <- 0
  # compute growth
  dzdt <- diag(z) %*% (r + A %*% z)
  return(list(dzdt))
}

# Interface to integrate GLV; computes equilibria and their stability
GLV_dynamics <- function(z0 = NULL, labels = NULL, A, r, maxtime = 1000, bytime = 0.1){
  # z0: initial conditions
  # labels: species names
  # A, r: model parameters
  # maxtime, bytime: total integration time and time step
  n <- nrow(A)
  # if z0 not provide, sample fron U[0,1]
  if (is.null(z0)) z0 <- runif(n)
  # if labels are not provided, use generic ones
  if (is.null(labels)) labels <- paste("x", 1:n, sep = "")
  # now integrate
  times <- seq(0, maxtime, by = bytime)
  p <- list(A = A, r = r)
  ts <- ode(y = z0, func = GLV, times = times, parms = p, method = "ode45")
  # make time series into tidy form for plotting
  ts <- as.matrix(as.data.frame(ts))
  ts[ts < THRESH] <- 0
  # if the system has diverged...
  ts[is.infinite(ts)] <- 0
  ts[is.nan(ts)] <- 0
  # pruned community
  x_t <- ts[nrow(ts), -1]
  colnames(ts) <- c("time", labels)
  ts <- ts %>% as_tibble() %>% 
    pivot_longer(-time, names_to = "species", values_to = "density")
  # compute equilibria and stability
  # a) of the whole community
  xs <- solve(A, -r)
  feasible <- all(xs > 0)
  stable <- FALSE
  if (feasible) stable <- max(Re(eigen(diag(xs, nrow = length(xs)) %*% A, 
                                       only.values = TRUE)$values)) < 0
  # b) of the pruned community 
  pruned_A <- A[x_t > 0, x_t > 0, drop = FALSE]
  pruned_r <- r[x_t > 0]
  pruned_xs <- NA
  pruned_feasible <- FALSE
  pruned_stable <- FALSE
  if (length(pruned_r) > 0){
    pruned_xs <- solve(pruned_A, -pruned_r)
    pruned_feasible <- all(pruned_xs > 0)
    pruned_stable <- NA
    if (pruned_feasible) pruned_stable <- max(Re(eigen(diag(pruned_xs, nrow = length(pruned_xs)) %*% 
                                                       pruned_A, only.values = TRUE)$values)) < 0
    pruned_xs <- tibble(species = labels[x_t > 0], density = pruned_xs)
  }
  # make equilibria in tidy for for plotting
  xs <- tibble(species = labels, density = xs)
  # return a list with all the info
  return(list(
   n = n,
   r = r, 
   A = A,
   xs = xs,
   feasible = feasible,
   stable = stable,
   
   pruned_r = pruned_r, 
   pruned_A = pruned_A,
   pruned_xs = pruned_xs,
   pruned_feasible = pruned_feasible,
   pruned_stable = pruned_stable,
   
   ts = ts
  ))
}

# Integrate the model with z0; at time 0 introduce the invasion by summing
invasion_dynamics <- function(z0, invasion, A, r, timebefore = 1000, timeafter = 1000, bytime = 0.1){
  before <- GLV_dynamics(z0 = z0, labels = NULL, A = A, r = r, maxtime = timebefore, bytime = bytime)
  n <- length(r)
  state_at_invasion <- before$ts$density[(nrow(before$ts) - n+1):nrow(before$ts)]
  z_inv <- state_at_invasion + invasion
  after <- GLV_dynamics(z0 = z_inv, labels = NULL, A = A, r = r, maxtime = timebefore, bytime = bytime)
  after$ts$time <- after$ts$time + timebefore
  after$ts <- rbind(before$ts, after$ts)
  return(after)
}

# Plot the time series
plot_dynamics <- function(dyn, add_feasible_equilibria = TRUE){
  pl <- ggplot(dyn$ts) + 
    aes(x = time, y = density, colour = species) + 
    geom_line() + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    scale_y_sqrt()
  if (add_feasible_equilibria) {
    pl <- pl + geom_hline(data = dyn$pruned_xs, 
                          aes(yintercept = density, colour = species),
                          linetype = 2)
  }
  return(pl)
}