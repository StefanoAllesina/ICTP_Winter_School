build_competitive_unstable <- function(n){
  A <- matrix(0, n, n)
  A[upper.tri(A)] <- -abs(rnorm(n * (n - 1) / 2))
  # add a diagonal
  diag(A) <- -abs(rnorm(n))
  # make symmetric
  A <- A + t(A)
  # now find the largest eigenvalue
  l1A <- max(eigen(A, only.values = TRUE, symmetric = TRUE)$values)
  if (l1A < 0){
    # shift the diagonal to make A stable
    diag(A) <- diag(A) + (l1A * 1.1)
  }
  # only return competitive matrices
  if (all(diag(A) <=0)) return(A)
  # if it failed, call recursively
  return(build_competitive_unstable(A))
}

A <- build_competitive_unstable(4)
r <- runif(4)

#build_assembly_graph_unstable <- function(r, A, num_invasions = 1){
  n <- nrow(A)
  labels <- 0:(2^n-1) 
  # states contains 1) labels, 2) feas/stable, 3) num_species, 4) label it collapses to, 5) biomass/energy
  states <- matrix(0, 2^n, 5)
  states[,1] <- labels
  # step 1: check if feasible/stable, count species, and determine biomass/energy
  for (i in 1:nrow(states)){
    my_label <- labels[i]
    presence <- (intToBits(my_label)[1:n]) > 0
    my_num_species <- sum(presence)
    # empty state, treat separately
    if (my_num_species == 0){
      my_stable <- 1
      my_collapse <- 0
      my_energy <- 0
    } else {
      # single species, treat separately
      if (my_num_species == 1){
        my_species <- which(presence == TRUE)
        if (r[my_species] > 0){
          # the species can invade the empty state
          my_stable <- 1
          my_energy <- -r[my_species]^2 / A[my_species, my_species]
          my_collapse <- my_label
        } else {
          my_stable <- 0
          my_energy <- 0
          my_collapse <- 0
        }
      } else {
        # multispecies community
        my_stable <- 0
        r_comm <- r[presence]
        A_comm <- A[presence, presence]
        xstar <- solve(A_comm, -r_comm)
        if (all(xstar > 0)){
          # feasible
          is_stable <- (max(eigen(diag(xstar) %*% A_comm, 
                           symmetric = TRUE, only.values = TRUE)$values) < 0) 
          if (is_stable){
            # also stable
            my_stable <- 1
            my_energy <- sum(r_comm * xstar)
            my_collapse <- my_label
          } else {
            # unstable---need to find the collapse
            my_stable <- 0
            my_energy <- 0
            my_collapse <- NA
          }
        } else {
          # unfeasible
          my_stable <- 0
          my_energy <- 0
          my_collapse <- NA
        }
      }
      states[i,] <- c(my_label, my_stable, my_num_species, my_collapse, my_energy)
    }
  }
  # step 2: for every unstable community, select the subcommunity with the highest energy
  for (i in 1:nrow(states)){
    # if it will collapse
    if(is.na(states[i,4])){
      my_label <- states[i, 1]
      # find all subcommunities
      is_subcomm <- bitwAnd(labels, my_label) == labels
      print(my_label)
      print(is_subcomm)
      energy_subcomm <- is_subcomm * states[,5]
      # the subcommunity it collapses to is the one with the highest energy
      states[i,4] <- labels[which.max(energy_subcomm)]
    }
  }
#}