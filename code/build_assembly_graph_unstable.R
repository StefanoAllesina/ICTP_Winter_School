build_competitive_unstable <- function(n){
  A <- matrix(0, n, n)
  A[upper.tri(A)] <- -2 * abs(rnorm(n * (n - 1) / 2))
  # add a diagonal
  diag(A) <- -abs(rnorm(n))
  # make symmetric
  A <- A + t(A)
  # now find the largest eigenvalue
  l1A <- max(eigen(A, only.values = TRUE, symmetric = TRUE)$values)
  if (l1A < 0){
    # shift the diagonal to make A unstable
    diag(A) <- diag(A) - (l1A + 1)
  }
  # only return competitive matrices
  if (all(diag(A) <= 0)) return(A)
  # if it failed, call recursively
  return(build_competitive_unstable(A))
}

plot_assembly_graph <- function(graph, info, plot_biom = FALSE){
  # this works well if the graph is a DAG
  coords <- layout_with_sugiyama(graph, 
                                 attributes="all", 
                                 hgap=10, vgap=10)$layout
  if (plot_biom) {
    V(graph)$color <- grey(info$biomass / max(info$biomass))
  } else {
    V(graph)$color <- "SkyBlue"
  }
  plot(graph, layout = coords, vertex.label.color = "orange")
}

build_assembly_graph_unstable <- function(r, A, num_invasions = 1){
  n <- nrow(A)
  labels <- 0:(2^n-1) 
  # states contains 1) labels, 2) feas/stable, 3) num_species, 4) label it collapses to, 5) biomass/energy
  states <- matrix(0, 2^n, 5)
  states[,1] <- labels
  # step 1: check if feasible/stable, count species, and determine biomass/energy
  for (i in 1:nrow(states)){
    my_label <- labels[i]
    my_collapse <- NA
    my_energy <- NA
    my_stable <- NA
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
        is_stable <- NA
        r_comm <- r[presence]
        A_comm <- A[presence, presence]
        xstar <- solve(A_comm, -r_comm)
        if (all(xstar > 0)){
          # feasible
          is_stable <- (max(eigen(diag(xstar) %*% A_comm, 
                           only.values = TRUE)$values) < 0)
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
    }
    states[i,] <- c(my_label, my_stable, my_num_species, my_collapse, my_energy)
  }
  # step 2: for every unstable community, select the subcommunity with the highest energy
  for (i in 1:nrow(states)){
    # if it will collapse
    if(is.na(states[i,4])){
      my_label <- states[i, 1]
      # find all subcommunities
      is_subcomm <- bitwAnd(labels, my_label) == labels
      energy_subcomm <- is_subcomm * states[,5]
      # the subcommunity it collapses to is the one with the highest energy
      states[i, 4] <- labels[which.max(energy_subcomm)]
    }
  }
  # step 3: build the graph
  states <- as.data.frame(states)
  colnames(states) <- c("labels", "stable", "num_species", "attractor", "biomass")
  edges <- matrix(0, 0, 2)
  # for each node, determine neighbors
  for (i in 1:nrow(states)){
    my_num_spp <- states$num_species[i]
    my_label <- states$labels[i]
    my_stable <- states$stable[i]
    if (my_stable == 1){
      # compute x_bar
      x_bar <- rep(0, n)
      presence <- (intToBits(my_label)[1:n]) > 0
      if (sum(presence) > 0){
        r_comm <- r[presence]  
        A_comm <- A[presence, presence, drop = FALSE]
        xstar <- solve(A_comm, -r_comm)
        x_bar[presence] <- xstar
      }
      per_capita <- r + A %*% x_bar
      possible_neighbors <- states[states$num_species <= my_num_spp + num_invasions,]
      if (nrow(possible_neighbors) > 0){
        for (j in 1:nrow(possible_neighbors)){
          # check if it is really a neighbor
          if (possible_neighbors$labels[j] != my_label){
            if (bitwAnd(possible_neighbors$labels[j], my_label) == my_label){
              # if it can invade
              # compute per-capita growth rate
              to_check <- bitwXor(my_label, possible_neighbors$labels[j])
              invasion_growth_rates <- per_capita[(intToBits(to_check)[1:n]) > 0]
              if (all(invasion_growth_rates > 0)) {
                if (possible_neighbors$stable[j] == 1){
                  edges <- rbind(edges, c(my_label, possible_neighbors$labels[j]))
                } else {
                  edges <- rbind(edges, c(my_label, possible_neighbors$attractor[j]))
                }
              }
            }
          }
        }
      }
    }
  }
  edges <- unique(edges)
  edges <- edges[edges[,1] != edges[,2], ]
  gg <- graph_from_edgelist(cbind(as.character(edges[,1]),as.character(edges[,2])), directed = TRUE)
  
  # extract info on states
  info <- states[as.numeric(names(V(gg))) + 1, ]
  return(list(graph = gg, info = info))
}
