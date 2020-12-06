build_competitive_stable <- function(n){
  A <- matrix(0, n, n)
  A[upper.tri(A)] <- -abs(rnorm(n * (n - 1) / 2))
  # add a diagonal
  diag(A) <- -abs(rnorm(n))
  # make symmetric
  A <- A + t(A)
  # now find the largest eigenvalue
  l1A <- max(eigen(A, only.values = TRUE, symmetric = TRUE)$values)
  if (l1A > 0){
    # shift the diagonal to make A stable
    diag(A) <- diag(A) - (l1A * 1.1)
  }
  return(A)
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


build_assembly_graph <- function(r, A, num_invasions = 1){
  n <- nrow(A)
  # strategy: 
  # 1) go through each possible state and determine
  # a) is it feasible/stable
  # b) how many species
  # c) if not feasible/stable, what it falls down to
  # d) what is the biomass
  # 
  # 2) build the graph
  # a) go through every stable state and determine their neighbor
  # b) if the neighbor is feasible/stable, add the edge
  # c) if not, add an edge to the state the neighbor collapses to
  
  
  n <- nrow(A)
  labels <- 0:(2^n - 1)
  sp_base2 <- 2^(0:(n-1))
  empty <- rep(0, n)
  states <- data.frame(labels = labels, 
                   stable = 0, 
                   num_species = 0,
                   attractor = NA,
                   biomass = NA
                   )
  
  # the first state is the empty state
  states[1, ] <- c(0, 1, 0, 0, 0)
  for (i in 2:(2^n)){
    state <- as.numeric(intToBits(i - 1)[1:n])
    # single-species states
    if (sum(state) == 1){
      my_spp <- which(state == 1)
      if ((r[my_spp] > 0) & (A[my_spp, my_spp] < 0)){
        states[i,] <- c(i-1, 1, 1, i-1, -r[my_spp]/A[my_spp, my_spp])
      } else {
        states[i,] <- c(i-1, 0, 1, 0, NA)
      }
    } else {
      present <- state > 0
      Ared <- A[present, present, drop = FALSE]
      rred <- r[present]
      x_star <- get_final_composition(Ared, rred)
      x_bar <- empty
      x_bar[present] <- x_star
      to_state <- sum(sp_base2 * (x_bar > 0))
      if (to_state == i-1){
        states[i,] <- c(i-1, 1, sum(state), i-1, sum(x_bar * r))
      } else {
        states[i,] <- c(i-1, 0, sum(state), to_state, NA)
      }
    }
  }
  
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
      possible_neighbors <- states[states$num_species <= my_num_spp + num_invasions, ]
      if (nrow(possible_neighbors) > 0){
        for (j in 1:nrow(possible_neighbors)){
          # check if it is really a neighbor
          if (possible_neighbors$labels[j] != my_label){
            if (bitwAnd(possible_neighbors$labels[j], my_label) == my_label){
              # if it can invade
              # compute per-capita growth rate
              per_capita <- r + A %*% x_bar
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

