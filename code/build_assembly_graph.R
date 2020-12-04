source("general_code.R")
source("L-H.R")
n <- 5
A <- -matrix(runif(n * n), n, n)
# make symmetric
A <- A + t(A)
# positive growth rates
r <- runif(n)

# strategy: 
# 1) go through each possible state and determine
# a) is it feasible/stable
# b) how many species
# c) if not feasible/stable, what it falls down to
# 
# 2) build the graph
# a) go through every stable state and determine their neighbor
# b) if the neighbor is feasible/stable, add the edge
# c) if not, add an edge to the state the neighbor collapses to


n <- nrow(A)
states <- data.frame(labels = 0:(2^n - 1), 
                 stable = 0, 
                 num_species = 0,
                 attractor = NA
                 )

# the first state is the empty state
states[1, ] <- c(0, 1, 0, 0)
for (i in 2:(2^n)){
  state <- as.numeric(intToBits(i - 1)[1:n])
  # single-species states
  if (sum(state) == 1){
    my_spp <- which(state == 1)
    if ((r[my_spp] > 0) & (A[my_spp, my_spp] < 0)){
      states[i,] <- c(i-1, 1, 1, i-1)
    } else {
      states[i,] <- c(i-1, 0, 1, 0)
    }
  } else {
    present <- state > 0
    Ared <- A[present, present, drop = FALSE]
    rred <- r[present]
    to_state <- (get_final_composition(Ared, rred) > 0) * 1
    attractor <- sum(2^to_state)
    if (attractor == i-1){
      states[i,] <- c(i-1, 1, sum(state), i-1)
    } else {
      states[i,] <- c(i-1, 0, sum(state), attractor)
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
    possible_neighbors <- states[states$num_species == my_num_spp + 1,]
    if (nrow(possible_neighbors) > 0){
      for (j in 1:nrow(possible_neighbors)){
        # check if it is really a neighbor
        if (bitwAnd(possible_neighbors$labels[j], my_label) == my_label){
          # check if it is a stable state
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
edges <- unique(edges)
edges <- edges[edges[,1] != edges[,2], ]

library(igraph)
# build igraph
gg <- graph_from_edgelist(cbind(as.character(edges[,1]),as.character(edges[,2])), directed = TRUE)
# level
