A <- -matrix(runif(100), 10, 10)
A <- A + t(A)
r <- runif(10)

# cycle through all possible compositions and check whether it is a node in the graph
n <- nrow(A)
states <- tibble(label = 0, state = list(rep(0, n)))
edges <- tibble(from = numeric(0), to = numeric(0))
for (i in 1:(2^n - 1)){
  label <- i
  state <- as.numeric(intToBits(i)[1:n])
  # single-species states
  if (sum(state) == 1){
    my_spp <- which(state == 1)
    if ((r[my_spp] > 0) & (A[my_spp, my_spp] < 0)){
      states <- bind_rows(states, tibble(label = i, state = list(state)))
      edges <- bind_rows(edges, tibble(from = i, to = i))
    } else {
      edges <- bind_rows(edges, tibble(from = i, to = 0))
    }
  } else {
    present <- state > 0
    Ared <- A[present, present, drop = FALSE]
    rred <- r[present]
    to_state <- (get_final_composition(Ared, rred) > 0) * 1
    if (all(to_state > 0)){
      states <- bind_rows(states, tibble(label = i, state = list(state)))
    }
    edges <- bind_rows(edges, tibble(from = i, to = sum(2^to_state)))
  }
}
