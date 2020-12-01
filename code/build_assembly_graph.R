source("general_code.R")
source("L-H.R")
n <- 5
A <- -matrix(runif(n * n), n, n)
# make symmetric
A <- A + t(A)
# positive growth rates
r <- runif(n)

n <- nrow(A)
edges <- tibble(from = numeric(0), to = numeric(0))
for (i in 1:(2^n - 1)){
  label <- i
  state <- as.numeric(intToBits(i)[1:n])
  # single-species states
  if (sum(state) == 1){
    my_spp <- which(state == 1)
    if ((r[my_spp] > 0) & (A[my_spp, my_spp] < 0)){
      edges <- bind_rows(edges, tibble(from = 0, to = i))
      edges <- bind_rows(edges, tibble(from = i, to = i))
    } else {
      edges <- bind_rows(edges, tibble(from = i, to = 0))
    }
  } else {
    present <- state > 0
    Ared <- A[present, present, drop = FALSE]
    rred <- r[present]
    to_state <- (get_final_composition(Ared, rred) > 0) * 1
    edges <- bind_rows(edges, tibble(from = i, to = sum(2^to_state)))
  }
}

# for each node, determine who is one invasion above

