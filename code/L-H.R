## ==========================================================================================
## Implementation of the Lemke-Howson algorithm for symmetric games using a Tableaux method,
## following  https://github.com/s3rvac/lemke-howson
## ==========================================================================================

makePivotingStep <- function(t, n, e){
  ne = ifelse(e > 0, e + 2, abs(e) + n + 2)
  L = get_leaving_variable(t, ne)
  l = L[["l"]]
  lindex = L[["index"]]
  nl = ifelse(l > 0, l + 2, abs(l) + n + 2)
  t = update_tableux(t, e, ne, nl, lindex)
  return (list("left" = l, "t" = t))
}


lemkeHowson_symmetric <- function(M){
  ## Before we start, we need to normalize both matrices
  ## to ensure some assumptions about values in both matrices
  normM <- normalizeMatrix(M)
  ## Create the tableaux that will be used in the pivoting procedure
  t <- createTableaux(normM)
  
  ## Make pivoting steps until the equilibrium is found
  ## (the variable that left the basis is the same (in absolute value)
  ## as the variable that we used as an initial pivot)
  n = nrow(normM)
  init = n
  U = makePivotingStep(t, n, init)
  left = U[["left"]]
  t = U[["t"]]
  while (left != -init){
    U = makePivotingStep(t, n, -left)
    left = U[["left"]]
    t = U[["t"]]
  }    
  ## Get the equilibrium from the resulting tableaux,
  ## normalize it and return it
  eq <- getEquilibrium(t, n)
  S <- t[,1][t[,1] > 0]
  S <- S[order(S)]
  S <- S[-length(S)]
  return (list(eq = eq / sum(eq), subset = S))
}


#### get the final equilibrium
getEquilibrium <- function(t, n){
  eq <- rep(0, n)
  for (i in 1:n){
    p <- t[i,1]
    if (p > 0){
      eq[p] <- t[i, 2]
    }    
  }
  return (eq)
}

### Convert from Replicator to Lotka-Volterra
toLV <- function(x, n){
  return (x[1:n]/x[n+1])
}


### normalize the matrix to have all positive values
normalizeMatrix <- function(M){
  ## Check for the least value in both matrices
  lowestVal <- min(M)
  n <- nrow(M)
  normM <- t(t(M) + seq(1, n))
  
  if (lowestVal <= 0){
    normM <- normM + abs(lowestVal)
  }
  
  return (normM)
}


### Initialize the tableaux 
createTableaux <- function(M){
  S <- nrow(M)
  N <- 2 * S
  t <- matrix(0, S, N + 2)
  basic <- seq(-1, -S, -1)
  t[, 1] <- basic
  t[, 2] <- 1
  t[, 3:(S + 2)] <- -M
  return (t)
}

### Get the leaving variable following a min rule
get_leaving_variable <- function(t, e){
  n <- nrow(t)
  vals <- ifelse(t[,e] >=0, Inf, -t[,2]/t[,e])
  index <- which.min(vals)
  l <- t[index, 1]
  
  return (list("l" = l, "index" = index))
}

### Update tableaux after pivoting step
update_tableux <- function(t, e, ne, l, lindex){
  t[lindex, 1] <- e
  n <- nrow(t)
  m <- ncol(t)
  ## subset of no zero coefficients
  P <- as.logical((t[,ne]!=0) * (seq(1,n) != lindex))
  t[lindex, l] <- -1
  target <- t[lindex, 2:m] / -t[lindex, ne]
  t[lindex, 2:m] <- target
  c <- t[P, ne]
  W <- c %*% t(target)
  t[P, 2:m] <- t[P, 2:m] + W
  t[P, ne] <- 0
  t[lindex, ne] <- 0
  return (t)
}

## =================================================
## Translate from LV to Replicator dynamics
## =================================================

build.rep <- function(A, r, n){
  M <- matrix(0, n+1, n+1)
  M[1:n, n+1] <- r
  M[1:n, 1:n] <- A
  return (M)
}

get_final_composition <- function(A, r){
  # A is symmetric and D-Stable
  # r is a vector of growth rates
  # return the final composition of the community using 
  # Lemke-Howson algorithm
  n <- nrow(A)
  M <- build.rep(A, r, n)                                     ## set payoff matrix for replicator dynamics
  
  ### Run dynamics
  y <- lemkeHowson_symmetric(M)  ## get nash equilibrium of the game
  S <- y$subset                  ## Subset of coexisting species
  x <- toLV(y$eq, n)          ## equilibrium of LV model
  return(x)
}