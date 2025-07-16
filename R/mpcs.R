
# cosine similarity helper functions
# 1: fast cosine similarity calculation between vector and vector
# 2: fast cosine similarity calculation between matrix and vector (with max by default)
vectorCosSim <- function(v1, v2){
  z <- v1 %*% v2
  n <- sqrt(crossprod(v1)*crossprod(v2))
  return( z/n )
}
matrixCosSim <- function(v, m, max=T){
  z <- m %*% v
  n <- sqrt(apply(m, 1, crossprod) %*% crossprod(v))
  
  if(max) return( max(z/n) )
  else return( z/n )
}

# precursor to the maximum pairwise cosine similarity 
# (one direction = not symmetric, m1 -> m2)
preMPCS <- function(m1, m2, agg="mean"){
  # convert matrices (or vectors) to a list
  if(is.vector(m1)) m1 <- t(m1)
  if(is.vector(m2)) m2 <- t(m2)
  
  # check if any row is only 0s, should not be possible
  if( any(rowSums(m1)==0) | any(rowSums(m2)==0) ) return()
  
  # check if the permutations have the same length (symptoms)
  if(ncol(m1) != ncol(m2)) return()
  
  # maximum cos sim from each row of m1 to each row of m2
  sim <- sapply(1:nrow(m1), function(i) matrixCosSim(m1[i,], m2))
  agg_sim <- switch(agg,
             "mean"= mean(sim),
             "med" = median(sim),
             "max" = max(sim),
             "min" = min(sim),
             "topx"= mean( sort(sim[1:(floor(nrow*0.05))], decreasing=T) ),
             "This aggregate doesnt exist for this function! (DPS)"
  ) # wrong name input case at the end; "med", "min", "topx" unused/experimental
  
  return(agg_sim)
}

# maximum pairwise cosine similarity (symmetric)
MPCS <- function(m1, m2, agg="mean"){
  sym_sim <- max(preMPCS(m1, m2, agg), preMPCS(m2, m1, agg))
  return(sym_sim)
}

# cosine similarity as a distance object (=> cosine distance/angular distance)
distCosSim <- function(m, ang=T, distObj=T){
  if(is.vector(m)) m <- t(m)
  
  # check if any row is only 0s, should not be possible
  # if( any(rowSums(m1)==0) | any(rowSums(m2)==0) ) return()
  
  sims <- sapply(1:nrow(m), function(i) matrixCosSim(m[i,], m, max=F))
  
  if(ang) sims <- 2*acos(sims)/pi #angular distance
  else    sims <- 1-sims #cosine distance
  
  if(distObj) return( as.dist(sims) )
  else return( sims )
}
