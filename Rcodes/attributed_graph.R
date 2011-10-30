library("evd")

kidney.egg.attributed <- function(n,m,pi0,piA){
  
  pi00 <- pi0*pi0
  pi0A <- pi0*piA; piA0 <- pi0A
  piAA <- piA*piA
  
  runif.mat <- matrix(0,nrow = n, ncol = n)
  runif.mat[col(runif.mat) > row(runif.mat)] <- runif(n*(n-1)/2)
  runif.mat <- (runif.mat + t(runif.mat))

  edge.mat <- matrix(0, nrow = n, ncol =n)
  edge.mat[1:m,1:m] <- sum(piAA)
  edge.mat[1:m,(m+1):n] <- sum(pi0A) 
  edge.mat[(m+1):n,1:m] <- sum(pi0A)
  edge.mat[(m+1):n,(m+1):n] <- sum(pi00)

  red.mat <- matrix(0, nrow = n, ncol = n)
  red.mat[1:m,1:m] <- piAA[1]
  red.mat[1:m,(m+1):n] <- pi0A[1]
  red.mat[(m+1):n,1:m] <- pi0A[1]
  red.mat[(m+1):n,(m+1):n] <- pi00[1]

  adjacency.mat <- matrix(0,nrow = n, ncol = n)
  red.idx <- (runif.mat <= red.mat)
  green.idx <- (runif.mat > red.mat & runif.mat <= edge.mat)

  adjacency.mat[red.idx] <- 1
  adjacency.mat[green.idx] <- 2
  diag(adjacency.mat) <- 0

  vertex.colors <- c(rep("red",length.out=m),
                     rep("green",length.out=(n-m)))

  return(list(adjacency=adjacency.mat,
              v.colors = vertex.colors))
}

maxd.statistic <- function(x,g){
  t1 <- rowSums(g$adjacency == 1)
  t2 <- rowSums(g$adjacency == 2)
  t <- t1*x[1] + t2*x[2]
  return(max(t))
}

num.triangle.trace <- function(w, g){
  A <- g$adjacency
  A1 <- (A == 1) + 0
  A2 <- (A == 2) + 0
  A.weighted <- w[1]*A1 + w[2]*A2

  return(sum(diag(A.weighted %*% A.weighted %*% A.weighted))/6)
}

scan.statistic <- function(x,g){
  x.red <- x[1]
  x.green <- x[2]
  nvertex <- nrow(g$adjacency)
  locality.vec <- seq(0,0,length.out = nvertex)
  for(i in 1:nvertex){
    i.neighbours <- which(g$adjacency[i,] %in% c(1,2))
    idx <- c(i,i.neighbours)
    s.red <- sum(g$adjacency[idx,idx] == 1)
    s.green <- sum(g$adjacency[idx,idx] == 2)
    locality.vec[i] <- (s.red*x.red + s.green*x.green)/2
  }
  return(max(locality.vec))
}

maxd.mc <- function(theta, nvertex, nvertex.abnormal, pi0, piA, mc){
  maxd.vec <- numeric(mc)
  x <- c(cos(theta), sin(theta))
  for(i in 1:mc){
    g.i <- kidney.egg.attributed(nvertex, nvertex.abnormal, pi0, piA)
    maxd.vec[i] <- maxd.statistic(x,g.i)
  }
  return(mean(maxd.vec))
}

maxd.first.approximation <- function(nvertex, nvertex.abnormal, pi0, piA, mc = 1000){
  theta.seq <- seq(from = 0, to = pi/2, by = 0.01)
  maxd.vec <- numeric(length(theta.seq))
  for(i in 1:length(theta.seq)){
    maxd.vec[i] <- maxd.mc(theta.seq[i], nvertex, nvertex.abnormal, pi0, piA, mc)
  }
  return(maxd.vec)
}
    
scan.mc <- function(theta,nvertex, nvertex.abnormal, pi0, piA, mc){
  scan.vec <- numeric(mc)
  x <- c(cos(theta), sin(theta))
  for(i in 1:mc){
    g.i <- kidney.egg.attributed(nvertex, nvertex.abnormal, pi0, piA)
    scan.vec[i] <- scan.statistic(x, g.i)
  }
  return(scan.vec)
}
