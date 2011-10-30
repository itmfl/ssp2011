kidney.egg.attributed <- function(n,m,theta,pi0,piA){
  
  pi00 <- pi0*pi0
  pi0A <- pi0*piA
  piA0 <- pi0A
  piAA <- piA*piA
  x.red <- cos(theta)
  x.green <- sin(theta)
  
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
  edge.idx <- (runif.mat < edge.mat)
  adjacency.mat[edge.idx] <- 1

  labels.mat <- matrix("none", nrow = n, ncol = n)
  red.idx <- (runif.mat <= red.mat)
  green.idx <- (runif.mat > red.mat & runif.mat <= edge.mat)

  labels.mat[red.idx] <- "red"
  labels.mat[green.idx] <- "green"
  diag(labels.mat) <- "none"
  diag(adjacency.mat) <- 0

  d.red <- x.red*apply(labels.mat,1, count.red)
  d.green <- x.green*apply(labels.mat, 1, count.green)
  deg <- d.red + d.green
  
  vertex.colors <- c(seq("red","red",length.out=m),
                     seq("green","green",length.out=(n-m)))

  return(list(adjacency=adjacency.mat,
              labels=labels.mat,
              dred=d.red,
              dgreen=d.green,
              degree=deg,
              v.colors = vertex.colors))
}

## Compute the scan statistics
## If the scan statistic belong
scan.statistic <- function(x,g){
  x.red <- x[1]
  x.green <- x[2]
  nvertex <- nrow(g$adjacency)
  locality.vec <- seq(0,0,length.out = nvertex)
  for(i in 1:nvertex){
    i.neighbours <- which(g$adjacency[i,] == 1)
    idx <- c(i,i.neighbours)
    s.red <- sum(g$labels[idx,idx] == "red")
    s.green <- sum(g$labels[idx,idx] == "green")
    locality.vec[i] <- s.red*x.red + s.green*x.green
  }
  return(max(locality.vec))
}

## Test the power for the first approximation at theta

power.test2 <- function(theta){
  power <- 0
  mc <-  1000
  
  g0.normalized<- seq(0,0,length.out = mc)
  gA.normalized <- seq(0,0,length.out = mc)
  
  for(i in 1:mc){
    tmp2 <- 0
    gseq <- seq(0,0,length.out = 10)
    for(j in 1:10){
      g <- kidney.egg.attributed(100,9,theta, pi0 = c(0.1,0.3), piA = c(0.1,0.3))
      gseq[j] <- scan.statistic(c(cos(theta),sin(theta)), g)
    }
    gseq.mean <- sum(gseq)/10
    gseq.sse <- sqrt(1/9*sum((gseq - gseq.mean)^2))

    g.new <- kidney.egg.attributed(100,9,theta, pi0 = c(0.1,0.3), piA = c(0.1,0.3))
    g.new.scan <- scan.statistic(c(cos(theta),sin(theta)), g.new)

    if(gseq.mean == 0){
      if(g.new.scan == 0)
        g0.normalized[i] <- 0
      else
        g0.normalized[i] <- Inf
    }
    else{
      g0.normalized[i] <- (g.new.scan - gseq.mean)/gseq.sse
    }
    
    g.new <- kidney.egg.attributed(100,9,theta, pi0 = c(0.1,0.3), piA = c(0.25,0.4))
    g.new.scan <- scan.statistic(c(cos(theta),sin(theta)), g.new)

    if(gseq.mean == 0){
      if(g.new.scan == 0)
        gA.normalized[i] <- 0
      else
        gA.normalized[i] <- Inf
    }
    else{
      gA.normalized[i] <- (g.new.scan - gseq.mean)/gseq.sse
    }
  }

  g0.normalized.sorted <- sort(g0.normalized)
  cv <- g0.normalized.sorted[floor(0.95 * mc)]
  power <- length(which(gA.normalized > cv))/mc

  return(list(cv=cv,power=power))
}

power.test1 <- function(){
  power <- 0
  mc <-  1000
  for(i in 1:mc){
    tmp2 <- 0
    gseq <- seq(0,0,length.out = 10)
    for(j in 1:10){
      g <- kidney.egg.attributed(100,9,theta = 0, pi0 = c(0.1,0), piA = c(0.1,0))
      gseq[j] <- locality(g)
    }
    gseq.mean <- sum(gseq)/10
    gseq.sse <- sqrt(1/9*sum((gseq - gseq.mean)^2))

    g.new <- kidney.egg.attributed(100,9,theta = 0, pi0 = c(0.1,0), piA = c(0.25,0))
    g.new.scan <- locality(g.new)

    calpha <- pt(0.05, 9, ncp = 0, lower.tail = FALSE)

    if(gseq.mean == 0){
      if(g.new.scan > 0)
        power <- power + 1
    }
    else{
      if((g.new.scan - gseq.mean)/gseq.sse > calpha)
        power <- power + 1
    }
  }
  return(power/mc)
}
 
