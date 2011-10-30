youngser1 <- function(tstat.mat, past.len, mc){
  gamma <- 0.577215665
  d1 <- dim(tstat.mat)[1]
  d2 <- dim(tstat.mat)[2]
  power.vec <- numeric(d2)
  cv.vec <- numeric(d2)
  for(i in 1:d2){
    a.i <- tstat.mat[,i,]
    a.i.null <- a.i[1:(d1-1),]
    a.i.null.mean <- mean(a.i.null)
    a.i.null.sd <- sd(as.vector(a.i.null))

    a.i.alt <- a.i[d1,]
    a.i.alt.mean <- mean(a.i.alt)
    a.i.alt.sd <- sd(a.i.alt)

    bn0 <- a.i.null.sd*sqrt(6)/pi
    an0 <- a.i.null.mean - bn0*gamma

    bnA <- a.i.alt.sd*sqrt(6)/pi
    anA <- a.i.alt.mean - bnA*gamma

    gpt <- gumbel.power.test(an0, bn0, anA, bnA, past.len, mc)
    power.vec[i] <- gpt$power
    cv.vec[i] <- gpt$cv
  }
  return(list(cv = cv.vec, power = power.vec))
}

pwr1.test <- function(){
  theta.seq <- seq(from = 0, to = pi/4, by = 0.01)
  power.vec <- seq(0,0,length.out = length(theta.seq))
  cv.vec <- seq(0,0,length.out = length(theta.seq))
  for(i in 1:length(theta.seq)){
    print(i)
    tmp <- power.test(theta.seq[i])
    power.vec[i] <- tmp$power
    cv.vec[i] <- tmp$cv
  }
  return(list(cv = cv.vec, power = power.vec))
}

scan.mc.mat <- function(nvertex,nvertex.abnormal,pi0.vec, piA.vec, mc){
  theta.seq <- seq(from = 0, to = pi/4, by = 0.01)
  g0.mat <- matrix(0,nrow=length(theta.seq), ncol = mc)
  gA.mat <- matrix(0,nrow=length(theta.seq), ncol = mc)

  print("Long ass computation: theta = ")
  for(i in 1:length(theta.seq)){
    print(theta.seq[i])
    g0.mat[i,] <- scan.mc(theta.seq[i], nvertex, nvertex.abnormal, pi0.vec, pi0.vec, mc)
    gA.mat[i,] <- scan.mc(theta.seq[i], nvertex, nvertex.abnormal, pi0.vec, piA.vec, mc)
  }
  return(list(g0 = g0.mat, gA = gA.mat))
}
    
## Test the power for the first approximation at theta
power.test <- function(theta,mc,past.len = 10){
  
  g0.normalized<- seq(0,0,length.out = mc)
  gA.normalized <- seq(0,0,length.out = mc)
  
  for(i in 1:mc){
    tmp2 <- 0
    g0seq <- seq(0,0,length.out = past.len + 1)
    for(j in 1:(past.len + 1)){
      g <- kidney.egg.attributed(100,9,theta, pi0 = c(0.1,0.3), piA = c(0.1,0.3))
      g0seq[j] <- scan.statistic(c(cos(theta),sin(theta)), g)
    }
    g0seq.mean1 <- mean(g0seq[1:past.len])
    g0seq.sd1 <- sd(g0seq[1:past.len])

    g0.normalized[i] <- (g0seq[past.len + 1] - g0seq.mean1)/g0seq.sd1
    
    gA.new <- kidney.egg.attributed(100,9,theta, pi0 = c(0.1,0.3), piA = c(0.25,0.4))
    gA.new.scan <- scan.statistic(c(cos(theta),sin(theta)), gA.new)

    g0seq.mean2 <- mean(g0seq[2:(past.len+1)])
    g0seq.sd2 <- sd(g0seq[2:(past.len + 1)])
    gA.normalized[i] <- (gA.new.scan - g0seq.mean2)/g0seq.sd2
  }

  cv <- quantile(g0.normalized, probs = c(0.95))
  power <- (sum(gA.normalized >= cv))/mc

  return(list(cv=cv,power=power))
}

tmp.fun <- function(g0.mat, gA.mat){
  d1 <- nrow(g0.mat)
  power.vec <- numeric(d1)
  cv.vec <- numeric(d1)
  for(i in 1:d1){
    d2 <- floor(length(g0.mat[i,])/11)
    g0.i <- g0.mat[i,1:(11*d2)] 
    g0.i.mat <- matrix(g0.i, ncol = 11)

    g0.i.mat.mean <- rowMeans(g0.i.mat[,1:10])
    g0.i.mat.sd <- sd(t(g0.i.mat[,1:10]))
    g0.i.vec <- (g0.i.mat[,11] - g0.i.mat.mean)/g0.i.mat.sd
    cv.vec[i] <- quantile(g0.i.vec, probs = c(0.95))
    
    power.vec[i] <- sum((gA.mat[i,1:(11*d2)] - g0.i.mat.mean)/g0.i.mat.sd > cv.vec[i])/(11*d2)
  }
  return(list(cv = cv.vec, power = power.vec))
}

gumbel.power.test <- function(loc0, scale0, locA, scaleA, past.len, mc, cv){
  
  g0.normalized<- seq(0,0,length.out = mc)
  gA.normalized<- seq(0,0,length.out = mc)
  #cv <- seq(0,0,length.out = mc)
  #power <- seq(0,0,length.out = mc)
  
  for(i in 1:mc){
    g0seq <- rgumbel(past.len+1,loc0,scale0)
    g0seq.mean1 <- mean(g0seq[1:past.len])
    g0seq.sd1 <- max(sd(g0seq[1:past.len]),1)
    
    g0.normalized[i] <- (g0seq[past.len + 1] - g0seq.mean1)/g0seq.sd1

    ## loc0.adj <- (loc0 - g0seq.mean1)/g0seq.sd1
    ## scale0.adj <- scale0/g0seq.sd1

    ## locA.adj <- (locA - g0seq.mean2)/g0seq.sd2
    ## scaleA.adj <- scaleA/g0seq.sd2

    ## cv <- qgumbel(0.95, loc0.adj, scale0.adj, lower.tail = TRUE)
    ## power <- pgumbel(cv, locA.adj, scaleA.adj, lower.tail = FALSE)
    
    ## g0.new <- rgumbel(1,loc0,scale0)
    ## gA.new <- rgumbel(1,locA,scaleA)
    gA <- rgumbel(1,locA,scaleA)
    
    g0seq.mean2 <- mean(g0seq[2:(past.len+1)])
    g0seq.sd2 <- max(sd(g0seq[2:(past.len+1)]),1)
    gA.normalized[i] <- (gA - g0seq.mean2)/g0seq.sd2
    }

  g0.normalized.sort <- sort(g0.normalized, decreasing = FALSE)
#  cv <- g0.normalized.sort[floor(0.95*mc)]
  power <- length(which(gA.normalized > cv))/mc

  return(list(cv=cv,power=power))
}

cv.vals <- function(nvertex, nvertex.abnormal, pi0, piA, past.len, mc){
  
  theta <- seq(0,pi/2, by = 0.01)
  alpha <- 0.05
  cv <- seq(0,pi/2, by = 0.01)

  for(i in 1:length(theta)){
    scan.i <- scan.mc(theta[i], nvertex, nvertex.abnormal, pi0, piA, (past.len + 1)*mc)

    xyz <- matrix(scan.i, ncol = (past.len + 1))
    xyz1 <- xyz[,1:past.len]
    xyz2 <- xyz[,past.len + 1]

    mean.vec <- rowMeans(xyz)
    sd.vec <- sd(t(xyz))

    norm.seq <- (xyz2 - mean.vec)/sd.vec

    cv[i] <- quantile(norm.seq, probs = c(1 - alpha))
  }
  return(cv)
}

gumbel.kidney.egg <- function(nvertex, nvertex.abnormal, theta, pi0, piA, num.sample, from.egg = TRUE){

  pi00 <- c(1 - sum((pi0*pi0)[2:3]), (pi0*pi0)[2:3])
  pi0A <- c(1 - sum((pi0*piA)[2:3]), (pi0*piA)[2:3])
  piA0 <- pi0A
  piAA <- c(1 - sum((piA*piA)[2:3]), (piA*piA)[2:3])

  w <- c(0,cos(theta), sin(theta))

  loc.scale <- null.alternative.params(nvertex, nvertex.abnormal, pi00, pi0A, piAA, w, from.egg)
  anA <- loc.scale$anA
  bnA <- loc.scale$bnA
  return(rgumbel(num.sample, anA, bnA))
}

null.alternative.params <- function(nvertex, nvertex.abnormal, pi00, pi0A, piAA, w, from.egg = TRUE){
  w.pi00 <- sum(w*pi00)
  w.pi0A <- sum(w*pi0A); w.piA0 <- w.pi0A
  w.piAA <- sum(w*piAA)
  wsq.pi00 <- sum(w^2*pi00)
  wsq.pi0A <- sum(w^2*pi0A); wsq.piA0 <- wsq.pi0A
  wsq.piAA <- sum(w^2*piAA)
  s00 <- sum(pi00[2:3])
  s0A <- sum(pi0A[2:3])
  sA0 <- s0A
  sAA <- sum(piAA[2:3])

  mu0 <- (nvertex - 1)*s00
  sigma0 <- sqrt((nvertex - 1)*s00*(1 - s00))
  zn <- sqrt(2*log(nvertex))*(1 - (log(log(nvertex)) + log(4*pi))/(4*log(nvertex)))
    
  N0 <- mu0 + zn*sigma0
  an01 <- N0*w.pi00/s00
  an02 <- w.pi00*N0*(N0)/2
  an0 <- an01 + an02
  bn0 <- (1 + w.pi00*N0)*sigma0/sqrt(2*log(nvertex))
  
  mu.B <- nvertex.abnormal * s0A
  mu.C <- (nvertex - nvertex.abnormal - 1)*s00
  mu.E <- (nvertex.abnormal - 1)*sAA
  mu.F <- (nvertex - nvertex.abnormal)*sA0
  
  sigma.B <- sqrt((nvertex.abnormal)*s0A*(1 - s0A))
  sigma.C <- sqrt((nvertex - nvertex.abnormal -1)*s00*(1 - s00))
  
  sigma.E <- sqrt((nvertex.abnormal - 1)*sAA*(1 - sAA))
  sigma.F <- sqrt((nvertex - nvertex.abnormal)*sA0*(1 - sA0))

  if(from.egg == TRUE){
    muA<- mu.E + mu.F
    sigmaA <- sqrt(sigma.E^2 + sigma.F^2)
    
    zm <- sqrt(2*log(nvertex.abnormal))*(1 - (log(log(nvertex.abnormal))
                                              + log(4*pi))/(4*log( nvertex.abnormal)))
    Nkappa <- muA + zm*sigmaA
    anA1 <- Nkappa*w.piAA/sAA
    anA2 <- w.pi00*Nkappa*(Nkappa - 1)/2
    anA3 <- (w.piAA - w.pi00)*mu.E*(mu.E - 1)/2
    anA4 <- (w.piA0 - w.pi00)*mu.E*mu.F
  
    anA <- anA1 + anA2 + anA3 + anA4
    bnA <- (1  + w.pi00* Nkappa)*sigmaA/sqrt(2*log(nvertex.abnormal))
  }
  else{
    muA<- mu.B + mu.C
    sigmaA <- sqrt(sigma.B^2 + sigma.C^2)
    
    zm <- sqrt(2*log(nvertex - nvertex.abnormal))*(1 - (log(log(nvertex - nvertex.abnormal))
                                                        + log(4*pi))/(4*log(nvertex - nvertex.abnormal)))
    
    Nkappa <- muA + zm*sigmaA
    anA1 <- Nkappa*w.pi00/s00
    anA2 <- w.pi00*Nkappa*(Nkappa - 1)/2
    anA3 <- (w.piAA - w.pi00)*mu.B*(mu.B - 1)/2
    anA4 <- (w.piA0 - w.pi00)*mu.B*mu.C
    
    anA <- anA1 + anA2 + anA3 + anA4
    bnA <- (1  + w.pi00* Nkappa)*sigmaA/sqrt(2*log(nvertex - nvertex.abnormal))
  }
  return(list(an0 = an0, bn0 = bn0, anA = anA, bnA = bnA))
}

scan.asy.egg <- function(nvertex, nvertex.abnormal,  pi0, piA, past.len, cv.vec, mc){
  
  pi00 <- c(1 - sum((pi0*pi0)[2:3]), (pi0*pi0)[2:3])
  pi0A <- c(1 - sum((pi0*piA)[2:3]), (pi0*piA)[2:3])
  piA0 <- pi0A
  piAA <- c(1 - sum((piA*piA)[2:3]), (piA*piA)[2:3])

  gamma <- 0.577215665
  alpha <- 0.05
  x1seq <- seq(0,pi/2, by = 0.01)
  W <- matrix(0,nrow = length(x1seq), ncol = 3)
  
  for(i in 1:length(x1seq))
    W[i,] <- c(0,cos(x1seq[i]), sin(x1seq[i]))
  
  power.vec <- seq(0,0,length.out = length(x1seq))
  ##   cv.vec <- seq(0,0,length.out = length(x1seq))
 
  for(i in 1:length(x1seq)){
    loc.scale.list <- null.alternative.params(nvertex, nvertex.abnormal, pi00, pi0A, piAA, W[i,], from.egg = TRUE)
    an0 <- loc.scale.list$an0
    anA <- loc.scale.list$anA
    bn0 <- loc.scale.list$bn0
    bnA <- loc.scale.list$bnA

#    cv <- qt(alpha, df = past.len - 1, lower.tail = FALSE)

#    scale <- bnA/(pi/sqrt(6)*bn0)
#    ncp <- (anA - an0 - bn0*gamma)/(pi/sqrt(6)*bn0)

#    power.vec[i] <- pgumbel(cv.vec[i], ncp, scale, lower.tail = FALSE)
    
     tmp <- gumbel.power.test(an0, bn0, anA, bnA, past.len, mc)
     power.vec[i] <- tmp$power
    ## tmp <- gumbel.power.test(an0.vec[i], bn0.vec[i], anA.vec[i], bnA.vec[i], past.len, mc)
    ## power.vec[i] <- tmp$power
     cv.vec[i] <- tmp$cv
  }
  plot(x1seq, power.vec, type = "l", xlab = "theta", ylab = "beta")
  return(list(cv = cv.vec, power = power.vec))
}

mean.sd.vec <- function(){
  theta.seq <- seq(from = 0, to = pi/4, by = 0.01)
  mean0 <- seq(0,0,length.out = length(theta.seq))
  sd0 <- seq(0,0,length.out = length(theta.seq))
  meanA <- seq(0,0,length.out = length(theta.seq))
  sdA <- seq(0,0,length.out = length(theta.seq))

  for(i in 1:length(theta.seq)){
    print(i)
    xyz <- scan.mc(theta.seq[i], 100, 9, c(0.1,0.3), c(0.1,0.3), 1000)
    mean0[i] <- mean(xyz)
    sd0[i] <- sd(xyz)
    xyz <- scan.mc(theta.seq[i], 100, 9, c(0.1,0.3), c(0.25,0.4), 1000)
    meanA[i] <- mean(xyz)
    sdA[i] <- sd(xyz)
  }
  return(list(mean0 = mean0, sd0 = sd0, meanA = meanA, sdA = sdA))
}

pwr.asy.egg <- function(an0, bn0, anA, bnA, past.len, mc, adjustment){
  power.vec <- seq(0,0,length.out = length(an0))
  for(i in 1:length(an0)){
    power.vec[i] <- gumbel.power.test(an0[i],bn0[i],anA[i] + adjustment,bnA[i], past.len, mc)$power
  }
  plot(seq(from = 0, to = pi/4, by = 0.01), power.vec, type = "l", xlab = "theta", ylab = "beta")
  return(power.vec)
}

scan.asy.kidney <- function(nvertex, nvertex.abnormal, pi0, piA, cv.vec, past.len, mc){
  
  pi00 <- c(1 - sum((pi0*pi0)[2:3]), (pi0*pi0)[2:3])
  pi0A <- c(1 - sum((pi0*piA)[2:3]), (pi0*piA)[2:3])
  piA0 <- pi0A
  piAA <- c(1 - sum((piA*piA)[2:3]), (piA*piA)[2:3])

  gamma <- 0.577215665
  alpha <- 0.05
   
  x1seq <- seq(0,pi/2, by = 0.01)
  W <- matrix(0,nrow = length(x1seq), ncol = 3)
  
  for(i in 1:length(x1seq))
    W[i,] <- c(0,cos(x1seq[i]), sin(x1seq[i]))

  
  power.vec <- seq(0,0,length.out = length(x1seq))
  ## cv.vec <- seq(0,0,length.out = length(x1seq))
 
  for(i in 1:length(x1seq)){
    loc.scale.list <- null.alternative.params(nvertex, nvertex.abnormal, pi00, pi0A, piAA, W[i,], from.egg = FALSE)
    an0 <- loc.scale.list$an0
    anA <- loc.scale.list$anA
    bn0 <- loc.scale.list$bn0
    bnA <- loc.scale.list$bnA
  
    tmp <- gumbel.power.test(an0, bn0, anA, bnA, past.len, mc)
    power.vec[i] <- tmp$power
    ## tmp <- gumbel.power.test(an0.vec[i], bn0.vec[i], anA.vec[i], bnA.vec[i], past.len, mc)
    ## power.vec[i] <- tmp$power
    cv.vec[i] <- tmp$cv
    ## scale <- bnA/(pi/sqrt(6)*bn0)
    ## ncp <- (anA - an0 - bn0*gamma)/(pi/sqrt(6)*bn0)

    ## power.vec[i] <- pgumbel(cv, ncp, scale, lower.tail = FALSE)
    
    ## tmp <- gumbel.power.test(an0.vec[i], bn0.vec[i], anA.vec[i], bnA.vec[i], past.len, mc)
    ## power.vec[i] <- tmp$power
    ## cv.vec[i] <- tmp$cv
  }
  plot(x1seq, power.vec, type = "l", xlab = "theta", ylab = "beta")
  return(list(cv = cv.vec, power = power.vec))
}

weights.seq.scan <- function(nvertex = 100, nvertex.abnormal = 9,
                             pi0 = c(0.6,0.1,0.3), piA = c(0.35,0.25,0.4),
                             num.sample){

  pi00 <- c(1 - sum((pi0*pi0)[2:3]), (pi0*pi0)[2:3])
  pi0A <- c(1 - sum((pi0*piA)[2:3]), (pi0*piA)[2:3])
  piA0 <- pi0A
  piAA <- c(1 - sum((piA*piA)[2:3]), (piA*piA)[2:3])

  alpha <- 0.05
  x1seq <- seq(0,pi/2, by = 0.01)
  W <- matrix(0,nrow = length(x1seq), ncol = 3)
  
  for(i in 1:length(x1seq))
    W[i,] <- c(0,cos(x1seq[i]), sin(x1seq[i]))

  wseq <- seq(0,0,length.out = length(x1seq))
  
  ## power.vec <- seq(0,0,length.out = length(x1seq))
  ## cv.vec <- seq(0,0,length.out = length(x1seq))
 
  for(i in 1:length(x1seq)){

    loc.scale.egg <- null.alternative.params(nvertex, nvertex.abnormal, pi00, pi0A, piAA, W[i,], from.egg = TRUE)
    anA.egg <- loc.scale.egg$anA
    bnA.egg <- loc.scale.egg$bnA

    loc.scale.kidney <- null.alternative.params(nvertex, nvertex.abnormal, pi00, pi0A, piAA, W[i,], from.egg = FALSE)
    anA.kidney <- loc.scale.kidney$anA
    bnA.kidney <- loc.scale.kidney$bnA
    
    g1 <- rgumbel(num.sample, anA.kidney, bnA.kidney)
    g2 <- rgumbel(num.sample, anA.egg, bnA.egg)

    wseq[i] <- sum(g1 >= g2)/num.sample
  }
  return(wseq)
}



cv.vals2 <- function(n = 100, m = 9,
                     pi0 = c(0.6,0.1,0.3), piA = c(0.35,0.25,0.4),
                     past.len = 10,
                     num.sample = 1000){

  pi00 <- c(1 - sum((pi0*pi0)[2:3]), (pi0*pi0)[2:3])
  pi0A <- c(1 - sum((pi0*piA)[2:3]), (pi0*piA)[2:3])
  piA0 <- pi0A
  piAA <- c(1 - sum((piA*piA)[2:3]), (piA*piA)[2:3])

  alpha <- 0.05
  x1seq <- seq(0,pi/2, by = 0.01)
  W <- matrix(0,nrow = length(x1seq), ncol = 3)
  
  for(i in 1:length(x1seq))
    W[i,] <- c(0,cos(x1seq[i]), sin(x1seq[i]))

  cv.vec <- seq(0,0,length.out = length(x1seq))
 
  for(i in 1:length(x1seq)){
    loc.scale.0 <- null.alternative.params(nvertex, nvertex.abnormal, pi00, pi0A, piAA, W[i,])
    an0 <- loc.scale.0$an0
    bn0 <- loc.scale.0$bn0
    cv.vec[i] <- gumbel.cv(an0, bn0, past.len, num.sample)
  }
  return(cv.vec)
}

gumbel.cv <- function(loc, scale, past.len, mc){
  alpha <-  0.05
  
  xyz.sample <- rgumbel((past.len + 1) * mc, loc, scale)
  xyz <- matrix(xyz.sample, ncol = past.len + 1)
  xyz1 <- xyz[,1:past.len]
  xyz2 <- xyz[,past.len + 1]
  
  mean.vec <- rowMeans(xyz)
  sd.vec <- sd(t(xyz))

  cv <- quantile((xyz2 - mean.vec)/sd.vec, probs = c(0.95))
  return(cv)
}
