normal.power.test <- function(mu0,sd0,mu1,sd1,past.len,mc,cv){
  g0.normalized <- seq(0,0,length.out = mc)
  gA.normalized <- seq(0,0,length.out = mc)

  g0mat <- matrix(rnorm(mc*(past.len + 1), mu0, sd0), nrow = mc, ncol = past.len + 1)
  g0seq.mean1 <- rowMeans(g0mat[,1:past.len])
  g0seq.sd1 <- sd(t(g0mat[,1:past.len]))
  g0.normalized <-  (g0mat[,past.len+1] - g0seq.mean1)/g0seq.sd1

  gAvec <- rnorm(mc,mu1,sd1)
  g0seq.mean2 <- rowMeans(g0mat[,2:(past.len+1)])
  g0seq.sd2 <- sd(t(g0mat[,2:(past.len+1)]))
  gA.normalized <- (gAvec - g0seq.mean2)/g0seq.sd2

  power <- length(which(gA.normalized > cv))/mc
  return(list(cv = cv, power = power))
}

## numtri.pwr <- numeric(158)
## numtri.cv <- numeric(158)
## for(i in 1:158){
##   tmp <-  normal.power.test(numtri.null.mean[i],numtri.null.sd[i],numtri.alt.mean[i],
##                             numtri.alt.sd[i], 10,10000,numtri.stat$cv[i])
##   numtri.cv[i] <- tmp$cv
##   numtri.pwr[i] <- tmp$power
## }
  
num.triangles.asy <- function(nvertex, nvertex.abnormal, pi0, piA,past.len){
  
  pi00 <- c(1 - sum((pi0*pi0)[-1]), (pi0*pi0)[-1])
  pi0A <- c(1 - sum((pi0*piA)[-1]), (pi0*piA)[-1])
  piA0 <- pi0A
  piAA <- c(1 - sum((piA*piA)[-1]), (piA*piA)[-1])

  alpha <- 0.05
  x1seq <- seq(0, pi/2, by = 0.01)
  W <- matrix(0, nrow = length(x1seq), ncol = 3)
  cv <- numeric(length(x1seq))
  power <- numeric(length(x1seq))

  for(i in 1:length(x1seq))
    W[i,] <- c(0, cos(x1seq[i]), sin(x1seq[i]))

  for(i in 1:length(x1seq)){
    w.pi00 <- sum(W[i,]*pi00)
    wsq.pi00 <- sum((W[i,])^2*pi00)
    w.rho00 <- wsq.pi00 - (w.pi00)^2
    
    mu0 <- choose(nvertex,3)*(w.pi00)^3
    sigma0 <- (nvertex - 2)*(w.pi00)^2*sqrt(choose(nvertex,2)*w.rho00)

    w.pi0A <- sum(W[i,]*pi0A); w.piA0 <- w.pi0A;
    wsq.pi0A <- sum((W[i,])^2*pi0A); wsq.piA0 <- wsq.pi0A
    w.rho0A <- wsq.pi0A - (w.pi0A)^2; w.rhoA0 <- w.rho0A

    w.piAA <- sum(W[i,]*piAA)
    wsq.piAA <- sum((W[i,])^2*piAA)
    w.rhoAA <- wsq.piAA - (w.piAA)^2

    muA.1 <- choose(nvertex.abnormal,3)*(w.piAA)^3
    muA.2 <- choose(nvertex.abnormal,2)*(nvertex - nvertex.abnormal)*w.piAA*(w.pi0A)^2
    muA.3 <- nvertex.abnormal*choose(nvertex - nvertex.abnormal, 2)*(w.pi0A)^2*w.pi00
    muA.4 <- choose(nvertex - nvertex.abnormal,3)*(w.pi00)^3
    muA <- muA.1 + muA.2 + muA.3 + muA.4

    S1 <- ((nvertex.abnormal - 2)*(w.piAA)^2 + (nvertex - nvertex.abnormal)*(w.pi0A)^2)^2
    S2 <- ((nvertex.abnormal - 1)*w.piAA*w.pi0A + (nvertex - nvertex.abnormal - 1)*w.pi00*w.pi0A)^2
    S3 <- (nvertex.abnormal*(w.pi0A)^2 + (nvertex - nvertex.abnormal - 2)*(w.pi00)^2)^2

    sigmaA.1 <- choose(nvertex.abnormal,2)*w.rhoAA*S1
    sigmaA.2 <- nvertex.abnormal*(nvertex - nvertex.abnormal)*w.rho0A*S2
    sigmaA.3 <- choose(nvertex - nvertex.abnormal,2)*w.rho00*S3
    sigmaA <- sqrt(sigmaA.1 + sigmaA.2 + sigmaA.3)

    
    sigma.ratio <- (sigmaA/sigma0)

    quantile.1 <- qt(alpha, past.len - 1, ncp = 0, lower.tail = FALSE) 
    quantile.2 <- (past.len/past.len+1)*(sigma.ratio^2 - 1) + 1
    quantile <- quantile.1*(1/sqrt(quantile.2)) 

    cv[i] <- quantile

    ncpA <- (muA - mu0)/sqrt(sigmaA^2 + sigma0^2/past.len)

    power[i] <- pt(quantile, past.len - 1, ncp = ncpA, lower.tail = FALSE) 

  }
  return(list(cv = cv, power = power))
}

                              

                              
                    
  
  
