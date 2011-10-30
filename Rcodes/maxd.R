max.asy <- function(nvertex, nvertex.abnormal, pi0, piA, past.len, mc,
                        from.egg = TRUE, adjustment = 0){
  
  pi00 <- c(1 - sum((pi0*pi0)[2:3]), (pi0*pi0)[2:3])
  pi0A <- c(1 - sum((pi0*piA)[2:3]), (pi0*piA)[2:3])
  piA0 <- pi0A
  piAA <- c(1 - sum((piA*piA)[2:3]), (piA*piA)[2:3])

  gamma <- 0.577215665
  alpha <- 0.05
  x1seq <- seq(0,pi/2, by = 0.01)
  W <- matrix(0,nrow = length(x1seq), ncol = 3)
  cv <- numeric(length(x1seq))
  power <- numeric(length(x1seq))
  
  for(i in 1:length(x1seq))
    W[i,] <- c(0,cos(x1seq[i]), sin(x1seq[i]))

  for(i in 1:length(x1seq)){
    params.list <- maxd.0A.params(nvertex, nvertex.abnormal, pi00, pi0A,
                                  piAA, W[i,], from.egg)

    an0 <- params.list$an0
    anA <- params.list$anA
    bn0 <- params.list$bn0
    bnA <- params.list$bnA

    ## sigma0 <- pi/sqrt(6)*bn0
    ## sigmaA <- pi/sqrt(6)*bnA

    ## ncp <- (anA - an0 + (bnA - bn0)*gamma)/sqrt(sigma0^2/past.len + sigmaA^2)
    ## scale <- bnA/(pi/sqrt(6)*bn0)

    ## quartile.1<- qt(alpha, df = past.len - 1, ncp = 0, lower.tail = FALSE)
    ## quartile.2 <- (1 + (past.len)/(past.len + 1)*(sigmaA^2/sigma0^2 - 1))
    ## quartile <- quartile.1/sqrt(quartile.2) + adjustment

    tmp <- gumbel.power.test(an0, bn0, anA, bnA, past.len, mc)
    cv[i] <- tmp$cv
    #power[i] <- pt(quartile, df = past.len - 1, ncp = ncp, lower.tail = FALSE)
    power[i] <- tmp$power
  }
#  plot(x1seq, power, type = "l", xlab = "theta", ylab = "beta")
  return(list(cv = cv, power = power))
}

maxd.power <- function(nvertex, nvertex.abnormal, pi0, piA, wseq,
                       past.len, mc, adjustment){

  maxd.power.egg <- max.asy(nvertex, nvertex.abnormal, pi0, piA, past.len, mc,
                            from.egg = TRUE, adjustment)
  maxd.power.kidney <- max.asy(nvertex, nvertex.abnormal, pi0, piA, past.len, mc,
                               from.egg = FALSE, adjustment = 0)
  power <- wseq*maxd.power.egg$power + (1 - wseq)*maxd.power.kidney$power
  plot(seq(from = 0, to = pi/2, by = 0.01), power, type = "l")
  return(power)
}

maxd.0A.params <- function(nvertex, nvertex.abnormal, pi00, pi0A, piAA, w, from.egg = TRUE){

  w.pi00 <- sum(w*pi00)
  w.pi0A <- sum(w*pi0A); w.piA0 <- w.pi0A
  w.piAA <- sum(w*piAA)
  wsq.pi00 <- sum(w^2*pi00)
  wsq.pi0A <- sum(w^2*pi0A); wsq.piA0 <- wsq.pi0A
  wsq.piAA <- sum(w^2*piAA)

  w.rho00 <- wsq.pi00 - w.pi00^2
  w.rho0A <- wsq.pi0A - w.pi0A^2 ; w.rhoA0 <- w.rho0A
  w.rhoAA <- wsq.piAA - w.piAA^2

  zn0 <- sqrt(2*log(nvertex))*(1 - (log(log(nvertex)) + log(4*pi))/(4*log(nvertex)))
  mu0 <- (nvertex - 1)*w.pi00
  sigma0 <- sqrt((nvertex-1)*w.rho00)
  an0 <- mu0 + zn0*sigma0
  bn0 <- 1/sqrt(2*log(nvertex))*sigma0

  muB <- nvertex.abnormal*w.pi0A
  muC <- (nvertex - nvertex.abnormal - 1)*w.pi00
  muE <- (nvertex.abnormal - 1)*w.piAA
  muF <- (nvertex - nvertex.abnormal)*w.piA0

  sigmaB <- sqrt(nvertex.abnormal*w.rho0A)
  sigmaC <- sqrt((nvertex - nvertex.abnormal - 1)*w.rho00)
  sigmaE <- sqrt((nvertex.abnormal-1)*w.rhoAA)
  sigmaF <- sqrt((nvertex - nvertex.abnormal)*w.rhoA0)
  
  if(from.egg == TRUE){
    muA <- muE + muF
    sigmaA <- sqrt(sigmaE^2 + sigmaF^2)
    znA <- sqrt(2*log(nvertex.abnormal))*(1 - (log(log(nvertex.abnormal)) +
                                               log(4*pi))/(4*log(nvertex.abnormal)))

    anA <- muA + znA*sigmaA
    bnA <- 1/sqrt(2*log(nvertex.abnormal))*sigmaA
  } else{
    muA <- muB + muC
    sigmaA <- sqrt(sigmaB^2 + sigmaC^2)
    l <- nvertex - nvertex.abnormal
    znA <- sqrt(2*log(l))* (1 - (log(log(l)) + log(4*pi))/(4*log(l)))

    anA <- muA + znA*sigmaA
    bnA <- 1/sqrt(2*log(nvertex - nvertex.abnormal))*sigmaA
  }

  return(list(an0 = an0, bn0 = bn0, anA = anA, bnA = bnA))
}

## maxasy1 <- function(a=c(0.6,0.1,0.3),b=c(0.35,0.25,0.4))
## {
## # oripi <- pi
## oripi <- 3.141592654
## #  a = c(0.6, 0.1, 0.3)
## #  b = c(0.35, 0.25, 0.40)

##  pi = c(1-sum((a * a)[2:3]), (a * a)[2:3])
##  phi = c(1-sum((a * b)[2:3]), (a * b)[2:3])
##  nu = c(1-sum((b * b)[2:3]), (b * b)[2:3])

##  W = rbind(
##    c(0,1,1),
##    c(0,1,0),
##    c(0,2,1)
##    )

##  x1seq = seq(0,2*3.14,by=0.01)
##  W=matrix(0,nrow=length(x1seq),ncol=3)

##  for(i in 1:length(x1seq))
##    W[i,] = c(0,cos(x1seq[i]),sin(x1seq[i]))


##  for( scale in 1:1 ) {

##    n =10
##    nvertex = 100*scale
##    abnorm = 9
##    chart = numeric(abnorm)

##    my = matrix(nrow=abnorm,ncol=length(x1seq))

##    for (row in 1:length(x1seq)) {
##      for ( i in abnorm:abnorm) {
##        alpha = 0.05
##        m = i

##        xi_pi = t(W[row,])%*% out(pi) %*% W[row,]
##        xi_phi = t(W[row,])%*% out(phi) %*% W[row,]
##        xi_nu = t(W[row,])%*% out(nu) %*% W[row,]

##        An = sqrt(2*log(nvertex - m))*(1 - log(log(nvertex - m))/(4*log(nvertex - m)) - log(2*sqrt(oripi))/(2*log(nvertex - m)))
##        Bn = 1/sqrt(2*log(nvertex -m))
##        em.constant = 0.57721

##        quartile_1 = qt(alpha,n-1,ncp = 0, lower.tail = FALSE)
##        quartile_2 = 1 + (n/(n+1))*((m)/(nvertex-1)*(xi_phi/xi_pi - 1) +
##          (nvertex - m - 1)/(nvertex-1)*(xi_pi/xi_pi - 1))
##        quartile = quartile_1/sqrt(quartile_2)

## #       kappa_num = (An + Bn*em.constant)*(sqrt((m)*xi_phi + (nvertex - m - 1)*xi_pi) - sqrt((nvertex - 1)*xi_pi)) +
## #        (m) * W[row,] %*% (phi-pi) + (nvertex - m - 1 )* W[row,] %*% (pi-pi)
##         kappa_num =  (m) * W[row,] %*% (phi-pi) + (nvertex - m - 1 )* W[row,] %*% (pi-pi)
##        kappa_den_1 = (1+1/n)*(nvertex - 1)* xi_pi
##        kappa_den_2 = (m)* ( xi_phi - xi_pi )
##        kappa_den_3 = (nvertex - m - 1)*( xi_pi - xi_pi)
## #       kappa = kappa_num/(oripi/sqrt(6)*Bn*sqrt(kappa_den_1+kappa_den_2 + kappa_den_3))
##        kappa = kappa_num/(oripi/sqrt(6)*Bn*sqrt(kappa_den_1+kappa_den_2 + kappa_den_3))
## #	kappa = 0
##        chart[i]=(pt(quartile, n-1, ncp = kappa , lower.tail = FALSE))
##    #    print(c(scale,i,kappa_num^2,kappa_den_1,kappa_den_2,kappa_den_3,kappa, An, Bn))
##      }

##      my[,row] = chart[1:abnorm]
##    }
##  }

## #  pdf("namasy.pdf")
##  plot(x1seq[1:length(x1seq)/4],my[9,1:length(x1seq)/4],type="l",xlab="theta",ylab="beta")
##  abline(h=0.05)
##  abline(v=0.465)
##  abline(v=0)
##  abline(v=0.785)
## #  dev.off()
##  return(my)
## }

## maxasy2 <- function(a=c(0.6,0.1,0.3),b=c(0.35,0.25,0.4))
## {
## # oripi <- pi
## oripi <- 3.141592654
## #  a = c(0.6, 0.1, 0.3)
## #  b = c(0.35, 0.25, 0.40)

##  pi = c(1-sum((a * a)[2:3]), (a * a)[2:3])
##  phi = c(1-sum((a * b)[2:3]), (a * b)[2:3])
##  nu = c(1-sum((b * b)[2:3]), (b * b)[2:3])

##  W = rbind(
##    c(0,1,1),
##    c(0,1,0),
##    c(0,2,1)
##    )

##  x1seq = seq(0,2*3.14,by=0.01)
##  W=matrix(0,nrow=length(x1seq),ncol=3)

##  for(i in 1:length(x1seq))
##    W[i,] = c(0,cos(x1seq[i]),sin(x1seq[i]))

##  for( scale in 1:1 ) {

##    n =10
##    nvertex = 100*scale
##    abnorm = 9
##    chart = numeric(abnorm)

##    my = matrix(nrow=abnorm,ncol=length(x1seq))

##    for (row in 1:length(x1seq)) {
##      for ( i in abnorm:abnorm) {
##        alpha = 0.05
##        m = i

##        xi_pi = t(W[row,])%*% out(pi) %*% W[row,]
##        xi_phi = t(W[row,])%*% out(phi) %*% W[row,]
##        xi_nu = t(W[row,])%*% out(nu) %*% W[row,]

##        An = sqrt(2*log(m))*(1 - log(log(m))/(4*log(m)) - log(2*sqrt(oripi))/(2*log(m)))
##        Bn = 1/sqrt(2*log(m))
##        em.constant = 0.57721

##        quartile_1 = qt(alpha,n-1,ncp = 0, lower.tail = FALSE)
##        quartile_2 = 1 + (n/(n+1))*((m-1)/(nvertex-1)*(xi_nu/xi_pi - 1) +
##          (nvertex - m)/(nvertex-1)*(xi_phi/xi_pi - 1))
##        quartile = quartile_1/sqrt(quartile_2) + Bn

## #       kappa_num = (An + Bn*em.constant)*(sqrt((m-1)*xi_nu + (nvertex - m)*xi_phi) - sqrt((nvertex - 1)*xi_pi)) +
## #        (m - 1) * W[row,] %*% (nu-pi) + (nvertex - m )* W[row,] %*% (phi-pi)
##        kappa_num =  (m - 1) * W[row,] %*% (nu-pi) + (nvertex - m )* W[row,] %*% (phi-pi)
##        kappa_den_1 = (1+1/n)*(nvertex - 1)* xi_pi
##        kappa_den_2 = (m - 1)* ( xi_nu - xi_pi )
##        kappa_den_3 = (nvertex -m)*( xi_phi - xi_pi)
##        kappa = kappa_num/(oripi/sqrt(6)*Bn*sqrt(kappa_den_1+kappa_den_2 + kappa_den_3))
##        #kappa = kappa_num/(oripi/sqrt(6)*Bn*sqrt(kappa_den_1+kappa_den_2 + kappa_den_3))

##        chart[i]=(pt(quartile, n-1, ncp = kappa , lower.tail = FALSE))
##    #    print(c(scale,i,kappa_num^2,kappa_den_1,kappa_den_2,kappa_den_3,kappa, An, Bn))
##      }

##      my[,row] = chart[1:abnorm]
##    }
##  }

## #  pdf("namasy.pdf")
##  plot(x1seq[1:length(x1seq)/4],my[9,1:length(x1seq)/4],type="l",xlab="theta",ylab="beta")
##  abline(h=0.05)
##  abline(v=0.465)
##  abline(v=0)
##  abline(v=0.785)
## #  dev.off()
##  return(my)
## }

## out <-function(x) {
##   rbind(
##         c(x[1]*(1-x[1]),-x[1]*x[2],-x[1]*x[3]),
##         c(-x[1]*x[2],x[2]*(1-x[2]),-x[2]*x[3]),
##         c(-x[1]*x[3],-x[2]*x[3],x[3]*(1-x[3]))
##         )
## }

## w.mc.multinom <- function(n,m,x,pi,phi,nu,mc){

##   xyz <- 0
  
##   for(i in 1:mc){
##     x1 <- t(x)%*%(rmultinom(n-m,size=m,phi) + rmultinom(n-m,size=n-m-1,pi))
##     x2 <- t(x)%*%(rmultinom(m,size=m-1,nu) + rmultinom(m,size=n-m,phi))

##     x1.max <- max(x1)
##     x1.max.num <- length(which(x1 == x1.max))
##     x2.max <- max(x2)
##     x2.max.num <- length(which(x2 == x2.max))

##     if(x1.max > x2.max)
##       xyz <- xyz + 1
##     if(x1.max == x2.max)
##       xyz <- xyz + x1.max.num/(x1.max.num + x2.max.num)
##   }

##   return(xyz/mc)
## }

## w.scan.mvn <- function(n,m,x,a,b,mc){
  
##   pi.00 <- c(1 - sum((a*a)[2:3]), (a*a)[2:3])
##   pi.01 <- c(1 - sum((a*b)[2:3]), (a*b)[2:3])
##   pi.10 <- pi.01
##   pi.11 <- c(1 - sum((b*b)[2:3]), (b*b)[2:3])

##   w <- c(0,cos(x), sin(x))

##   x.pi00 <- sum(w*pi.00)
##   x.pi01 <- sum(w*pi.01)
##   x.pi10 <- x.pi01
##   x.pi11 <- sum(w*pi.11)
##   xsq.pi00 <- sum(w^2 * pi.00)
##   xsq.pi01 <- sum(w^2 * pi.01)
##   xsq.pi10 <- xsq.pi01
##   xsq.pi11 <- sum(w^2 * pi.11)

##   p00 <- (x.pi00)^2/(xsq.pi00)
##   p01 <- (x.pi01)^2/(xsq.pi01)
##   p10 <- p01
##   p11 <- (x.pi11)^2/(xsq.pi11)

##   mu1 <- m*p10 + (n-m-1)*p00
##   sigma1 <- sqrt(m*p10*(1-p10) + (n-m-1)*p00*(1-p00))

##   mu2 <- (m-1)*p11 + (n - m)*p10
##   sigma2 <- sqrt((m-1)*p11*(1-p11) + (n - m)*p10*(1-p10))
  
##   xyz <- 0
  
##   for(i in 1:mc){
##     x1 <- rnorm(n-m,mu1,sigma1)
##     x2 <- rnorm(m,mu2,sigma2)
	
##     x1 <- round(x1 + 0.5)
##     x2 <- round(x2 + 0.5)

##     x1.max <- max(x1)
##     x1.max.num <- length(which(x1 == x1.max))
##     x2.max <- max(x2)
##     x2.max.num <- length(which(x2 == x2.max))

##     if(x1.max > x2.max)
##       xyz <- xyz + 1
##     if(x1.max == x2.max)
##       xyz <- xyz + x1.max.num/(x1.max.num + x2.max.num)
##   }
##   return(xyz/mc)
## }

## weights.seq.scan <- function(a = c(0.6,0.1,0.3), b = c(0.35,0.25,0.4), n = 100, m = 9, mc = 10000){

##   n <- 100
##   m <- 9

##   x1seq <- seq(0,2*pi,by=0.01)
##   W <- matrix(0,nrow=length(x1seq),ncol=3)

##   for(i in 1:length(x1seq))
##     W[i,] <- c(0,cos(x1seq[i]),sin(x1seq[i]))

##   wseq <- seq(0,0,length.out = length(x1seq))
  
##   for(i in 1:length(x1seq)){
##       wseq[i] <- w.scan.mvn(n,m,x1seq[i],a,b,mc)
##   }
##   return(wseq)
## }

## w.mc.mvn <- function(n,m,x,pi,phi,nu,mc){

##   xi_pi <-  t(x) %*% out(pi) %*% x
##   xi_phi <- t(x) %*% out(phi) %*% x
##   xi_nu <-  t(x) %*% out(nu) %*% x
  
##   mu1 <- (m)*t(x)%*%phi + (n-m-1)*t(x)%*%pi
##   sigma1 <- sqrt((m)*xi_phi + (n-m-1)*xi_pi)

##   mu2 <- (m-1)*t(x)%*%nu + (n - m)*t(x)%*%phi
##   sigma2 <- sqrt((m-1)*xi_nu + (n - m)*xi_phi)

##   xyz <- 0
  
##   for(i in 1:mc){
##     x1 <- rnorm(n-m,mu1,sigma1)
##     x2 <- rnorm(m,mu2,sigma2)
	
##     x1 <- round(x1 + 0.5)
##     x2 <- round(x2 + 0.5)

##     x1.max <- max(x1)
##     x1.max.num <- length(which(x1 == x1.max))
##     x2.max <- max(x2)
##     x2.max.num <- length(which(x2 == x2.max))

##     if(x1.max > x2.max)
##       xyz <- xyz + 1
##     if(x1.max == x2.max)
##       xyz <- xyz + x1.max.num/(x1.max.num + x2.max.num)
##   }

##   return(xyz/mc)
## }

## weights.seq <- function(type="mvn"){

##   a <- c(0.6,0.1,0.3)
##   b <- c(0.35,0.25,0.4)
##   n <- 100
##   m <- 9
##   mc <- 10000

##   oripi <- 3.141592654
##   pi <- c(1-sum((a * a)[2:3]), (a * a)[2:3])
##   phi <- c(1-sum((a * b)[2:3]), (a * b)[2:3])
##   nu <- c(1-sum((b * b)[2:3]), (b * b)[2:3])
  
##   x1seq <- seq(0,2*oripi,by=0.01)
##   W <- matrix(0,nrow=length(x1seq),ncol=3)

##   for(i in 1:length(x1seq))
##     W[i,] <- c(0,cos(x1seq[i]),sin(x1seq[i]))

##   wseq <- seq(0,0,length.out = length(x1seq))
  
##   for(i in 1:length(x1seq)){
##     if(type == "mvn")
##       wseq[i] <- w.mc.mvn(n,m,W[i,],pi,phi,nu,mc)
##     if(type == "multinom")
##       wseq[i] <- w.mc.multinom(n,m,W[i,],pi,phi,nu,mc)
##   }

##   return(wseq)
## }

## power.fun <- function(){
##   abc1 <- weights.seq(type="mvn")
##   abc2 <- weights.seq(type="multinom")
##   pow1 <- maxasy1()
##   pow2 <- maxasy2()
##   powabc1 <- abc1*pow1 + (1 - abc1)*pow2
##   powabc2 <- abc2*pow1 + (1 - abc2)*pow2
##   return(list(pow1 = powabc1, pow2 = powabc2))
## }

## w.mc <- function(mc,n=100,m=9,x=c(1,0),pi0 = c(0.1,0.3),piA = c(0.25,0.4)){

##   total <- 0
##   g.list <- list()
##   for(i in 1:mc){
##     g <- kidney.egg.attributed(n,m,x,pi0,piA)
##     l <- which(g$degree == max(g$degree))
##     u <- which(l <= m)
##     total <- total + length(u)/length(l)
##     g.list[[i]] <- g
##   }
##   return(list(p=total/mc,glist = g.list))
## }

## abc <- function(mc,n,p){

##   x <- seq(0,0,length.out = mc)
##   y <- rbinom(mc,n-1,p)
##   for(i in 1:mc){
##     x[i] <- rbinom(1,(y[i])*(y[i]-1)/2,p)
##   }
##   hist(x)
##   return(x)
## }
  
## locality <- function(g){
##   nvertex <- nrow(g$adjacency)
##   locality.vec <- seq(0,0,length.out = nvertex)
##   for(i in 1:nvertex){
##     i.neighbours <- which(g$adjacency[i,] == 1)
##     A.tmp <- g$adjacency[c(i,i.neighbours),c(i,i.neighbours)]
##     locality.vec[i] <- sum(A.tmp)
##   }
##   tmp1 <- which(locality.vec == max(locality.vec))
##   tmp2 <- which(g$v.colors[tmp1] == "red")
##   return(length(tmp2)/length(tmp1))
## }
## tmp1 <- 0
## for(i in 1:10000){
##   g <- kidney.egg.attributed(100,9,c(1,0),c(0.1,0),c(0.25,0))
##   tmp1 <- tmp1 + locality(g)
## }


  
       
