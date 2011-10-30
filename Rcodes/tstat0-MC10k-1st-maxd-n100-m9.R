read.maxd <- function(){

  load("tstat0-MC10k-1st-maxd-n100-m9.RData")

  maxd0.mat <- matrix(0, ncol = 11*length(tstat0), nrow = length(tstat0[[1]]))

  for(i in 1:nrow(maxd0.mat)){
    for(j in 1:length(tstat0)){
      maxd0.mat[i,((j-1)*11+1):(11*j)] <- tstat0[[j]][[i]][1:11]
    }
  }

  maxdA.mat <- matrix(0, ncol = length(tstat0), nrow = length(tstat0[[1]]))
  for(i in 1:nrow(maxdA.mat)){
    for(j in 1:length(tstat0)){
      maxdA.mat[i,j] <- tstat0[[j]][[i]][12]
    }
  }
}

read.scan <- function(){
  
  load("tstat1-MC10k-1st-scan-n100-m9.Rbin")
  scan.null.mat <- matrix(0, ncol = 11*length(tstat1), nrow = length(tstat1[[1]]))

  for(i in 1:nrow(scan.null.mat)){
    for(j in 1:(ncol(scan.null.mat)/11)){
      scan.null.mat[i,((j-1)*11+1):(j*11)] <- tstat1[[j]][[i]][1:11]
    }
  }

  scan.alt.mat <- matrix(0, ncol = length(tstat1), nrow = length(tstat1[[1]]))
  for(i in 1:nrow(scan.alt.mat)){
    for(j in 1:ncol(scan.alt.mat)){
      scan.alt.mat[i,j] <- tstat1[[j]][[i]][12]
    }
  }
}


read.numtri <- function(){
  tstat3 <- dget("tstat3-MC10k-1st-ntri-n100-m9.RData.gz")
  numtri.null.mat <- matrix(0, ncol = 11*length(tstat2), nrow = length(tstat2[[1]]))

  for(i in 1:nrow(numtri.null.mat)){
    for(j in 1:(ncol(numtri.null.mat)/11)){
      numtri.null.mat[i,((j-1)*11+1):(j*11)] <- tstat2[[j]][[i]][1:11]
    }
  }

  numtri.alt.mat <- matrix(0, ncol = length(tstat2), nrow = length(tstat2[[1]]))
  for(i in 1:nrow(numtri.alt.mat)){
    for(j in 1:ncol(numtri.alt.mat)){
      numtri.alt.mat[i,j] <- tstat2[[j]][[i]][12]
    }
  }
}

compute.pwr <- function(null.mat,alt.mat){

  tstat.null.mat <- matrix(0, nrow = nrow(null.mat), ncol = ncol(null.mat)/11)
  tstat.alt.mat <- matrix(0, nrow = nrow(alt.mat), ncol = ncol(alt.mat))

  for(i in 1:nrow(tstat.null.mat)){
    for(j in 1:ncol(tstat.null.mat)){
      tmp1 <- null.mat[i,((j-1)*11+1):(j*11-1)]
      tmp2 <- mean(tmp1)
      tmp3 <- max(sd(tmp1),1)
      tstat.null.mat[i,j] <- (null.mat[i,11*j] - tmp2)/tmp3
    }
  }

  cv.vec <- apply(tstat.null.mat,1,quantile95)
  
  for(i in 1:nrow(tstat.alt.mat)){
    for(j in 1:ncol(tstat.alt.mat)){
      tmp1 <- null.mat[i,((j-1)*11+2):(j*11)]
      tmp2 <- mean(tmp1)
      tmp3 <- max(sd(tmp1),1)
      tstat.alt.mat[i,j] <- (alt.mat[i,j] - tmp2)/tmp3
    }
  }

  pwr.vec <- numeric(nrow(tstat.alt.mat))

  for(i in 1:nrow(tstat.alt.mat)){
    pwr.vec[i] <- sum(tstat.alt.mat[i,] > cv.vec[i],na.rm = TRUE)/ncol(tstat.alt.mat)
  }
  
    # tstat.mat

  return(list(cv = cv.vec, pwr = pwr.vec))
}

gumbel.scan <- function(){
 scan.pwr <- numeric(158)
 scan.cv <- numeric(158)
 for(i in 1:158){
   tmp <-  gumbel.power.test(loc0.scan[i],scale0.scan[i],loc1.scan[i],
                               scale1.scan[i], 10,10000,scan.stat$cv[i])
   scan.cv[i] <- tmp$cv
   scan.pwr[i] <- tmp$power
   }
}
  gumbel.maxd <- function(){
 scan.pwr <- numeric(158)
 scan.cv <- numeric(158)
 for(i in 1:158){
   tmp <-  gumbel.power.test(loc0.maxd[i],scale0.maxd[i],loc1.maxd[i],
                               scale1.maxd[i], 10,10000,scan.stat$cv[i])
   scan.cv[i] <- tmp$cv
   scan.pwr[i] <- tmp$power
   }
}
  
quantile95 <- function(x){
  quantile(x, probs = c(0.95),na.rm = TRUE)
}

youngser.pwr <- function(){
 nummc <- length(tstat1)
 numseq <- length(tstat1[[1]])
 numt <- length(tstat1[[1]][[1]])
 tstat2 <- array(0,dim=c(nummc,numseq,numt))
 for (i in 1:nummc) {
   for (j in 1:numseq) {
     for (k in 1:numt) {
       tstat2[i,j,k] <- tstat1[[i]][[j]][k]
     }
   }
 }
 rm(tstat1)

 ## calc power
 xseq <- seq(0,pi/2,by=0.01)
 pwr <- rep(0,numseq)
#  norm <- T
 for (i in 1:numseq) {
   tstat3 <- tstat2[,i,] # 10000 x 12

   ## normalize it
   if (norm==T) {
     mu3 <- matrix(0,ncol=ncol(tstat3),nrow=nrow(tstat3))
     sd3 <- matrix(1,ncol=ncol(tstat3),nrow=nrow(tstat3))
     for (tt in theta:ncol(tstat3)) {
       mu3[,tt] <- apply(tstat3[,(tt-tau):(tt-1)],1,mean)
       sd3[,tt] <- apply(tstat3[,(tt-tau):(tt-1)],1,sd)
     }
     ##        tstat3n <- (tstat3-mu3)/sd3
     tstat3n <- (tstat3-mu3)/pmax(sd3,1)
   }
   else {
     tstat3n <- tstat3
   }

   ## calc powers
   T03 <- sort(tstat3n[,theta],dec=T)
   TA3 <- sort(tstat3n[,(theta+1)],dec=T)
   cv3 <- T03[nrow(tstat3n)*alpha+1]
   sz3 <- sum(T03>cv3)/nrow(tstat3n)
   pwr[i] <- sum(TA3>cv3)/nrow(tstat3n)
   cat("angle = ", xseq[i],", cv =", cv3, ", beta =", pwr[i],"\n")
 }
  
  




      
    
