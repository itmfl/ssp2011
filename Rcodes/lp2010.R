#library("igraph")

namasy <- function(a=c(0.6,0.1,0.3),b=c(0.35,0.25,0.4))
{
# oripi <- pi
#  a = c(0.6, 0.1, 0.3)
#  b = c(0.35, 0.25, 0.40)

 pi = c(1-sum((a * a)[2:3]), (a * a)[2:3])
 phi = c(1-sum((a * b)[2:3]), (a * b)[2:3])
 nu = c(1-sum((b * b)[2:3]), (b * b)[2:3])

 W = rbind(
   c(0,1,1),
   c(0,1,0),
   c(0,2,1)
   )

 x1seq = seq(0,2*3.14,by=0.01)
 W=matrix(0,nrow=length(x1seq),ncol=3)

 for(i in 1:length(x1seq))
   W[i,] = c(0,cos(x1seq[i]),sin(x1seq[i]))


 for( scale in 1:1 ) {

   n =10
   nvertex = 100*scale
   abnorm = 9#*scale
   chart = numeric(abnorm)

   my = matrix(nrow=abnorm,ncol=length(x1seq))

   for (row in 1:length(x1seq)) {
     for ( i in abnorm:abnorm) {
       alpha = 0.05
       m = i

       xi_pi = t(W[row,])%*% out(pi) %*% W[row,]
       xi_phi = t(W[row,])%*% out(phi) %*% W[row,]
       xi_nu = t(W[row,])%*% out(nu) %*% W[row,]

       quartile_1 = qt(alpha,n-1,ncp = 0, lower.tail = FALSE)
       quartile_2 = 1 + (n/(n+1))*(choose(m,2)/choose(nvertex,2)*(xi_nu/xi_pi - 1) +
         m*(nvertex - m)/choose(nvertex,2)*(xi_phi/xi_pi - 1))
       quartile = quartile_1 / sqrt(quartile_2)

       kappa_num = choose(m,2) * W[row,] %*% (nu-pi) + m * (nvertex -m) * W[row,] %*% (phi-pi)
       kappa_den_1 = (1+1/n)*choose(nvertex,2)* xi_pi
       kappa_den_2 = choose(m,2) * ( xi_nu - xi_pi )
       kappa_den_3 = m*(nvertex -m)*( xi_phi - xi_pi)
       kappa = kappa_num/sqrt(kappa_den_1+kappa_den_2 + kappa_den_3)

       chart[i]=(pt(quartile, n-1, ncp = kappa , lower.tail = FALSE))
       ##print(c(scale,i,kappa_num^2,kappa_den_1,kappa_den_2,kappa_den_3,kappa))
     }

     my[,row] = chart[1:abnorm]
   }
 }

#  pdf("namasy.pdf")
 plot(x1seq,my[9,],type="l",xlab="theta",ylab="beta")
 abline(h=0.05)
 abline(v=0.465)
 abline(v=0)
 abline(v=0.785)
#  dev.off()
 return(my)
}

