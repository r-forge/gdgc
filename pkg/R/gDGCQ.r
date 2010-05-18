#     MDGCQ
#
# Description:

#     This function calculates, by Monte Carlo simulation, the upper
#     percentile points of the distribution of the "distance" to the 
#     root node of a binary tree obtained by aglomerative 
#     cluster analysis algorithms applied to the mahalanobis distance 
#     matrix between instances to be classified.
#     It is supposed that there are replicates of p-variate observations
#     for each instance.

# Dependencies: 
#         library(mvtnorm)

# Usage:
#         
#         MDGCQ(N,n,p,linkage,densplot)
# Arguments:

#    N:         simulation cycles
#    n:         vector of length k containing the replicate in each of
#               the k groups
#    p:         number of variables
#    linkage:   "average","single","complete","ward","mcquitty","median"
#    densplot:  boolean argument indicating whether or not, the kernel 
#               estimated density of the distance to the root node 
#               of the binary will be displayed
#    diagonal:  boolean argument indicating whether or not the covariance 
#               matrix is diagonal 
#
#Value:
#     an object of class list containing the upper percentile points of 
#     the distribution of the of the distance to the root node of a binary
#     tree obtained by hclus when applied to the mahalanobis distance matrix
#     between instances to be classified. Percentil points are: 80,85,90,95,99
#     and 99.9%
    
#Author(s):

#     Julio A. Di Rienzo dirienzo@agro.uncor.edu and 
#     Silvia G. Valdano  svaldano@exa.unrc.edu.ar

#References:
#     Di Rienzo, J. A., Guzmán, A. W.  and Casanoves, F. (2002), 
#    "A Multiple Comparisons Method Based on the Distribution of the 
#     Root Node Distance of a Binary Tree," 
#     Journal of Agricultural, Biological and Environment Statistics, 
#     7(2): 129-142.
#
#
#     Valdano S. and Di Rienzo J. (2007). Discovering meaningful groups
#     in hierarchical cluster analysis. An extension to the multivariate
#     case of a multiple comparison method based on cluster analysis.
#     http://interstat.statjournals.net/YEAR/2007/abstracts/0704002.php

# Examples:

#  MDGCQ(500,rep(3,100),5,"average",TRUE)


MDGCQ<-function(N,n,p,linkage,densplot,diagonal)
{
library(mvtnorm)
#begin of Q_MDGC
css<-function(v) {diag((sd(v,na.rm=TRUE))^2,ncol(v)) }

Q_MDGC<-function(n,p,index,linkage)
{
 k<-length(n);
 x<-as.matrix(rmvnorm(sum(n),rep(0,length =p)))
 if (diagonal) {v<-by(x,index,css)} else {v<-by(x,index,var)}
 ##v<-by(x,index,var)
 pooled=as.matrix(as.data.frame(v[1]))%*%diag(rep(n[1]-1,p),p);
 for (i in (2:k)) {pooled=pooled+as.matrix(as.data.frame(v[i]))%*%diag(rep(n[i]-1,p),p)}

 gl<-1/(sum(n)-k)
 nm<-diag(rep(gl,p),p)
 pooled=pooled%*%nm
 pooled<-solve(pooled)

 m<-by(x,index,mean,na.action=na.omit)
 M<-t(as.matrix(as.data.frame(m[1])))
 for (i in (2:k)) {
                  a=t(as.matrix(as.data.frame(m[i])))
                  M=rbind(M,a)
                  }
 D<-matrix(nrow=k,ncol=k)                  
 for (i in (1:k)) {D[,i]<-mahalanobis(M,M[i,],pooled,inverted=T)}
 D<-sqrt(D)
 hclust(as.dist(D),method=linkage)
}
#end of Q_MDGC

Q<-c()
index<-c()
for (i in (1:length(n))) index=c(index,rep(i,length=n[i]))
index=as.factor(index)

for (i in (1:N)){hc<-Q_MDGC(n,p,index,linkage);Q<-c(Q,max(hc$height))}
if (densplot) {plot(density(Q),xlab="Q",main="Estimated density of Q")}
c(quantile(Q,0.80),quantile(Q,0.85),quantile(Q,0.90),quantile(Q,0.95),quantile(Q,0.99),quantile(Q,0.999))
}


                                         

                                           
               

 