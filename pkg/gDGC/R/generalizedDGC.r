# generalizedDGC 

# Description:

#        This function generates a dendrogram representation of
#        group mean vectors according to Mahalanobis distance
#        In case the dimension of data is 1, the distances are just the
#        Euclidean distances using scaled atributes (divided by their common sd) 

# Dependencies: 

#         library(mvtnorm)
#         function generalizedDGCQ

# Usage:
#         
# usage     generalizedDGC(x,linkage.method="average",diagonal=FALSE,B=500,nivel=5)

# Arguments:
#
# data:     data frame containing the data. It is supossed
#           that its first column is a factor identifying the groups to
#           which the observations belong. The rest of the columns 
#           contains the p-variates

#           It is supossed that there are replication for each group
#           although it admits different number of replications per group.
#           The number of replications should be enough to garantee the
#           non singularity of the pooled covariance matrix

# linkage.method: 
#            must be one of the following '"ward"', '"single"',
#           '"complete"', '"average"', '"mcquitty"', '"median"' or
#           '"centroid"'. Default="average"

# diagonal: boolean. If TRUE the covariance matrix is supposed to be diagonal 

# B:        integer. Is the number of times the null distribution of root 
#           node of the binary tree will be sample to obtain critical points. 
#           Default=1500
#           
# plot:     boolean indicating if a plot of the dendrogram a the cutting point
#           should be displayed.

#Value:
#     an object of class data frame containing the value of the parameters used
#     and the list of the upper percentile points of 
#     the distribution of the distance to the root node of a binary
#     tree obtained by hclus when applied to the mahalanobis distance matrix
#     between instances to be classified. Percentil points are: 80,85,90,95,99
#     and 99.9%
    
#Author(s):

#     Julio A. Di Rienzo dirienzo@agro.uncor.edu 

#References:
#     Di Rienzo, J. A., Guzmán, A. W.  and Casanoves, F. (2002), 
#    "A Multiple Comparisons Method Based on the Distribution of the 
#     Root Node Distance of a Binary Tree," 
#     Journal of Agricultural, Biological and Environment Statistics, 
#     7(2): 129-142.
#
#     Valdano S. and Di Rienzo J. (2007). Discovering meaningful groups
#     in hierarchical cluster analysis. An extension to the multivariate
#     case of a multiple comparison method based on cluster analysis.
#     http://interstat.statjournals.net/YEAR/2007/abstracts/0704002.php
 
# Examples:

# generalizedDGC(data,linkage.method="average",FALSE,500,TRUE,level)
# 
 
generalizedDGC<-function(data,linkage.method="average",diagonal=FALSE,B=1500,nivel=4)
{
#source("gDGCQ.r") 

n=c(as.data.frame(summary(as.factor(data[,1])))[,1])
p=ncol(data)-1;
k=length(levels(as.factor(data[,1])))
x=c()
x=data[,2:ncol(data)];
index=as.factor(data[,1])
v=by(x,index,var)
nm=diag(rep(n[1]-1,p),p)
pooled=as.matrix(as.data.frame(v[1]))%*%nm
#pooled
for (i in (2:k)) 
{
nm=diag(rep(n[1]-1,p),p)
pooled=pooled+as.matrix(as.data.frame(v[i]))%*%nm
}

gl=1/(sum(n)-k)
nm=diag(rep(gl,p),p)
pooled=pooled%*%nm

if (diagonal==TRUE) pooled=diag(diag(pooled))

pooled<-solve(pooled)

m=by(x,index,mean)
M=t(as.matrix(as.data.frame(m[1])))
for (i in (2:k)) {
                  a=t(as.matrix(as.data.frame(m[i])))
                  M=rbind(M,a)
                  }

D=matrix(nrow=k,ncol=k)                  
for (i in (1:k)) {D[,i]<-mahalanobis(M,M[i,],pooled,inverted=T)}
D=sqrt(D)
rownames(D)=levels(index)
hc=hclust(as.dist(D),method=linkage.method)
q<-gDGCQ(N=B,n,p,linkage=linkage.method,densplot=FALSE,diagonal=diagonal)
res<-as.data.frame(c(linkage.method,diagonal,B,as.matrix(q)))
rownames(res)=c("Linkage","Diagonal","B",names(q))
colnames(res)='-------------'
resultados=as.data.frame(cbind(as.factor(hc$labels),cutree(hc,h=q[nivel])))
resultados[,1]=rownames(resultados)
resultados=list("tabla"=resultados[order(resultados[,2]),],"q"=q,"hc"=hc)
#list("q"=q,"hc"=hc,"id"=cutree(hc,h=q[nivel]))
}
