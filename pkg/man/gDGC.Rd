\name{gDGC}
\alias{gDGC}
\title{gDGC}
\description{This function generates a dendrogram representation of group mean vectors according to Mahalanobis distance.
 	In case the dimension of data is 1, the distances are just the Euclidean distances using scaled atributes (divided by their common sd) 
\usage{
gDGC(data,linkage.method="average",diagonal=FALSE,B=500,level=5)
}
\arguments{
  \item{data}{
	Data frame containing the data. It is supossed that its first column is a factor identifying the groups to
        which the observations belong. The rest of the columns contains the p-variates.
        It is supossed that there are replication for each group although it admits different number of replications per group.
        The number of replications should be enough to garantee the non singularity of the pooled covariance matrix
}
  \item{linkage.method}{
	It must be one of the following '"ward"', '"single"','"complete"', '"average"', '"mcquitty"', '"median"' or '"centroid"'. 
	Default="average"
}
  \item{diagonal}{
	It is a boolean. If TRUE the covariance matrix is supposed to be diagonal. 

}
  \item{B}{
	Integer that represents the number of times the null distribution of root node of the binary tree will be sample to obtain 
	critical points. Default=1500
}

 \item{level}{
	Boolean indicating if a plot of the dendrogram a the cutting point should be displayed.

}
}
\details{

}
\value{

	An object of class data frame containing the value of the parameters used and the list of the upper percentile points of 
	the distribution of the distance to the root node of a binary tree obtained by hclus when applied to the mahalanobis 
	distance matrix between instances to be classified. Percentil points are: 80,85,90,95,99 and 99.9%

}
\references{
	Di Rienzo, J. A., Guzmán, A. W.  and Casanoves, F. (2002), "A Multiple Comparisons Method Based on the Distribution of 
	the Root Node Distance of a Binary Tree," Journal of Agricultural, Biological and Environment Statistics, 7(2): 129-142.

	Valdano S. and Di Rienzo J. (2007). Discovering meaningful groups in hierarchical cluster analysis. 
	An extension to the multivariate case of a multiple comparison method based on cluster analysis.
	http://interstat.statjournals.net/YEAR/2007/abstracts/0704002.php
 
}
\author{
	Julio A. Di Rienzo  <email: dirienzo@agro.unc.edu.ar>
	Fernando Casanoves  <email: casanoves@catie.ac.cr>
       }
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{gDGCQ}}
}
\examples{
	gDGC(data,linkage.method="average",FALSE,500,TRUE,level)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
