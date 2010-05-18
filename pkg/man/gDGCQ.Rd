\name{gDGCQ}
\alias{gDGCQ}
\title{gDGCQ}
\description{This function calculates, by Monte Carlo simulation, the upper percentile points of the distribution of the "distance" 
	to the root node of a binary tree obtained by aglomerative cluster analysis algorithms applied to the mahalanobis distance 
	matrix between instances to be clasSified.
	It is supposed that there are replicates of p-variate observations for each instance.

\usage{
	MDGCQ(N,n,p,linkage,densplot,diagonal)
}
\arguments{
  \item{N}{
	Simulation cycles
}
  \item{n}{
	Vector of length k containing the replicate in each of the k groups
}
  \item{p}{
	Number of variables

}
  \item{linkage}{
	"average","single","complete","ward","mcquitty","median"
}

 \item{densplot}{
	Boolean argument indicating whether or not, the kernel estimated density of the distance to the root node of the binary 
	will be displayed

}
 \item{diagonal}{
	Boolean argument indicating whether or not the covariance matrix is diagonal 

}
}
\details{

}
\value{

	An object of class list containing the upper percentile points of the distribution of the distance to the root node of a 
	binary tree obtained by hclus when applied to the mahalanobis distance matrix between instances to be classified. 
	Percentil points are: 80,85,90,95,99 and 99.9%

}
\references{
	Di Rienzo, J. A., Guzmán, A. W.  and Casanoves, F. (2002), "A Multiple Comparisons Method Based on the Distribution of 
	the Root Node Distance of a Binary Tree" Journal of Agricultural, Biological and Environment Statistics, 7(2): 129-142.

	Valdano S. and Di Rienzo J. (2007). Discovering meaningful groups in hierarchical cluster analysis. An extension to the 
	multivariate case of a multiple comparison method based on cluster analysis.
	http://interstat.statjournals.net/YEAR/2007/abstracts/0704002.php
 
}
\author{
	Julio A. Di Rienzo dirienzo@agro.uncor.edu and 
        Silvia G. Valdano  svaldano@exa.unrc.edu.ar
	Fernando Casanoves  <email: casanoves@catie.ac.cr>
       }
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	
}
\examples{
	MDGCQ(500,rep(3,100),5,"average",TRUE,TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
