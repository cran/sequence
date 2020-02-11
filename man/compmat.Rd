\name{compmat}
\alias{compmat}
\title{
Comparison of transition (or succession) matrices
}
\description{
 Compares statistically succession matrices by likelihood ratio tests. Performs also a cluster analysis of the sequences and a Pricipal Coordinates Analysis (PCA) on the distance matrix between them.
}
\usage{
compmat(serMat, alpha = 0.05, meth = "ward.D", printdata = FALSE, printdico = FALSE, 
printmat = FALSE, eps = 1e-07, clust = TRUE, pca = TRUE)
}
\arguments{
  \item{serMat}{ List of data.frames. Each of them must contain a matrix of identical dimension, with the same row.names.}
  \item{alpha}{ numeric: global risk threshold for pairwise comparisons. Default = 0.05}
  \item{meth}{character: Clustering method. cf \code{hclust}.}
  \item{printdata}{ Boolean: Print original list of matrices.}
   \item{printdico}{ Boolean:Print the dictionnary of states from \code{ser}.}
  \item{printmat}{ Boolean: print all transition matrices and the consensus matrix.}
  \item{eps}{ numeric: precision for the convergence of \code{cmdscale}.}
  \item{clust}{ Boolean: performs cluster analysis.}
  \item{pca}{ Boolean: performs a Principal Coordinates Analysis.}
}
\details{
The log likelihood ratio times -2 is used both for tests (Chi-Square approximation followed by Bonferroni post hoc tests) and as a distance to cluster the sequences and to represent them on factorial plans (Principal Coordinates Analysis). Warning: not a metric distance. Susceptible to give incoherent clustering with some methods (\code{meth}). 
}
\value{
an object of class compseq with attributes
\itemize{
 \item{dico}{ Dictionnary of states }
 \item{mdist}{ Matrix of pairwise distances between sequences }
 \item{msign}{ Matrix of pairwise significance levels between sequences }
 \item{mcom}{ Common or consensus transition matrix}}
}

\references{a
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.\cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\author{ Jean-Sebastien Pierre \email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{hclust}}, \code{\link{cmdscale}}, \code{\link{ca}} }
\examples{
# Compares 10 transition matrices in \code{aphmat} 
data(aphmat)
compmat(aphmat,clust=FALSE,pca=FALSE)

}
\keyword{ multivariate }

