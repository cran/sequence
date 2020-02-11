\name{compseq}
\alias{compseq}
\title{ Comparison of behavioral (or any) sequences }
\description{
  Compares statistically sequences of states (behavior, texts, molecular data) by likelihood ratio tests on their markovian transition matrices. Performs also a cluster analysis of the sequences and a Principal Coordinates Analysis on the distance matrix between them.
}
\usage{
compseq(ser,alpha=0.05,meth="ward.D",printdata=FALSE,printdico=TRUE,printmat=FALSE,
eps=1e-07,clust=TRUE,pca=TRUE)
}
\arguments{
  \item{ser}{ list of list: set of sequences }
  \item{alpha}{ numeric: global risk threshold for pairwise comparisons. }
  \item{meth}{ character:Clustering method. cf \code{hclust}. }
  \item{printdata}{ Boolean:Print original data.}
  \item{printdico}{ Boolean:Print the dictionnary of states from \code{ser}.}
  \item{printmat}{ Boolean: print all transition matrices and the consensus matrix. }
  \item{eps}{ numeric: precision for the convergence of \code{cmdscale}.}
  \item{clust}{ Boolean: do the cluster analysis.}
  \item{pca}{ Boolean: do the principal coordinates analysis.}
}
\details{
The log likelihood ratio times -2 is used both for tests (Chi-Square approximation followed by Bonferroni post hoc tests) and as a distance to cluster the sequences and to represent them on factorial plans (Principal Coordinates Analysis). Warning: not a metric distance. Susceptible to give incoherent clustering with some methods (\code{meth}).\cr
This function does essentially the same work as \code{\link{compmat}} but with matrices instead of sequences entry.  
}
\value{
an object of class compseq with attributes
\itemize{
 \item{dico}{ Dictionnary of states }
 \item{mdist}{ Matrix of pairwise distances between sequences }
 \item{msign}{ Matrix of pairwise significance levels between sequences }
 \item{mcom}{ Common or consensus transition matrix}}
}

\references{
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.\cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\author{ Jean-Sebastien Pierre\cr
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{ \code{\link{hclust}}, \code{\link{cmdscale}}, \code{\link{ca}}, \code{\link{flux}} }
\examples{
	data(seriseq)
	compseq(seriseq)
}
\keyword{ multivariate }

