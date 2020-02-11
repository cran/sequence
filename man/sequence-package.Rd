\name{sequence-package}
\alias{sequence-package}
\docType{package}
\title{
Analysis of behavioral (or of any other type) sequences.
}
\description{
Builds transition/succession matrices (lag analysis), compares several matrices. Plane drawing of flow charts obtained from those matrices by a method based on a modified Correspondence Analysis (CA). The modification consists in symmetrising the matrix and loading its diagonal following Foucart (1985) and Van der Heijden (1986). This leads to a 'nice' positionning of the edges of the graph, representing the different items in succession. 'Smart' arrows are drawn as vertices.
}
\details{
\code{\link{read.series.seq}} reads a series of sequences, each of them considered as an identifiers list. \code{\link{compseq}} compares the sequences by likelihood ratio tests (LRT) with a Bonferroni correction. The LRT serve to build a distance matrix of the sequences two-by-two. This matrix of distance is used both for a clustering of them and for a mapping by Principal Coordinates Analysis (PCA) also called Multidimensional Scaling (MDS). 

From the sequences a transition matrix is obtained. \code{\link{symet}} symmetrises and/or loads its diagonal. \code{\link{ca}} makes a Correspondence Analysis and saves the factorial coordinates as an object. \code{\link{flux}} draws the transition graph. The complete analysis chain is the following:\cr
	read.series.seq -> compseq -> symet -> ca -> flux \cr
A sketch of the method is given in Pierre and Kasper (1990).  
}
\author{ Jean-Sebastien Pierre \email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\references{
Pierre, J. S. and C. Kasper (1990). "The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship." Biology of Behaviour 15(3-4): 125-151.\cr
Foucart, T. (1985). "Tableaux symetriques et tableaux d'echanges." Revue statistiques appliquee 33: 37-53.\cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\keyword{ package }
