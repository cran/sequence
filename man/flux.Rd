\name{flux}
\alias{flux}
\title{  Draws a flow chart from a transition matrix }
\description{
 Draws a flow chart from a transition matrix on the projection plans of a CA (Correspondence Analysis). requires two objects of class \code{data.frame}: \code{tabcoord} a table of axis coordinates on up to 7 columns, and \code{tabtr} the transition matrix.
}
\usage{
flux(tabcoord, tabtr, dic = NULL, fac = c(1, 2), Sort = NULL, threshold = 0, 
scale = 0.1, cscale = 0.1,main="Transition graph", ...)
}
\arguments{
  \item{tabcoord}{ data.frame: nfac + 3 columns 
	\itemize{ 
                 \item{id}{ factor: identifier of individuals }
		  \item{w}{ numeric: weight }
                  \item{f1}{ numeric: coordinates on axis 1 }
                  \item{f2,...}{numeric: next coordinates until nfact }
                  \item{type}{ factor: "pri" for principal, "ill" for illustrative (or 'supplementary') }
}
}
  \item{tabtr}{ data.frame: transition matrix. Either frequencies or probabilities. \strong{Caution} this data.frame must contain the row identifiers as first column }
  \item{dic}{ character vector: list of items labels. If NULL (default) must be extractible from \code{tabtr} column names }
  \item{fac}{ vector: two elements: the numbers of the two factors to use in \code{tabcoord} }
  \item{Sort}{ numeric: number of the factor used to sort the drawing of the elements.\cr Default = NULL (sorting in the order of the dictionary. Otherwise must be comprised between 1 et nfac. \cr If positive, sort in ascending order. If negative, sort in decreasing order  }
  \item{threshold}{numeric: minimal threshold to draw the arrows. Represent a oroportion of the weight of the element used as origin of the arrow. }
  \item{scale}{numeric: Controls the width of the arrows  }
  \item{cscale}{numeric: Controls the radius of the circles.  }
  \item{main}{character:title of the plot}
  \item{\dots}{ \code{\dots} supplementary arguments for \code{\link{smartArrow}} }
}
\note{ The correct definition of fac is under the user's responsability. \strong{Must know} how many coordinates exist in the object tabcoord. Use \code{summary(<tabcoord>)} in case of doubt.
}
\details{
  \code{tabcoord} must be an the attribute $fl or $fc of an object of class \code{\link{ca}}
}
\value{
An object of class \code{graph}.
}
\references{  
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.(in french) \cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\author{ Jean-Sebastien Pierre\cr 
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
%	 \url{http://www.perso.univ-rennes1.fr/jean-sebastien.pierre}
}
\seealso{\code{\link{compseq}} to create transition matrices, \code{\link{compmat}} to compare transition matrices, \cr\code{\link{read.series.seq}} to read a file of sequences,\code{\link{ca}} to perform a Correspondence Analysis }
\examples{

# Analysis of Wiepkema's data
data(wiepkema)
# Symmetrisation and diagonal loading
symet(wiepkema,ident=TRUE)->swiepkema
# Correspondence Analysis
ca(wiepkema,isup=c(5,11),jsup=c(5,11))->afcwiepk
# Flow chart / transition graph
flux(tabcoord=afcwiepk$fr,tabtr=wiepkema,threshold=0.025,scale=2)

}
\keyword{ multivariate }
