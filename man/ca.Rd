\name{ca}
\alias{ca}
\title{Correspondence analysis}
\description{
  Achieves a Correspondence Analysis (CA) on a numeric table of class data.frame
}
\usage{
ca(x, nfac = 3, isup = 0, jsup = 0, histev = FALSE, grr = FALSE, grc = FALSE, 
grrc = FALSE, grlist = rbind(c(1, 2), c(1, 3), c(2, 3)), prtm = FALSE, 
prtevr = FALSE, prtevc = FALSE, eps = 1e-09)
}
\arguments{
  \item{x}{data.frame minimal dimension 4 x 3. The first column \strong{must contain the character strings of the identifiers of raws}
           any other type, class or dimension results in an error and in the program break.
}
  \item{nfac}{ Number of factors to retain (maximum 7)
}
  \item{isup}{list of illustrative rows. 0 = no illustrative rows (default) 
}
  \item{jsup}{List of illustrative columns. Same as isup.
}
  \item{histev}{Boolean : whether to plot or not the histogram of eigenvalues.
}
  \item{grr}{Boolean : plot the graph of rows on the axes defined by grlist. 
}
  \item{grc}{Boolean : Plot the graph of columns on the axes defined by grlist.
}
  \item{grrc}{Boolean : Plot the simultaneous graph of rows and columns  on the axes defined by grlist. Labels of rows in black, labels of columns in red. 
}
  \item{grlist}{matrix: defines the factorial plans to plot. See details for an example.
}
  \item{prtm}{Boolean: Print or not the data frame. Default = FALSE
}
  \item{prtevr}{Boolean: Print or not the rows eigenvectors. Default = FALSE
}
  \item{prtevc}{ Boolean: Print or not the columns eigenvectors. Default = FALSE
}
  \item{eps}{ numeric: (tolerance) Precision for null eigenvalues. Default = 10E-09
}
}
\details{
\strong{grlist:} the successive plots to draw are defined by a matrix of dimension k,2. k = number of plans to plot. Example: to plot the plans 1-2, 1-3 and 2-3 enter sometning as matrix(1,2,1,3,2,3,nrow=3,ncol=2,byrow=2) or rbind(c(1,2),c(1,3),c(2,3)). 
\strong{Markovian matrix:} In the case of a Markovian or of a transition matrix, one can symetrise (X + t(X)) and load it (sum of the margins added to the diagonal, before applying CA (cf \code{See Also}).\cr   
  In the case of a markovian square matrix (succession or transition matrix) one can symmetrize and load it (\code{\link{symet}}) before representing it by a graph (\code{\link{flux}})
}

\value{
 An object of class ca with attributes
  \item{fr }{data.frame: weight and factorial coordinates of each row (principal and illustrative). The attribute \code{type} has the value "pri" for principal and "ill" for illustrative}
  \item{fc }{data.frame: weight and factorial coordinates of each column (principal and illustrative). \code{type} as in fr}
}
\references{ Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers. }
\author{ Jean-Sebastien Pierre 	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}
\seealso{
  \code{\link{princomp}}, \code{\link{compseq}} to build a transition matrix, \cr \code{\link{symet}} to modify it (symmetrization and diagonal loading), \code{\link{flux}} for the design of a graph.
}
\examples{
# On Csa data (xcsa)
library(sequence)
   data(xcsa)
   ca(xcsa)
}
\keyword{multivariate}
