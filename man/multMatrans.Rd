\name{multMatrans}
\alias{multMatrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts a list of sequences into a list of transition matrices
}
\description{
Converts a list of sequences into a list of transition matrices. Matrices are square, of same size,  and provided as \code{data.frame}.
}
\usage{
multMatrans(ser, printdata = FALSE, printdico = TRUE, printmat = FALSE)
}
\arguments{
  \item{ser}{
The list of sequences 
}
  \item{printdata}{
: print  the input sequences 
}
  \item{printdico}{
: Prints the dictionary of items
}
  \item{printmat}{
: Print the list of transition matrices
}
}
\value{
An object of class \code{lmat} 
	\item{dico}{: dictionnary of sequences components}
        \item{listmat}{The list of data.frames} 
}
\references{
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. 
Biology of Behaviour 15(3-4): 125-151.(in french) \cr
}
\author{ Jean-Sebastien Pierre\cr
	\email{Jean-sebastien.pierre@univ-rennes1.fr}\cr
}

\examples{
data(aphelinus)
aphMat=multMatrans(aphelinus)
aphMat$listmat
}

\keyword{ multivariate}
