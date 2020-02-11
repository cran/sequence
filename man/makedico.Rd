\name{makedico}
\alias{makedico}
\title{ Build the dictionary of the sequences   }
\description{
  Assesses the dictionnary common to various sequences, counts the occurences of each word (item, state) 
}
\usage{
makedico(x, printdata = FALSE, printdico = FALSE)
}
\arguments{
  \item{x}{ sequence: an object of class sequence }
  \item{printdata}{ \code{printdata} logical: print sequences or not.}
  \item{printdico}{ \code{printdico} logical: print or not the resulting dictionnary.}
}
\value{A list of character strings}
\author{Jean-Sebastien Pierre email: Jean-sebastien.pierre@univ-rennes1.fr} 
\examples{
data(seriseq)
dic=makedico(seriseq,printdico=TRUE)
}
\keyword{ multivariate }
