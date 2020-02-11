\name{aphmat}
\alias{aphmat}
\docType{data}
\title{
Transition matrices of the flight behaviour of 10 individuals of \emph{Aphelinus abdominalis} 
}
\description{
10 succession matrices extracted from the object \code{\link{aphelinus}}. Matrices are provided as \code{data.frames}
}
\usage{data(aphmat)}
\format{
  The format is:
List of 10 \code{data.frames}\cr
 aphmat[[1]] :'data.frame':	5 obs. of  5 variables:\cr
   B : num [1:5] 0 0 0 0 0\cr
   E0: num [1:5] 0 0 0 0 0\cr
   TD: num [1:5] 1 0 0 0 0\cr
   TH: num [1:5] 0 0 1 0 0\cr
    V: num [1:5] 0 0 0 0 0\cr
 and so on until aphmat[[10]]
}

\source{
Personal communication by courtesy of Anne Le Ralec
}
\references{
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. 
Biology of Behaviour 15(3-4): 125-151.(in French) \cr
}
\examples{
data(aphmat)
compmat(aphmat)
}
\keyword{datasets}
