\name{aphelinus}
\alias{aphelinus}
\docType{data}
\title{ Flight behaviour of \emph{Aphelinus abdominalis}}
\description{
Flight behaviour of the parasitic hymenoptera \emph{Aphelinus abdominalis} in flight tunnel.
}
\usage{data(aphelinus)}
\format{
  The format is:
List of 10 sequences on the following pattern:\cr
 aphelinus[[1]] =list("\%1","B","TD","TH")\cr
 and so on until aphelinus[[10]]
}
\details{ A sequences of 10 individuals}
\source{Anne Le Ralec, personnal communication.}
\references{Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.(in French)}

\examples{
data(aphelinus)
u=compseq(aphelinus)
mat=u$mcom
print(mat)
}
\keyword{datasets}
