\name{makeSeries}
\alias{makeSeries}
\title{
Builds a list of list from a vector and a factor
}
\description{
Builds an object of class list of lists, from two columns, of a \code{data.frame}. The first one contains the succession of codes to use, the second is a factor adequate to split the series into separate bouts.
}
\usage{
makeSeries(x, fac)
}
\arguments{
  \item{x}{\code{character vector}: the record of item codes under study
}
  \item{fac}{
\code{factor}: a factor, of same length as x, indicating how to split the sequences into separate bouts.
}
}
\value{
a list of lists, one per level of the factor \code{fac}
}
\details{
Suppose we observe the succession a,c,d,a,a,f on the subject A and b,d,f,e,c,c,a,d on the subject B. This can be coded as the concatenation \cr
  x: a,c,d,a,a,f,b,d,f,e,c,c,a,d \cr
fac: A,A,A,A,A,A,B,B,B,B,B,B,B,B \cr

fac indicating that the 6 first itms were successively emitted by A and the 8 following by B. See example.
}
\author{
Jean-Sebastien Pierre
}


\seealso{
\code{\link{read.old.seq}}, \code{\link{read.series.seq}},\code{\link{compseq}}
}
\examples{

ser=sample(letters[1:6],100,replace=TRUE)
fact=as.factor(c(rep("one",30),rep("two",20),rep("three",10),rep("four",40)))
print(rbind(ser,fact)) # Facultative to see the structure
seq=makeSeries(ser,fact)
print(seq)
compseq(seq) # Facultative, compares the sequences obtained.
}
\keyword{ multivariate }
