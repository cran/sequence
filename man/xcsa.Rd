\name{xcsa}
\alias{xcsa}
\docType{data}
\title{ French CSA (1976) data}
\description{
  Medical attendence depending on socio-professional categories. 
}
\usage{data(xcsa)}
\format{
  A data frame with 15 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{factor: Socio-professional categories: levels: \code{acm} \code{fcm} \code{fco} \code{fem} \code{foo} \code{foq} \code{fpl} \code{fse} \code{hco} \code{hem} \code{hma} \code{hoq} \code{hos} \code{hpl} \code{hse} \code{sup} \code{hcm} \code{fos}}
    \item{\code{ge1}}{ visited a general practitioner during the previous year}
    \item{\code{ge0}}{ did not visit a general practitioner during the previous year}
    \item{\code{sp1}}{ visited a specialist during the previous year}
    \item{\code{sp0}}{ did not visite a specialist during the previous year}
    \item{\code{de1}}{ visited a dentist during the previous year}
    \item{\code{de0}}{ did not visite a dentist during the previous year}
  }
}
\details{ Detail of socio-professional categories codes :  \cr
	First letter = h 	:	man \cr
	First letter = f	:	woman \cr
		cm		:	middle manager \cr
		co		:       manager, unprecised status \cr
		pl		:	Senior executive or liberal profession 	 \cr
		oq		: 	qualified worker \cr
		se		:	unemployed \cr
		os		:	specialised worker \cr
		ma		:	unskilled worker
		
}
\source{
  Ronald Cehessat (1976) exercices de statistique et informatique appliquees Dunod (Paris)
}
\examples{
data(xcsa)
ca(xcsa)
}
\keyword{datasets}
