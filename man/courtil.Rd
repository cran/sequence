\name{courtil}
\alias{courtil}
\docType{data}
\title{ Sexual courtship of the male molecricket }
\description{
  Succession matrix of the patterns of male molecricket (\emph{Gryllotalpa gryllotalpa} orthoptera)
}
\usage{data(courtil)}
\format{
  A data frame with 28 observations on the following 29 variables.
  \describe{
    \item{\code{id}}{a factor with levels \cr 
\code{AA} \code{AA+} \code{AAB} \code{ACC} \code{ACQ} \code{BAL} \code{CH} \code{CP1} \code{DIV} \code{E} \code{FR} \code{FUI} \code{IMM} \code{INV} \code{LA} \code{M} \code{MB} \code{P1} \code{PLP} \code{PRJ} \code{R} \code{RET} \code{RV} \code{SP1} \code{TR} \code{VAV} \code{VC} \code{VCH}}
    \item{\code{DIV}}{diverse (walking, feeding, grooming)}
    \item{\code{VAV}}{go to}
    \item{\code{AA.}}{lowers wings}
    \item{\code{ACQ}}{antenna-cerques contact}
    \item{\code{P1}}{rise P1}
    \item{\code{CP1}}{kicking with forelegs (P1)}
    \item{\code{SP1}}{joins and rises forelegs }
    \item{\code{LA}}{rises antennae and rehects them backward}
    \item{\code{PRJ}}{projection of liquid from anus}
    \item{\code{RET}}{returns}
    \item{\code{FUI}}{escapes}
    \item{\code{IMM}}{freezes}
    \item{\code{R}}{steps back}
    \item{\code{BAL}}{rocking}
    \item{\code{FR}}{wings rustling}
    \item{\code{TR}}{short and acute sound ("treets")}
    \item{\code{CH}}{call song}
    \item{\code{PLP}}{palpation of the female}
    \item{\code{E}}{avoids}
    \item{\code{INV}}{turns around and approaches}
    \item{\code{AA}}{lowers wings at maximum}
    \item{\code{VCH}}{vibers and sings}
    \item{\code{AAB}}{lowers wings and abdomen}
    \item{\code{VC}}{vibration of cerques}
    \item{\code{RV}}{steps back lowering wings and abdomen}
    \item{\code{M}}{inserts backward under the female which climbs on the mail}
    \item{\code{ACC}}{mating}
    \item{\code{MB}}{'masturbation' - rubbing of genitalia with p3 tarsa}
  }
}
\details{
  Matrix of succession of patterns, provided as a \code{data.frame}.  The code of the rows, character values of the levels of the factor \code{id} are the same as the names of the columns after \code{id}. The elements of the matrix. The value of each element reflects the number of times the pattern \code{id[i]} is followed by the pattern names\code{(courtil)[j+1]}. 
}
\references{
Pierre, J. S. and C. Kasper (1990). The Design of Ethological Flow-Charts on Factorial Analysis Representations - an Application to the Study of the Male Mole-Cricket Sexual Courtship. Biology of Behaviour 15(3-4): 125-151.(in french) \cr
Van der Heijden, P. G. M. 1986. Transition matrices, model fitting and correspondence analysis. In: Data Analysis and Informatics IV (Ed. by E. Diday), pp. 221-226. Elsevier Science Publishers.
}
\examples{
data(courtil)
symet(courtil,ident=TRUE)->scourtil
ca(scourtil,grr=TRUE)
}
\keyword{datasets}
