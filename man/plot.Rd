\name{plot}
\alias{plot.GPvam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for GPvam
}
\description{
Plot teacher effects and residuals.
}
\usage{
\method{plot}{GPvam}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class GPvam
}
  \item{\dots}{
other arguments
}

}

\value{
returns residual plots
}
\references{
Karl, A., Yang, Y. and Lohr, S. (2012) Efficient Maximum Likelihood Estimation of the Generalized Persistence Value-Added Model.\emph{Submitted.} 

Mariano, L., McCaffrey, D. and Lockwood, J. (2010) A Model for Teacher Effects From Longitudinal Data Without Assuming Vertical Scaling. \emph{Journal of Educational and Behavioral Statistics} \bold{35}, 
253--279.

McCaffrey, D. and Lockwood, J. (2011) Missing Data in Value-Added Modeling of Teavher Effects," \emph{Annals of Applied Statistics} \bold{5}, 773--797
}
\author{
Andrew Karl \email{akarl@asu.edu}
Yan Yang
Sharon Lohr
}


\seealso{
 \code{\link{summary.GPvam}}
}
\examples{
result<-GPvam(vam_data)
\dontrun{
plot(result)
}
}
\keyword{ regression }

