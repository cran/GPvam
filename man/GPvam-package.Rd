\name{GPvam-package}
\alias{GPvam-package}
\docType{package}
\title{
Maximum Likelihood Estimation of the Generalized Persistence Model
}
\description{
An EM algorithm is used to estimate the generalized persistence model
}
\details{
\tabular{ll}{
Package: \tab GPvam\cr
Type: \tab Package\cr
Version: \tab 1.0-0\cr
Date: \tab 2012-02-19\cr
License: \tab GPL-2\cr
}
how to use package
}
\author{
Andrew Karl, Yan Yang, and Sharon Lohr

Maintainer: Andrew Karl <akarl@asu.edu>
}
\references{
Karl, A., Yang, Y. and Lohr, S. (2012) Efficient Maximum Likelihood Estimation of the Generalized Persistence Value-Added Model.\emph{Submitted.} 

Mariano, L., McCaffrey, D. and Lockwood, J. (2010) A Model for Teacher Effects From Longitudinal Data Without Assuming Vertical Scaling. \emph{Journal of Educational and Behavioral Statistics} \bold{35}, 
253--279.

McCaffrey, D. and Lockwood, J. (2011) Missing Data in Value-Added Modeling of Teavher Effects," \emph{Annals of Applied Statistics} \bold{5}, 773--797
}

\examples{
result<-GPvam(vam_data)
summary(result)
\dontrun{
plot(result)
}
}
