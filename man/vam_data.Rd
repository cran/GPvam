\name{vam_data}
\alias{vam_data}
\docType{data}
\title{
Simulated Data
}
\description{
A small, simulated unbalanced data set used to illustrate the functionality of the package.
}
\usage{data(vam_data)}
\format{
  A data frame with 720 observations on the following 5 variables.
  \describe{
    \item{\code{y}}{a numeric vector representing the student score}
    \item{\code{student}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{teacher}}{a numeric vector}
    \item{\code{cont_var}}{a numeric vector representing a continuous covariate}
  }
}

\examples{
data(vam_data)
print(vam_data[1,])
}
\keyword{datasets}
