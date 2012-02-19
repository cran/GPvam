\name{GPvam}
\alias{GPvam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting the Generalized Persistence Model
}
\description{
GPvam is used to fit the value-added model developed by Mariano et al. (2010) via ML estimation using an EM algorithm (Karl et al. 2012).
}
\usage{
GPvam(vam_data, fixed_effects = formula(~as.factor(year) + 0), 
   student.side = "R", max.iter.EM = 1000, tol1 = 1e-07, 
   hessian = TRUE, hes.method = "simple", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vam_data}{
a data frame that contains at least a column "y" containing the student scores, a 
column "student" containing unique student ID's, a column "teacher" containing the 
teacher ID's, and a column "year" which contains the year (or semester, etc.) of the time period. The "y" and "year" variables needs to be numeric. If other variables are to be included 
as fixed effects, they should also be included in vam_data. See 'Note' for further discussion.
}
  \item{fixed_effects}{
an object of class \code{formula} describing the structure of the fixed effects. Categorical variables
should be wrapped in an \code{as.factor} statement.
}
  \item{student.side}{
a character. Choices are "G" or "R". See section 'Details'.
}


  \item{max.iter.EM}{
the maximum number of EM iterations
}
  \item{tol1}{
convergence tolerance for EM algorithm. The convergence criterion is specified under 'Details'.
}
  \item{hessian}{
logical indicating whether the Hessian of the variance parameters should be calculated after convergence of the EM algorithm. Standard errors for the fixed and random effects are calculated by default.
}
  \item{hes.method}{
a character string indicating the method of numerical differentiation used
            to calculate the Hessian of the variance parameters. Options are "simple" or "richardson".
}
  \item{verbose}{
logical. If TRUE, model information will be printed at each iteration.
}
}
\details{
The design for the random teacher effects according to the generalized persistence model
of Mariano et al. (2010) is built into the function. The model includes correlated current- and
future-year effects for each teacher. By setting \code{student.side="R"}, the intra-student correlation is 
modeled via an unstructured, block-diagonal error covariance matrix, as specified by Mariano et al. (2010). Setting 
\code{student.side="G"} keeps the same teacher structure, but models intra-student correlation via 
random student effects. This is similar to the model used by McCaffrey and Lockwood (2011), and is appropriate when the testing 
scale is the same across years. In this case, the error covariance matrix is diagonal, although a separate variance is calcualted for each year.
From a computational perspective, the model estimating the R-side student effects has better scalability properties, although the G-side function is faster (Karl et al. 2012).

Convergence is declared when \eqn{(l_k-l_{k-1})/l_k < 1E-07}, where \eqn{l_k} is the log-likelihood 
at iteration \eqn{k}.


The model is estimated via an EM algorithm. For details, see Karl et al. (2012). The model was estimated through Bayesian computation in Mariano et al. (2010).

Note: When student.side="R" is selected, the first two iterations of the EM algorithm will take longer than subsequent iterations. This is a result of the hybrid gradient-ascent/Newton-Raphson method used in the M-step for the R matrix in the first two iterations (Karl et al. 2012).
}
\value{
GPvam returns an object of class \code{GPvam}



An object of class \code{GPvam} is a list containing the following components:

\item{loglik}{the maximized log-likelihood at convergence of the EM algorithm}
\item{teach.effects}{a data frame containing the predicted teacher effects and standard errors}
\item{parameters}{a matrix of estimated model parameters and standard errors}
\item{Hessian}{if requested, the Hessian of the variance parameters}
\item{R_i}{(only when \code{student_side} is set to 'R') a matrix containing the error covariance matrix of a student}
\item{teach.cov}{a list containing the unique blocks of the covariance matrix of teacher effects}
\item{mresid}{a vector of the raw marginal residuals}
\item{cresid}{a vector of the raw conditional residuals}
\item{sresid}{a vector of the scaled conditional residuals}
\item{yhat}{a vector of the predicted values}

The function  \code{\link{summary}} provides a summary of the results. This includes the estimated model parameters and standard errors, along with the correlation matrices corresponding to the estimated correlation matrices. Summary information about scaled and raw residuals is reported.
}
\references{
Karl, A., Yang, Y. and Lohr, S. (2012) Efficient Maximum Likelihood Estimation of the Generalized Persistence Value-Added Model.\emph{Submitted.} 

Mariano, L., McCaffrey, D. and Lockwood, J. (2010) A Model for Teacher Effects From Longitudinal Data Without Assuming Vertical Scaling. \emph{Journal of Educational and Behavioral Statistics} \bold{35}, 
253--279.

McCaffrey, D. and Lockwood, J. (2011) Missing Data in Value-Added Modeling of Teacher Effects," \emph{Annals of Applied Statistics} \bold{5}, 773--797

}



\author{
Andrew Karl \email{akarl@asu.edu}, Yan Yang, Sharon Lohr
}
\note{
The model assumes that each teacher teaches only one year. If, for example, a teacher teaches in years 1 and 2, his/her first year performance is modeled independently of the second year performance.
To keep these effects seperate, the progam appends "(year i)" to each teacher name, where i is the year in which the teacher taught.

The \code{fixed_effects} argument of \code{GPvam} utilizes the functionality of R's \code{formula} class. In the statement \code{fixed_effects=formula(~as.factor(year)+cont_var+0))}, \code{as.factor(year)} identifies \code{year} as a categorical variable. \code{+0} indicates that no intercept is to be fitted, and \code{+cont_var} indicates that a seperate effect is to be fitted for the continuous variable "cont_var." An interaction between "year" and "cont_var" could be specified by \code{~as.factor(year)*cont_var+0}, or equivalently, \code{~as.factor(year)+cont_var+as.factor(year):cont_var+0}. See \code{\link{formula}} for more details.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{plot}},
 \code{\link{summary}},
 \code{\link{vam_data}}
}
\examples{
data(vam_data)
 result <- GPvam(vam_data,student.side="G",
fixed_effects=formula(~as.factor(year)+cont_var+0))
 summary(result)
\dontrun{
 plot(result)
}
}

\keyword{ regression }

