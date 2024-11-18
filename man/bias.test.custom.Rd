\name{bias.test.custom}
\alias{bias.test.custom}
\title{
Permutation Tests for Fixed Effects Bias Assessment
}
\description{
Performs permutation tests on fixed effects within a linear mixed model to assess the bias of fixed effect parameters or contrasts. The function allows for both standard basis vectors and custom vectors to define the effects being tested. See Karl and Zimmerman (2021) <doi:10.1016/j.jspi.2020.06.004>.
}
\usage{
bias.test.custom(result, 
                 k_vectors = NULL, 
                 n_perms = 1e5, 
                 save_plots = FALSE, 
                 output_folder = ".")
}
\arguments{
  \item{result}{
    An object containing GPvam results, including the fixed effects matrix (\code{X}), random effects design matrix (\code{Z}), inverse variance matrix (\code{vinv}), estimated random effects (\code{eta.hat}), variance components matrix (\code{G}), number of teachers per group (\code{num.teach}), and persistence type (\code{persistence}). The object must contain the following components:
    \describe{
      \item{\code{X}}{Fixed effects matrix.}
      \item{\code{Z}}{Random effects design matrix.}
      \item{\code{vinv}}{Inverse variance matrix.}
      \item{\code{eta.hat}}{Estimated random effects.}
      \item{\code{G}}{Variance components matrix for random effects.}
      \item{\code{num.teach}}{Vector indicating the number of teachers (random effects) per group.}
      \item{\code{persistence}}{Persistence type, must be either \code{"CP"} or \code{"VP"} or \code{"ZP"}.}
    }
  }
  
  \item{k_vectors}{
    (Optional) A list of numeric vectors specifying custom \eqn{k} vectors for combined fixed effects. Each vector should be the same length as the number of fixed effects in the model. If \code{NULL}, the function generates standard basis vectors (one-hot vectors) to perform permutation tests for each fixed effect individually.
  }
  
  \item{n_perms}{
    (Optional) The number of permutations to perform for each \eqn{k} vector. A higher number of permutations increases the accuracy of the p-value estimates but also increases computation time. Default is \code{1e5}.
  }
  
  \item{save_plots}{
    (Optional) Logical flag indicating whether to save the permutation test histograms as PNG files. If set to \code{TRUE}, histograms for each \eqn{k} vector will be saved in the specified \code{output_folder}. Default is \code{FALSE}.
  }
  
  \item{output_folder}{
    (Optional) Directory path where plots and results will be saved if \code{save_plots} is \code{TRUE}. The directory will be created if it does not exist. Default is the current working directory (\code{"."}).
  }
}
\value{
A list containing:
\item{permutation_results}{A data frame with the following columns:
  \describe{
    \item{\code{Fixed_Effect}}{Name of the fixed effect or custom contrast tested.}
    \item{\code{Nu_Prime_Eta}}{The observed value of \eqn{\nu' \hat{\eta}}.}
    \item{\code{Permutation_P_Value}}{Permutation p-value for the test of the fixed effect bias.}
  }
}
\item{plot_list}{A list of ggplot2 objects for the permutation histograms.}
}
\examples{
\dontrun{
# Assuming 'result' is your GPvam object

# Perform bias test for all fixed effects
test_results <- bias.test.custom(result)

# Perform bias test including a custom contrast
k_custom <- c(1, -1, 0, 0)  # Contrast between first and second fixed effects
test_results <- bias.test.custom(result, k_vectors = list(k_custom))
}
}
