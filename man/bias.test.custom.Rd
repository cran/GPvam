\name{bias.test.custom}
\alias{bias.test.custom}
\title{
Permutation Tests for Fixed Effects Bias Assessment
}
\description{
Performs permutation tests on fixed effects within a linear mixed model to assess the significance of each effect and potential biases. The function allows for both standard basis vectors and custom vectors to define the effects being tested. See Karl and Zimmerman (2021) <doi:10.1016/j.jspi.2020.06.004>.


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
    An object containing model results, including the fixed effects matrix (`X`), random effects design matrix (`Z`), inverse variance matrix (`vinv`), and estimated random effects (`eta.hat`). The object must contain the following components:
    \describe{
      \item{\code{X}}{Fixed effects matrix.}
      \item{\code{Z}}{Random effects design matrix.}
      \item{\code{vinv}}{Inverse variance matrix.}
      \item{\code{eta.hat}}{Estimated random effects.}
    }
  }
  
  \item{k_vectors}{
    (optional) A list of numeric vectors specifying custom `k` vectors for combined fixed effects. Each vector should be the same length as the number of fixed effects in the model. If `NULL`, the function generates standard basis vectors (one-hot vectors) to perform permutation tests for each fixed effect individually.
  }
  
  \item{n_perms}{
    (optional) The number of permutations to perform for each `k` vector. A higher number of permutations increases the accuracy of the percentile and p-value estimates but also increases computation time. Default is \code{1e6}.
  }
  
  \item{save_plots}{
    (optional) Logical flag indicating whether to save the permutation test histograms as PNG files. If set to \code{TRUE}, histograms for each `k` vector will be saved in the specified \code{output_folder}. Default is \code{FALSE}.
  }
  
  \item{output_folder}{
    (optional) Directory path where plots and results will be saved if \code{save_plots} is \code{TRUE}. The directory will be created if it does not exist. Default is the current working directory (\code{"."}).
  }
}
