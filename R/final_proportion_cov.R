#' Computes the frequency covariance matrix given stage sample size and observed frequency. This only works for one strata. For multiple stratum, this function is applied recursively, assuming independence between stratum.
#'
#' @param n_list Integer, sample size at each stage. Each Each element denotes the sample size of a stage. The number of elements denote the number of stages.
#' @param p Numeric, vector of observed frequencies, should add up to 1.
#'
#' @return A list. The first element is a numeric matrix, which is the frequency variance covariance matrix. The second element is numeric, its the equivalent sample size across stages. The third element is a numeric matrix, which is the the unmodified frequency variance covariance.
#'
#' @export
#' 

final_proportion_cov = function(n_list, p) {
  # frequency covariance matrix for multinomial
  # proportion instead of counts
  # only for one strata
  # n_list: numeric, integer vector of stage sizes (length S)
  # p: numeric, initial probability vector p^(1) (length K)
  
  S = length(n_list)
  #number of stages
  nS = n_list[S]
  #last stage samsize
  
  # closed form for a_S
  prod_term <- prod(1 - 1 / n_list[1:S])
  a_S <- nS^2 * (1 - prod_term)
  
  # V = diag(p) - p %*% t(p)
  V <- diag(p) - outer(p, p)
  
  # covariance for proportions
  Sigma_p <- (a_S / nS^2) * V
  
  return(list(Sigma = Sigma_p, a_S = a_S, V = V))
  #sigma is the variance covariance of frequencies
}
