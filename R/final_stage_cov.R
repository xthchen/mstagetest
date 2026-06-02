#' Computes the sample size covariance matrix given stage sample size and observed frequency. This only works for one strata. For multiple stratum, this function is applied recursively, assuming independence between stratum.
#'
#' @param n_list Integer, sample size at each stage. Each Each element denotes the sample size of a stage. The number of elements denote the number of stages.
#' @param p Numeric, vector of observed frequencies, should add up to 1.
#'
#' @return A list. The first element is a numeric matrix, which is the sample size variance covariance matrix. The second element is numeric, its the equivalent sample size across stages. The third element is a numeric matrix, which is the the unmodified frequency variance covariance.
#'
#' @export
#' 


final_stage_cov = function(n_list, p) {
  
  S = length(n_list)
  #number of sampling stages
  if (S < 1) stop("n_list must have at least one stage")
  if (abs(sum(p) - 1) > 1e-12) stop("p must sum to 1")
  K <- length(p)
  #number of categories
  
  #finding the effective sample size a_S
  a_prev = n_list[1]
  #first stage sample size
  if (S == 1) {
    #if only one sampling step
    a_S = a_prev
  } else {
    for (t in 2:S) {
      n_t = n_list[t]
      #sample size of t-th sampling step
      n_prev = n_list[t-1]
      #sampling size of previous sampling step
      a_t = n_t + (n_t^2 - n_t) * (a_prev / (n_prev^2))
      #see section 3.1
      a_prev = a_t
      #n1(n1+n0-1)/n0
    }
    a_S = a_prev
  }
  
  # construct V = diag(p) - p p^T
  V = diag(p, nrow = K, ncol = K) - outer(p, p)
  #diag(f)-fft
  # final covariance
  Sigma <- a_S * V
  return(list(Sigma = Sigma, a_S = a_S, V = V))
  #Sigma the sample size variance covariance, a_S the effective population size, 
  #V the unmodified frequency variance covariance.
}
