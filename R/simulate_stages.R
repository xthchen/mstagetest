#' Given starting frequencies, simulate frequencies after some stages using multinomial sampling.
#'
#' @param seq_counts Integer, sample size at each stage. Each Each element denotes the sample size of a stage. The number of elements denote the number of stages.
#' @param p1 Numeric, vector of starting frequencies. Should add up to 1.
#' @param seed Numeric, simulate using the given seed. Will not have a fix seed with no input.
#'
#' @return A list. First element is the simulated frequency. Second element the sample size of the last stage.

#' @export



simulate_stages = function(seq_counts, p1, seed = NULL) {
  #use to simulate the frequencies after some stages
  #seq_counts: Numeric, sample size at each stage. Each element denotes the sample size of a stage.
  #p1: Numeric, starting frequencies, add up to one.
  #seed: whether to fix a seed.
  if (!is.null(seed)) set.seed(seed)
  #set seed
  S = length(seq_counts)
  #number of stages
  K = length(p1)
  #number of categories(unused variable)
  
  p = p1
  #starting category probability
  for (t in 1:S) {
    x = as.vector(rmultinom(1, seq_counts[t], p))
    #multinomial sampling according to first sample size
    
    if (t < S) {
      p = x / seq_counts[t]   # update probabilities
      #for second sampling step
    }
  }
  return(list(x,seq_counts[t]))
  #returns probability of final stage and sample size of last stage
}
