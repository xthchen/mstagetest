#' Homogeneity test for multi-stage sampling.
#'
#' The test must be between two populations, but the number of stages does not need to be the same for the two populations.
#' @param n0 Numeric vector/matrix. The sample size of population 0. For one stratum, a numeric vector, with each element being the sample size of the stage. For multiple strata, a matrix, with the rows being the strata and columns the stages.
#' @param n1 Numeric vector/matrix. The sample size of population 1. Same form as n0.
#' @param f0 Numeric vector/matrix. Observed frequencies of population 0. For one stratum, a numeric vector. For multiple strata, a matrix with the rows being the strata and columns the frequencies.
#' @param f1 Numeric vector/matrix. Observed frequencies of population 1. Same form as f0.
#' @return A list with three elements. They are p-values from three different approaches. The first is the approach where there is no frequency pooling between strata/populations. The second with pooling between populations. The third with pooling between strata and populations.
#'
#'@import corpcor
#'
#' @export



mstagetest = function(n0, n1, f0, f1){
  #code is only written for 2 stages for now.
  #n0: the sample size of population 0. For one strata, a numeric vector, 
  #with each element being the sample size of the stage. 
  #For multiple strata, a matrix, with the rows being the strata and columns the stages.
  #n1: likewise, for population 1.
  #f0: Observed frequencies of population 0. For one strata, a numeric vector. 
  #For multiple strata, a matrix with the rows being the strata and columns the frequencies.
  #f1: Likewise, but for population 1.
  #check number of strata nrep.
  if (is.matrix(n0)==TRUE){
    nrep = nrow(n0)
  }else{
    nrep = 1
  }
  #one stratum case
  if (nrep == 1){
    df = length(f0)-1
    #no pooling at all
    tv = final_proportion_cov(n0, f0)$Sigma+final_proportion_cov(n1,f1)$Sigma
    #scalor factor ci
    c0 = 1/n0[2]*(1-1/n0[1])+1/n0[1]
    c1 = 1/n1[2]*(1-1/n1[1])+1/n1[1]
    #pooled frequency fp
    fp = (1/c0*f0+1/c1*f1)/(1/c0+1/c1)
    #pooling between populations
    tvp = final_proportion_cov(n0, fp)$Sigma+final_proportion_cov(n1,fp)$Sigma
    #test statistics
    ts = t(f1-f0)%*%pseudoinverse(tv)%*%t(t(f1-f0))
    tsp = t(f1-f0)%*%pseudoinverse(tvp)%*%t(t(f1-f0))
    #p-values
    pval = 1-pchisq(ts,df)
    pval_p = 1-pchisq(tsp,df)
    return(c(pval, pval_p))
  }else{
    #multiple strata
    df = ncol(f0)-1
    c = ncol(f0)
    
    tv = matrix(0,nrow = c, ncol = c)
    fv = rep(0,c)
    tvp = matrix(0,nrow = c, ncol = c)
    fvp = rep(0,c)
    
    c0 = rep(0,nrep)
    c1 = rep(0,nrep)
    
    pool_u = rep(0,c)
    pool_d = rep(0,c)
    for (i in 1:nrep){
      #no pooling at all
      Sv = final_proportion_cov(n0[i,], f0[i,])$Sigma+final_proportion_cov(n1[i,], f1[i,])$Sigma
      #weighting
      wv = 1/sum(diag(Sv))
      #weighted variance covariance
      tv = tv + wv^2*Sv
      #weighted frequency difference
      fv = fv + wv*(f1[i,]-f0[i,])
      #pooling between populations
      #scalar factor ci
      #c0[i] = 1/n0[i,2]*(1-1/n0[i,1])+1/n0[i,1]
      c0[i] = 1 - prod(1-1/n0[i,])
      #c1[i] = 1/n1[i,2]*(1-1/n1[i,1])+1/n1[i,1]
      c1[i] = 1 - prod(1-1/n1[i,])
      #pooled frequency fp
      fp = (f0[i,]/c0[i]+f1[i,]/c1[i])/(1/c0[i]+1/c1[i])
      #variance covariance
      Svp = final_proportion_cov(n0[i,], fp)$Sigma+final_proportion_cov(n1[i,],fp)$Sigma
      #weighting
      wvp = 1/sum(diag(Svp))
      #weighted variance covariance
      tvp = tvp + wvp^2*Svp
      #weighted frequency difference
      fvp = fvp + wvp*(f1[i,]-f0[i,])
      #pooling between populations and strata
      pool_u = pool_u + f0[i,]/c0[i] + f1[i,]/c1[i]
      pool_d = pool_d + 1/c0[i] + 1/c1[i]
    }
    fpool = pool_u/pool_d
    tvpool = matrix(0,nrow = c, ncol = c)
    fvpool = rep(0,c)
    wvpool0 = 0
    wvpool1 = 0
    for (j in 1:nrep){
      #variance covariance
      Svpool0 = final_proportion_cov(n0[j,], fpool)$Sigma
      Svpool1 = final_proportion_cov(n1[j,],fpool)$Sigma
      #weighting
      wvpool0 = wvpool0 + 1/sum(diag(Svpool0))
      wvpool1 = wvpool1 + 1/sum(diag(Svpool1))
      
    }
    wvpool0 = 1/wvpool0
    wvpool1 = 1/wvpool1
    for(j in 1:nrep){
      #variance covariance
      Svpool0 = final_proportion_cov(n0[j,], fpool)$Sigma
      Svpool1 = final_proportion_cov(n1[j,],fpool)$Sigma
      
      nw0 = wvpool0/sum(diag(Svpool0))
      nw1 = wvpool1/sum(diag(Svpool1))
      
      #weighted variance covariance
      tvpool = tvpool + (Svpool0* nw0^2+ Svpool1* nw1^2)
      #weighted frequency difference
      fvpool = fvpool + (f1[j,]*nw1 - f0[j,]*nw0)
    }
    #no pooling
    ts = t(fv)%*%pseudoinverse(tv)%*%t(t(fv))
    #pooled between populations
    tsp = t(fvp)%*%pseudoinverse(tvp)%*%t(t(fvp))
    #pooled between populations and strata
    tspool = t(fvpool)%*%pseudoinverse(tvpool)%*%t(t(fvpool))
    
    pval = 1-pchisq(ts,df)
    pvalp = 1-pchisq(tsp,df)
    pvalpool = 1-pchisq(tspool,df)
    return(c(pval, pvalp, pvalpool))
  }
  
}

