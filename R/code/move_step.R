#-----------------------------------------------------------------------------------------------#
# Script Name: move_step                                                                        #
# Author: Hayden Moffat                                                                         #
# email: hayden.moffat@hdr.qut.edu.au                                                           #
#                                                                                               #
# This R script contains code to conduct a MCMC move step on the particles after the resampling #
# step.                                                                                         #                           
#                                                                                               #                                                                                            #
#-----------------------------------------------------------------------------------------------#


print("MOVE STEP")

# If ESS is very low -> use previous covariance matrix
# Otherwise -> compute parameters of MCMC kernel qt, and override previous covariance matrix

if(ESS[M] > tol){
  c <- cov(theta[,,M])
  all_cov_matrices[[M]] <- c
}

cov_matrix <- all_cov_matrices[[M]]


# Move each particle theta with MCMC kernel for 1 iteration
n_moves = rep(0, N*K) %>% matrix(nrow = N, ncol = K) # variable that describes the number of particles that have moved after 1 iteration
for (l in 1:N){
  u <- runif(1)
  x <- theta[l,,M] %>% matrix(nrow = 1)
  logpri[l,M] = log_prior(x, models[M]) # determine log prior of initial theta values
  loglik[l,M] = loglikelihood_Hollings(exp(x), data_subset, time, models[M]) # determine log likelihood of initial theta values
  
  px[l,M] = loglik[l,M]+logpri[l,M]

  xs <- rmvnorm(1, x, cov_matrix) # new sample xs based on existing x from proposal pdf.
  
  # Determine log prior and log likelihood of proposed theta values
  logprior_xs <- log_prior(xs, models[M])
  loglik_xs <- loglikelihood_Hollings(exp(xs), data_subset, time, models[M])
  pxs <- loglik_xs +logprior_xs

  
  if (is.nan(pxs) == 0){
    if (u<min(1,exp(pxs-px[l,M]))){
      theta[l,,M] <- xs
      logpri[l,M] <- logprior_xs
      loglik[l,M] <- loglik_xs
      px[l,M] <- pxs
      n_moves[l,M] = n_moves[l,M]+1
    }
  }
 
}


# Determine acceptance probability and an appropriate number of iterations
# of the MCMC
prob = sum(n_moves[,M])/N
Rt = ceiling(log(0.01)/log(1-prob))

# Displays the acceptance probability and the number of iterations
# (Rt) for each move step
print(paste("Acceptance probability =", prob, ", Rt = ", Rt))

# Move particles the rest of the iterations
for (l in 1:N){
  for (q in 2:Rt){
    u <- runif(1)
    x <- theta[l,,M] %>% matrix(nrow = 1)
    xs <- rmvnorm(1, x, cov_matrix) # new sample xs based on existing x from proposal pdf.
    
    # Determine log prior and log likelihood of proposed theta values
    logprior_xs <- log_prior(xs, models[M])
    loglik_xs <- loglikelihood_Hollings(exp(xs), data_subset, time, models[M]);
    pxs <- loglik_xs +logprior_xs;
    
    if (u<min(1,exp(pxs-px[l,M]))){
      theta[l,,M] <- xs
      logpri[l,M] <- logprior_xs
      loglik[l,M] <- loglik_xs
      px[l,M] <- pxs
      n_moves[l,M] = n_moves[l,M]+1

    }
  }
}


