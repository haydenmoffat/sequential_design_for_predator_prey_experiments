#-----------------------------------------------------------------------------------------------#
# Script Name: total_entropy_utility                                                            #
# Author: Hayden Moffat                                                                         #
# email: hayden.moffat@hdr.qut.edu.au                                                           #
#                                                                                               #
# This R script determines the next design point from the discrete design space which maximises #
# the total entropy utility.                                                                    #                           
#                                                                                               #                                                                                            #
#-----------------------------------------------------------------------------------------------#

Nt <- Nmin:Nmax # possible designs
utility <- matrix(0, length(Nt), 1)  # initialise utility values


for (j in 1:length(Nt)){
  
  log_wsum = matrix(0, K, Nt[j]+1)
  B <- matrix(0, K, Nt[j]+1)
  
  
  for (mod in 1:K){
    
    # Compute the log likelihood for each of the possible y values using the set of particles from model mod
    parameters <- find_parameters(theta[,,mod], Nt[j], time, models[mod])
    y <- 0:Nt[j]
    llh <- log_lik(y, Nt[j], parameters[,1], parameters[,2], models[mod]) # Calculate the log likelihood of observing the datapoint [Nt(j), y]
    
    log_w_hat <- log(W[,mod]) + llh # updated log unnormalised weights of particles if y was the new observation
    log_wsum[mod,] = logsumexp(log_w_hat,1)
    log_W_hat = sweep(log_w_hat, 2, log_wsum[mod,]) # updated log normalised weights
      
    b = llh * exp(log_W_hat) # multiply log likelihood by normalised weights
    b[is.nan(b)] <- 0
    B[mod,] = colSums(b)
    
    rm(parameters, y, llh, log_w_hat, log_W_hat)
  }
  
  log_Z_n = log_Z - logsumexp(log_Z,0) # normalised evidence
  log_p_y = logsumexp(sweep(t(log_wsum), 2, log_Z_n, FUN = '+'), 2) # posterior predictive probabilities
  log_p_y = as.matrix(log_p_y - logsumexp(log_p_y,0), ncol = 1) # normalised posterior predictive probabilities
    
    
  # Determine the expected utility
  utility[j] = exp(log_Z_n)%*% matrix(rowSums(exp(log_wsum) * B), ncol = 1) - (t(log_p_y) %*% exp(log_p_y))
  rm(log_wsum, B, b, log_p_y, log_Z_n)
  
}

# Display the optimal design point
idx = which(utility == max(utility))
cat(paste0('Optimal design point at ', Nt[idx]))
data[i,1] <- Nt[idx]