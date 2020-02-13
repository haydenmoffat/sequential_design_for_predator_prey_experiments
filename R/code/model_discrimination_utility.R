#-----------------------------------------------------------------------------------------------#
# Script Name: model_discrimination_utility                                                     #
# Author: Hayden Moffat                                                                         #
# email: hayden.moffat@hdr.qut.edu.au                                                           #
#                                                                                               #
# This R script determines the next design point from the discrete design space which maximises #
# the model discrimination utility.                                                             #                           
#                                                                                               #                                                                                            #
#-----------------------------------------------------------------------------------------------#


Nt <- Nmin:Nmax # possible designs
utility <- matrix(0, length(Nt), 1) # initialise utility values

for (j in 1:length(Nt)){
  
  log_f <- matrix(0, K, Nt[j]+1)
  log_Z_new <- matrix(0, K, Nt[j]+1)
  
  for (mod in 1:K){
    
    parameters <- find_parameters(theta[,,mod], Nt[j], time, models[mod]) # determine the mu and lambda parameters for all the particles
    y <- 0:Nt[j] # define possible responses
    llh <- log_lik(y, Nt[j], parameters[,1], parameters[,2], models[mod]) # calculate the log likelihood of observing the datapoint [Nt(j), y]
    
    log_w_new <- log(W[,mod]) + llh # determine weights of particles if y was the new observation
    log_f[mod,] <- logsumexp(log_w_new, 1) # ratio of evidences
    log_Z_new[mod,] <- log_Z[mod] + log_f[mod,] # updated evidence
    
  }
  
  log_Z_n = log_Z - logsumexp(log_Z, 0) # normalise current evidence
  log_Z_new_n <- log_Z_new - repmat(logsumexp(log_Z_new, 1),K,1) # normalise updated evidence
  utility[j] = exp(log_Z_n) %*% (rowSums(log_Z_new_n * exp(log_f))) # calculate utility for design point Nt[j]
  rm(parameters, y, llh, log_w_new, log_f, log_Z_new,log_Z_new_n, log_Z_n)
}

# Display the optimal design point
idx = which(utility == max(utility))
cat(paste0('Optimal design point at ', Nt[idx]))
data[i,1] <- Nt[idx]
