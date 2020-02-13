#-----------------------------------------------------------------------------------------------#
# Function Name: generate_data                                                                  #
# Author: Hayden Moffat                                                                         #
# email: hayden.moffat@hdr.qut.edu.au                                                           #
#                                                                                               #
# This function generates random observations at specified design points from either the       #
# beta-binomial distribution or the binomial distribution (depending on the model of interest). #        
#                                                                                               #
# Inputs:                                                                                       #
# design - design at which we will generate the data                                            #
# theta -  true model parameters                                                                #
# observation_time - length of time that predator has access to prey in hours                   #
# model - model of interest                                                                     #
#                                                                                               #
# Ouputs:                                                                                       #
# data - the response generated from the true model                                             #
#                                                                                               #                                                                                            #
#-----------------------------------------------------------------------------------------------#


generate_data <- function(design, theta, observation_time, model){
  
  if (model %in% 1:2){
    
    # Beta-binomial distribution
    v <-  solve_ode_Hollings(exp(theta), N0= design, observation_time = observation_time, model = model) # solve ode
    
    # Determine parameters
    design_matrix <- matrix(rep(design, times = nrow(theta)), byrow = T, nrow = nrow(theta))
    mu <- as.matrix((design_matrix - v)/design_matrix)
    alpha <- as.matrix(apply(mu, 2, '/', exp(theta[,3])))
    beta <- as.matrix(apply((1 - mu), 2, '/', exp(theta[,3])))

    # Generate data
    final_data <- list()
    for (i in 1:nrow(theta)){
      data <- rep(NA, length(design))
      for (j in 1:length(design)){
        data[j] <- rbinom(1, design[j], rbeta(1, alpha[i,j], beta[i,j]))
      }
      final_data[[i]] <- data
      rm(data)
    }
  }
  
  else if (model %in% 3:4){
    # Binomial distribution
    v <- solve_ode_Hollings(exp(theta), N0= design, observation_time = observation_time, model= model) # solve ode 
    
    design_matrix <- matrix(rep(design, times = nrow(theta)), byrow = T, nrow = nrow(theta))
    mu <- (design_matrix - v)/design_matrix # determine parameter
    
    # Generate data
    final_data <- list()
    for (i in 1:nrow(theta)){
      data <- rep(NA, length(design))
      for (j in 1:length(design)){
        data[j] <- rbinom(1, design[j], mu[i,j])
      }
      final_data[[i]] <- data
      rm(data)
    }
  }
  
  else{
    stop("Non valid input for model")
  }
  
  return(final_data)
}

