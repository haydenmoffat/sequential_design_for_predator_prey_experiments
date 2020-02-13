#-----------------------------------------------------------------------------------------------#
# Script Name: SMC_algorithm                                                                    #
# Author: Hayden Moffat                                                                         #
# email: hayden.moffat@hdr.qut.edu.au                                                           #
#                                                                                               #
# This R script conducts sequential Bayesian experimental design for functional response        #
# experiments. The entire sequential experimental design algorithm can be run from this script. #                           
#                                                                                               #                                                                                            #
#-----------------------------------------------------------------------------------------------#


##############################################################
# WORKSPACE SETUP
##############################################################

# Clear workspace
rm(list = ls())

# Set seed for reproducibility
 set.seed(100)

# Load required packages
library(dplyr)
library(deSolve)
library(mvtnorm)
library(pracma)
library(reshape)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggthemes)
library(latex2exp)

# Load required functions
source('generate_data.R')
source('additional_functions.R')




##############################################################
# PROBLEM AND SMC SET UP (INPUT REQUIRED)
##############################################################

# Specify method of determining next design point
R <- 0

# Methods (R):
# 0 - selects design point randomly and generates data
# 1 - selects design point that maximises utility for parameter estimation and generates data
# 2 - selects design point that maximises utility for model discrimination and generates data
# 3 - selects design point that maximises total entropy utility and generates data
# 4 - selects both design and observation points from an example experiment
# 5 - selects design points used in an example experiment and generates data

# Specify models to be considered
models <- 1:4 %>% matrix(nrow = 1)
modeltype <- list("Beta Binomial type 2 functional response model",
                  "Beta Binomial type 3 functional response model",
                  "Binomial type 2 functional response model",
                  "Binomial type 3 functional response model")

# Load example dataset for methods that require it
if (R %in% 4:5){
  dataset <- read.table("papanikolau_data.txt", sep = " ", col.names = c("N", "n")) %>% as.matrix
}

#  Select "true values" of parameters for methods that require it
if (R %in% c(0:3,5)){
  a <- 0.5
  Th <- 0.7
  lambda <- 0.5
  parameter_set <- matrix(c(a, Th, lambda), nrow = 1) %>% log
  M_true <- 1
}

# Assign values to parameters (Input required for these values)
I <- 25 # number of data points
K <- length(models) # number of models
N <- 500 # number of particles for SMC
E <- N/2 # threshold for ESS for SMC
Nmin <- 1 # min value in the discrete design space
Nmax <- 300 # max value in the discrete design space
time <- 24 # length of time that predator has access to prey in hours
tol <- 2 # tolerance for ESS which indicates when the move step might fail (see move step script for more detail)






##############################################################
# CONDUCTING THE SMC (NO INPUT REQUIRED)
##############################################################

# Initialise other required quantities (Input not needed for these values)
data <- matrix(0, I, 2) # Set up data
log_Z = matrix(0, 1, K) # initialise log estimate of evidence
loglik = matrix(0, N, K) # initialise vector of log likelihood
logpri = matrix(0, N, K) # initialise vector of log prior
px = matrix(0, N, K) # initialise vector of the log posterior
w = matrix(0, N, K) # initialise vector of unnormalised weights of particles
ESS = matrix(0, 1, K) # initialise vector of effective sample size of each model

# Draw theta from prior distribution, set intial weighting, set covariance matrices
source('prior_sampling.R')
W <- rep(1/N, N * K) %>% matrix(nrow = N)

all_cov_matrices <- list()
for(M in 1:K){
  all_cov_matrices[[M]] = cov(theta[,,M])
}

## RUN SMC ALGORITHM

for (i in 1:I){
  
  # Print iteration number to console
  cat(paste("   **** Iteration number", i, "****   "), sep ="\n")
  cat(paste("--------------------------------------   "), sep ="\n")

  if (R %in% c(0, 1, 2, 3, 5)){
    
    if  (R == 0){
      # Select design points randomly
      data[i,1]= sample(Nmin:Nmax, 1, replace = T)
    }
    
    
    else if (R == 1){
      source('parameter_estimation_utility.R') # select design points by maximising the parameter estimation utility
      source('utility_plot.R') # plot utilty curve
      cat(" ", sep ="\n")  
    }
    
    
    else if (R == 2){
      source('model_discrimination_utility.R') # select design points by maximising the model discrimination utility
      source('utility_plot.R') # plot utility curve
      cat(" ", sep ="\n")  
    }
    
    
    else if (R == 3){
      source('total_entropy_utility.R') # select design points by maximising the total entropy utility
      source('utility_plot.R') # plot utility curve
      cat(" ", sep ="\n")  
    }
    
    
    else {
      data[i,1] <- dataset[i,1] # take design point from example experiment
    }
    
     # Generate observation from the design point
    data[i,2] <- generate_data(data[i,1], parameter_set, observation_time = time, M_true)[[1]]
    
  }
  
  else if(R == 4){
    data[i,] <- dataset[i,] # take all design point and observations from example dataset
  }

  data_subset <- data[1:i,, drop = FALSE] # current data available
  
  for (M in 1:K){
    
    #  Re-weight
    w[,M] = W[,M] * exp(loglikelihood_Hollings(exp(theta[,,M]), data[i,,drop = FALSE], observation_time = time, models[M]))
    log_Z[M] = log_Z[M] + log(sum(w[,M]))
                          
    # Normalise weights
    W[,M] = w[,M]/sum(w[,M])
                          
    # Compute ESS at time t
    ESS[M] = 1/ sum(W[,M]^2)
    cat(paste("Model", models[M], "ESS:", floor(ESS[M])), sep="\n")
    
    if (ESS[M] < E){
      source('resample_then_move.R')
    }
  }
  
cat(" ", sep ="\n")  

}




##############################################################
# FINAL RESAMPLING AND RESULTS
##############################################################


## Conduct a final resampling and move step to ensure all weights of particles are equal
for (M in 1:K){
  cat(paste("Model", models[M],":"), sep ="\n")
  
  if (ESS[M] != N){
    source('resample_then_move.R')
  }
  
  # Display log evidence
  cat(" ", sep ="\n")  
  cat(paste('log_Estimate_of_Evidence =', log_Z[M]), sep="\n") 
  
  # Calculate and display log Bayesian D-posterior precision for the parameters for each model
  covariance = cov(exp(theta[,,M]))
  if (models[M] %in% 3:4){covariance = covariance[1:2, 1:2]}
  log_D_post_prec = -log(det(covariance));
  cat(paste('log D-posterior precision =', log_D_post_prec), sep="\n")
  
  # Plot marginal posterior distributions
  source('SMC_plot.R')
  cat(" ", sep ="\n")  
}


# Calculate Posterior Model Probabilities
posterior_model_prob = exp(log_Z - logsumexp(log_Z,0))
for (M in 1:K){
  cat(paste(modeltype[[models[M]]], "probability:", posterior_model_prob[M]), sep="\n")
}
