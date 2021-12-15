# Author: Roberto Zuniga
# Date: December 14, 2021
# Description: This script partially replicates Blunder, Gowrisankaran and Blundell
# (2020). Specifically, It implements the random coefficients estimator.

library(tidyr)
library(ggplot2)

# Set working directories
current_dir <- "./IO_part_3"
setwd(current_dir)

outputs_dir <- paste0(current_dir,"Outputs")
bellman_data_file <- paste0(current_dir,"Data/data_for_bellman_computations.csv")
panel_data_file <- paste0(current_dir,"Data/analysis_data.csv")

source("Bellman_functions.R")

# Read Bellman data
data_bellman <- read.csv(bellman_data_file)
data <- subset(data_bellman,naics_recode == 1)
data$index <- 1:nrow(data)

# Read analysis data for lilkelihood computations
panel_data <- read.csv(panel_data_file)
panel_data <- subset(panel_data,naics_recode == 1)

# Define the working parameters
params <- list()
params$data <- data
params$panel_data <- panel_data 
params$N1 <- 161  # number of states per each omega1
params$n_omega1 <- length(unique((data$omega1))) # number of omega1 states
params$DAVgrid <- (0:19)/19*9.5  # DAV grid
params$n_DAV <- length(params$DAVgrid) # number of points in DAV grid
params$DAV_deprate <- 0.1  # DAV depreciation rate
params$beta <- 0.95^(1/4)  # discount factor
params$gamma <- 0.57721566490153*0 # Euler's gamma constant

params$tol <- 1E-6 # stopping condition for the Bellman equation solution
params$nregoutcomes <- 80 # number of regulatory outcomes
params$ntransitions <- 3  # number of states transitions from omega to omega_tilde
################################################################################
# Create structures for fast computation of Bellman equation
params <- AddParams(params)
# Add mapping from omega1 to (region,naics,gravity)
params <- GetOmega1Mapping(params)

################################################################################
# Solve the Bellman equation
coeff_BGL <- c(2.872, -0.049, -0.077, -5.980, -0.065) # BGL
resultsBGL <- Bellman(coeff_BGL, params)
################################################################################
# Compute transition probabilities

# Add helper matrices
params <- AddTransitionParams(params)

# Get Transition matrix
omega1 <- 1
P_omega1_BGL <- GetTransitionProb(params,omega1,resultsBGL$Investprob)

# check results
ExpectedResult1 <- paste0(current_dir,"Data/transitionprob_part3_problem1_DAVgrid2_thetaBGL.csv")
CompareTransitionProb(P_omega1_BGL,ExpectedResult1,params)

# Question 1
FileQ1 <- "Output/transProb_DAV4_P3Q1.csv"
Question_P3Q1(P_omega1_BGL,FileQ1,params)

################################################################################
# Computation of the steady state distribution
piSS_omega1_BGL <- SolveSS(P_omega1_BGL,params)

# check results
ExpectedResult2 <- paste0(current_dir,"Data/steadystate_part3_problem2_DAVgrid2_thetaBGL.csv")
CompareSS(piSS_omega1_BGL,ExpectedResult2,params)

# Question 2
FileQ2 <- "Output/steadyState_DAV4_P3Q2.csv"
Question_P3Q2(piSS_omega1_BGL,FileQ2,params)

################################################################################
# Computation of the model moment inputs across parameter grid values

# Read parameter grid
GridFile <- paste0(current_dir,"Data/parameter_grid_assignment.csv")
parameterGrid <- as.matrix(read.csv(GridFile))

# check results
moments_param2 <- ComputeMoments(parameterGrid[2,], params)
moments_param2 <- moments_param2$moments
ExpectedResult3 <- paste0(current_dir,"Data/modelmoments_part3_problem3_param2.csv")
CompareMoments(moments_param2,ExpectedResult3,params)

# Question 3
FileQ3 <- "Output/moment_mean_sd_P3Q3.csv"
Question_P3Q3(parameterGrid,params,FileQ3)

################################################################################
# Computation of the data moment inputs

dataMoments <- ComputeDataMoments(params)
# Save data moments
#write.csv(data.frame(dataMoments),"Output/all_data_moments.csv",row.names = FALSE)

# Check results
ExpectedResult4 <- paste0(current_dir,"Data/datamoments_part3_problem4_param2.csv")
CompareDataMoments(dataMoments,ExpectedResult4,params)

# Question 4
FileQ4 <- "Output/datamoment_mean_sd_P3Q4.csv"
Question_P3Q4(dataMoments,params,FileQ4)

################################################################################
# Generate moments file
File1_Q5 <- "Output/moment_data_class_assignment_use.csv"
File2_Q5 <- "Output/moments_full_info.csv"
GenerateMomentFiles(parameterGrid,dataMoments,params,File1_Q5,File2_Q5)
