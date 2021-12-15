#' Solve the Bellman equation
#' @param coeff The coefficient at which the value function will be estimated
#' It is assumed that coeff = (theta_X, theta_I,theta_V, theta_F, theta_H) 
#' @param params A list of estimation parameters
#' @param V0 An (optional) initial value for the value function
Bellman <- function(coeff, params,V0){
  
  # Construct matrices and indexes that are necessary for solving the Bellman equation 
  if (!exists("BellmanParams",where = params)){
    params$BellmanParams <- PackBellmanParams(params)
  }

  # Run the loop to find the fixed point of the Bellman equation
  if(!missing(V0)) {
    return(BellmanLoop(coeff, params,V0))
  }
  else {
    return(BellmanLoop(coeff, params))
  }
  
}

#' Preliminary computations for solving the Bellman equation. This function
#' creates matrices and mappings between Omega states and matrix indexes
#' These matrices (and indeces) allow for faster comptutations in the Bellman loop
#' @param params A list of estimation parameters
PackBellmanParams <- function(params){
  
  # Estimation parameters
  n_omega1 <- params$n_omega1          # number of omega1 states
  N1 <- params$N1                      # number of states for each omega1
  nstates <- n_omega1*N1               # Number of states
  n_regoutcomes <- params$nregoutcomes # Number of regulatory outcomes (80)
  n_transitions <- params$ntransitions # violator transtions (3)
  dep <- params$DAV_deprate            # depreciation rate
  beta <- params$beta                  # discount factor
  gamma <- params$gamma                # euler's constant
  
  # Initialize matrices 
  Vtilde <- matrix(0,nstates,1)
  NewV   <- matrix(0,nstates,1)
  OldV   <- matrix(0,nstates,1)
  Investprob <- matrix(0,nstates,1)
  
  # The logic of the BellmanLoop function is to work with one-dimensional vectors
  # representing V and Vtilde. This allows for faster computation. But also requires
  # defining a mapping from the multidimensional Omega/Omega_tilde states to a 
  # one-dimensional index. We will use indexes to access V and Vtilde. Transitions from 
  # Omega to Omega_tilde, and from Omega_tilde to Omega, involve changing particular
  # states, such as DAV, lag investment, HPV status, etc. In these cases, we 
  # recover the states from indexes, change the states, and then compute the new 
  # indexes corresponding to these new states.
  
  # Omega (or omega tilde) states are indexed from 1 to nstates
  index <- c(1:nstates)
  # Each index corresponds to a unique Omega/Omega_tilde state
  # states are (omega1,lagInv1,lagInv2,orderViol,DAV) for NewV or OldV
  # or         (omega1,lagInv1,currentV,orderViol,DAV) for vtilde

  # Here, IndexToValues maps an index to the corresponding 5-dimensional state
  states <- t(sapply(index,IndexToValues)) 
  
  # Get updated states when the plant invests (transition from omega_tilde to omega)
  statesI <- states
  statesI[,2] <- 1 # lag1 of investment
  statesI[,3] <- states[,2] # lag2 of investment
  
  # Get updated states when the plant does not invest
  statesNI <- states
  statesNI[,2] <- 0 # lag1 of investment
  statesNI[,3] <- states[,2] # lag2 of investment
  
  # Get updated index when the plant transitions into compliance
  statesComp <- states
  statesComp[,2:4]<-0
  statesComp[,5]<-1
  # Here, ValuesToIndex maps a 5-dimensional state vector into its corresponding
  # index
  indexComp <- apply(statesComp,1,ValuesToIndex)
  
  # Get states in compliance
  Comp <- states[,4] == 0
  NComp <- states[,4] > 0
  
  # Get DAV values
  DAVgrid <- params$DAVgrid  # retrieve DAV grid
  DAV <- DAVgrid[states[,5]] # compute DAV
  npointsDAV  <- length(DAVgrid) # number of points in the grid
  maxpointDAV <- DAVgrid[npointsDAV] # maximum DAV value
  widthDAV <- DAVgrid[2]-DAVgrid[1]  # grid width
  
  # Compute updated DAV after depreciation and current violation
  newDAV <- (1-dep)*DAV + states[,3]
  
  # Interpolate newDAV using the grid: find interpolation DAVgrid indexes and 
  # weights
  
  # DAVgrid indexes
  i_below <- pmin(1 + floor(newDAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  i_above <- pmin(2 + floor(newDAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  # weights
  w_below <- pmax((DAVgrid[i_above]-newDAV)/widthDAV,0)
  w_above <- 1-w_below
  
  # Find updated indexes (for DAV above and below) when the plant invests
  statesI_below <- statesI
  statesI_below[,5] <- i_below
  indexI_below <- apply(statesI_below,1,ValuesToIndex)
  statesI_above <- statesI
  statesI_above[,5] <- i_above
  indexI_above <- apply(statesI_above,1,ValuesToIndex)
  
  # Find updated indexes (for DAV above and below) when the plant does not invest
  statesNI_below <- statesNI
  statesNI_below[,5] <- i_below
  indexNI_below <- apply(statesNI_below,1,ValuesToIndex)
  statesNI_above <- statesNI
  statesNI_above[,5] <- i_above
  indexNI_above <- apply(statesNI_above,1,ValuesToIndex)
  
  # Recover inspections, fines and violations in each regulatory outcome
  # Also recover probability transitions
  
  # First, define names of columns for each variable-regulatory outcome case
  # Initialize as empty lists
  ins_list <- c()
  vio_list <- c()
  fine_list <- c()
  prob_list <- c()
  # Loop through regulatory outcomes to construct column names
  for (j_ro in 1:n_regoutcomes){
    
    ins_name <- paste0("inspection",j_ro)
    ins_list <- c(ins_list,ins_name)
    
    vio_name <- paste0("violation",j_ro)
    vio_list <- c(vio_list,vio_name)
    
    fine_name <- paste0("fine",j_ro)
    fine_list <- c(fine_list,fine_name)
    
    prob_name <- paste0("prob",j_ro)
    prob_list <- c(prob_list,prob_name)
  }
  
  # Map (regularotory outcome, transition) pairs into one-dimensional indexes
  jRo <- rep(NA, n_regoutcomes*n_transitions)
  jTr <- rep(NA, n_regoutcomes*n_transitions)
  for (j_ro in 1:n_regoutcomes){
    for (j_tr in 1:n_transitions){
      jRo[(j_tr-1)*n_regoutcomes+j_ro] <- j_ro
      jTr[(j_tr-1)*n_regoutcomes+j_ro] <- j_tr
    }
  }
  
  # Get inspections, violations, fines, and new HPV status from the data
  # Here, each "case" is a regulatory outcome - sate transition pair
  insCase  <- params$data[index,ins_list[jRo]]
  vioCase  <- params$data[index,vio_list[jRo]]
  fineCase <- params$data[index,fine_list[jRo]]
  ordVioCase  <- matrix(rep(jTr,nstates),nstates,n_regoutcomes*n_transitions,byrow = T)
  HPVCase  <- 1*(ordVioCase == 3) 
  
  # Compute the probability of each case
  prob_ro <- params$data[index,prob_list[jRo]]
  transition_name <- paste0("transition",jTr-1,"_",jRo)
  prob_tr  <- params$data[index,transition_name]
  CaseProb <- prob_ro*prob_tr
  
  # Reshape as a vector
  CaseProb <- c(as.matrix(CaseProb))
  
  # Transition from omega to omega_tilde: update violation and HPV status
  # Here we define possible transitions to the next omega_tilde state
  # each entry is a pair (current_violation,new_violator_status)
  # there are 5 possible transitions
  transitions <- matrix(c(0,0,0,1,1,0,1,2,1,2),5,2)
  # It is helpful to reshape CaseProb, vioCase and ordVioCase
  CaseProbMatrix  <- matrix(CaseProb,nstates,n_regoutcomes*n_transitions)
  vioCaseMatrix     <- matrix(c(as.matrix(vioCase)),nstates,n_regoutcomes*n_transitions)
  ordVioCaseMatrix  <- matrix(c(as.matrix(ordVioCase)),nstates,n_regoutcomes*n_transitions)
  # Get transition probabilities and indexes
  indexTransition <- matrix(0,nstates*5,1)
  ProbTransition <- matrix(0,nstates*5,1)
  for (i in 1:5){
    # Get the index of the new omegatilde state after the transitions
    indexTransition[(nstates*(i-1)+1):(nstates*i),1] <- getIndexOfTransition(states,
                                                  transitions[i,1],transitions[i,2])
    # Get the probability of this transition
    ProbTransition[(nstates*(i-1)+1):(nstates*i),1] <- CollapseProbByTransitions(CaseProbMatrix,
                                                        vioCaseMatrix, ordVioCaseMatrix,transitions[i,1],
                                                        transitions[i,2])
  }
  
  # Pack Bellman parameters
  paramsBellman <- list()
  paramsBellman$NewV              <- NewV          
  paramsBellman$OldV              <- OldV          
  paramsBellman$w_below           <- w_below       
  paramsBellman$w_above           <- w_above       
  paramsBellman$indexI_below      <- indexI_below  
  paramsBellman$indexNI_below     <- indexNI_below 
  paramsBellman$indexI_above      <- indexI_above  
  paramsBellman$indexNI_above     <- indexNI_above 
  paramsBellman$states            <- states        
  paramsBellman$CaseProb          <- CaseProb      
  paramsBellman$insCase           <- insCase       
  paramsBellman$vioCase           <- vioCase       
  paramsBellman$fineCase          <- fineCase      
  paramsBellman$HPVCase           <- HPVCase
  paramsBellman$indexComp         <- indexComp
  paramsBellman$Comp              <- Comp
  paramsBellman$NComp             <- NComp
  paramsBellman$indexTransition  <- indexTransition
  paramsBellman$ProbTransition   <- ProbTransition
  
  return(paramsBellman)
}

#' Solve the Bellman equation
#' Assumes that BellmanParams are already stored in the params argument
#' @param coeff The coefficient at which the value function will be estimated
#' It is assumed that coeff = (theta_X, theta_I,theta_V, theta_F, theta_H) 
#' @param params A list of estimation parameters
#' @param V0 An (optional) initial value for the value function
BellmanLoop <- function(coeff, params,V0){
  
  # Retrieve Bellman parameters
  paramsBellman <- params$BellmanParams
  
  # Unpack Bellman parameters
  NewV <- paramsBellman$NewV
  OldV <- paramsBellman$OldV
  w_below <- paramsBellman$w_below
  w_above <- paramsBellman$w_above
  indexI_below <- paramsBellman$indexI_below
  indexNI_below <- paramsBellman$indexNI_below
  indexI_above <- paramsBellman$indexI_above
  indexNI_above <- paramsBellman$indexNI_above
  states <- paramsBellman$states
  CaseProb <- paramsBellman$CaseProb
  insCase  <- paramsBellman$insCase
  vioCase  <- paramsBellman$vioCase
  fineCase <- paramsBellman$fineCase
  HPVCase  <- paramsBellman$HPVCase
  indexComp <- paramsBellman$indexComp
  Comp <- paramsBellman$Comp 
  NComp <- paramsBellman$NComp 
  indexTransition <- paramsBellman$indexTransition
  ProbTransition  <- paramsBellman$ProbTransition
  beta <- params$beta                  # discount factor
  gamma <- params$gamma                # euler's constant
  n_omega1 <- params$n_omega1          # number of omega1 states
  N1 <- params$N1                      # number of states for each omega1
  nstates <- n_omega1*N1               # Number of states
  n_regoutcomes <- params$nregoutcomes # Number of regulatory outcomes (80)
  n_transitions <- params$ntransitions # violator transitions (3)
  
  if(!missing(V0)) {
    OldV <- V0
  } 
  
  # Unpack coefficient (theta)
  theta_X <- coeff[1]
  theta_I <- coeff[2]
  theta_V <- coeff[3]
  theta_F <- coeff[4]
  theta_H <- coeff[5]
  
  # Per period utility
  U <- theta_I*insCase + theta_V*vioCase + theta_F*fineCase +
    theta_H*HPVCase
  U <- c(as.matrix(U))
  U <- rowSums(matrix(U*CaseProb,nstates,n_regoutcomes*n_transitions))
  
  # Loop until NewV converges to the fixed point
  iconv = 0 # counts the number of iterations
  norm_Vdiff <- 1
  while (norm_Vdiff>=params$tol){
    iconv = iconv + 1
    
    # 1) Update Vtilde and investment probabilities 

    # Solve OldV if investment
    OldV_I <- w_below*OldV[indexI_below] + w_above*OldV[indexI_above]
    
    # Solve OldV if no investment
    OldV_NI <- w_below*OldV[indexNI_below] + w_above*OldV[indexNI_above]
    
    # Use logit value to find Vtilde
    V_C  <- gamma + beta*OldV[indexComp] # value if in compliance
    V_NC <- gamma + log(exp(beta*OldV_NI)+exp(-theta_X+beta*OldV_I)) # not in compliance
    
    # Update Vtilde
    Vtilde <- Comp*V_C + NComp*V_NC 
  
    # 2) Update NewV 
    EVtilde <- Vtilde[indexTransition]*ProbTransition
    EVtilde <- matrix(EVtilde,nstates,5)
    NewV <- U + rowSums(EVtilde)
    
    # 3) Calculate norm of diff between OldV and NewV 
    norm_Vdiff <- max(abs(OldV-NewV))
    OldV <- NewV
    #print(paste0("Iteration # ",iconv,", norm = ",norm_Vdiff))
  }
  
  # 4) Pack and return main results 
  results <- list()
  
  # Save prob. of investment
  Investprob <- ((states[,4] > 0))*
    exp(-theta_X+beta*OldV_I)/(exp(beta*OldV_NI)+exp(-theta_X+beta*OldV_I))
  
  # If V contains NA (exploding case) just return 1
  if (any(is.na(NewV))){
    NewV <- paramsBellman$OldV+1
    Vtilde <- paramsBellman$OldV+1
    Investprob <- paramsBellman$OldV*0
  }
  
  results$NewV <- NewV
  results$Vtilde <- Vtilde
  results$Investprob <- Investprob
  
  return(results)
}

#' Compute the probability of each one of the five possible transitions for each
#' omega state. Transitions are pairs (current violation, violator status)
#' @param CaseProb Probability matrix for each (omega, case) pair
#' @param vioCase  Matrix of current violation for each (omega, case) pair
#' @param ordvioCase Matrix of violator status for each (omega, case) pair
#' @param cv current violation in this transition (either 0 or 1)
#' @param vs violator status in this transition (either 0, 1 or 2)
CollapseProbByTransitions <- function(CaseProb,vioCase,ordvioCase,cv,vs){
  
  Prob <- CaseProb
  Prob[ordvioCase != (vs+1)] <- 0
  if (vs > 0){
    Prob[vioCase != cv] <- 0
  }
  Prob <- rowSums(Prob)
  return(Prob)
}

#' Get the index of each omega_tilde state after the transition given by 
#' (current violation, violator status)
#' @param omega_states omega states before transition
#' @param cv current violation in this transition (either 0 or 1)
#' @param vs violator status in this transition (either 0, 1 or 2)
getIndexOfTransition <- function(omega_states,cv,vs){
  omega_states[,3] <- cv
  omega_states[,4] <- vs
  return(apply(omega_states,1,ValuesToIndex))
}

#' Checks the value function with the provided file
CompareValues <- function(outBellman,file,params){
  
  library(ggplot2)
  library(tidyr)
  
  N1 <- params$N1
  n_omega1 <- params$n_omega1
  nstates <- n_omega1*N1
  myResults <- outBellman$NewV
  
  expectedResults <- read.csv(file)
  expectedResults <- expectedResults[1:nstates,"value"]
  
  diff <- max(abs(myResults-expectedResults))
  
  print(paste0("the difference is ",diff))
  
  
  index <- c(1:nstates)
  plotData <- data.frame(myResults,expectedResults,index)
  
  plotData %>%
    gather(key,value, myResults, expectedResults) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
  
}

#' Checks the valueTilde function with the provided file
CompareVTilde <- function(outBellman,file,params){
  
  library(ggplot2)
  library(tidyr)
  
  N1 <- params$N1
  n_omega1 <- params$n_omega1
  nstates <- n_omega1*N1
  myResults <- outBellman$Vtilde
  
  expectedResults <- read.csv(file)
  expectedResults$compliance <- 1*(expectedResults$orderedvio>0)
  expectedResults <- expectedResults[order(expectedResults$NAICS,
                                           expectedResults$region,
                                           expectedResults$gravity,
                                           expectedResults$compliance,
                                           expectedResults$laginv,
                                           expectedResults$violation,
                                           expectedResults$orderedvio,
                                           expectedResults$DAVgrid),]
  expectedResults <- expectedResults[1:nstates,"value"]
  
  diff <- max(abs(myResults-expectedResults))
  
  print(paste0("the difference is ",diff))
  
  
  index <- c(1:nstates)
  plotData <- data.frame(myResults,expectedResults,index)
  
  plotData %>%
    gather(key,value, myResults, expectedResults) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
  
}

#' Checks the investment probabilities and likelihoods with the provided file
CompareProb <- function(LLdata,file){

  library(ggplot2)
  library(tidyr)
  
  LLdata$obs2 <- c(1:nrow(LLdata))
  LLdata <- LLdata[LLdata$viostatus>0,]
  
  
  expectedResults <- read.csv(file)
 
  merged <- cbind(LLdata,expectedResults)
  tol <- 1e-3
  #merged <- merge(x = expectedResults, y = LLdata, by = "obs", all.x = TRUE)
  merged$probdiff <- abs(merged$invprob - merged$prob)>tol
  merged <- merged[,c("obs","obs2","o1","li1","li2","viostatus","violation","DAV","x",
                      "invprob","prob","probdiff","loglike","l_i")]
  
  merged$DAV0 <- 1*(merged$DAV == 0)
  nobs <- nrow(merged)
  ndiff <- sum(merged$probdiff)             
  
  
  diff <- max(abs(merged$invprob-merged$prob))
  
  print(paste0(ndiff," observations out of ",nobs," have different prob. values"))
  print(paste0("the difference in probailities is ",diff))
  diff <- max(abs(merged$loglike-merged$l_i))
  print(paste0("the difference in log-likelihood is ",diff))
  
  p<-ggplot(merged,aes(x=invprob, y=prob)) +
    geom_point()+xlab("Investment probability (expected result)")+
    ylab("Investment probability (my result)") + xlim(0,1)+ylim(0,1)
  show(p)
  return(merged)
}

#' Recover the likelihood for states with DAV = 2
Question1C <- function(outBellman,file){
  
  index <- c(1:length(outBellman$Investprob))
  states <- t(sapply(index,IndexToValues)) 
  dataout <- data.frame(states)
  colnames(dataout) <- c("omega1","lag1_inv","violation","viostatus","intDAV")
  dataout$DAV <- (dataout$intDAV-1)*0.5 
  dataout$invprob <- outBellman$Investprob
  dataout <- dataout[dataout$DAV==2,]
  write.csv(dataout,file, row.names = FALSE)
}

#' Add BellmanParams and LogLikeParams structures to params for faster computations
AddParams <- function(params){
  params$BellmanParams <- PackBellmanParams(params)
  params$LogLikeParams <- LogLikeParams(params)
  return(params)
}

#' Estimates the variance of the ML estimator using the outer product approximation
#' @param thetaML ML estimator
#' @param params Estimation parameters
#' @DeltaTHeta Difference used in numerical approximation of the gradient
#' @param V0 Initial value for the value function
EstimateMLVariance <- function(thetaML, params,DeltaTheta,V0){
  
  # Construct matrices to compute the Bellman equation
  if (!exists("BellmanParams",where = params)){
    params$BellmanParams <- PackBellmanParams(params)
  }
  
  # Estimate likehood at thetaML
  if (missing(V0)){
    outML <- BellmanLoop(thetaML, params)  
  }
  else {
    outML <- BellmanLoop(thetaML, params,V0)
  }
  
  # Compute the loglikelihood at the ML parameter
  LLOutput <- LogLike(outML,params)
  LLData <- LLOutput$lldata
  LLData$LL_ML <- LLData$l_i  
  
  # Estimate likelihood at small deviations
  for (i in 1:length(thetaML)){
    theta_i <- thetaML
    theta_i[i] <- theta_i[i] - DeltaTheta # Delta-deviation
    out <- BellmanLoop(theta_i, params,outML$NewV)
    LLOutput <- LogLike(out,params)
    varname <- paste0("dlog_",i)
    LLData[,varname] <- LLData$LL_ML-LLOutput$lldata$l_i
  }
  
  # Compute the outer product of the gradient of the log likelihood function
  GradProd <- matrix(0,length(thetaML),length(thetaML))
  Ncols <- ncol(LLData)
  LLData<- as.matrix(LLData[,c((Ncols + 1 - length(thetaML)):Ncols)])
  
  for (i in 1:nrow(LLData)){
    GradProd = GradProd + (t(t(LLData[i,]))%*%LLData[i,])/(DeltaTheta^2)
  }
  return(solve(GradProd))
}
 
#' Computes the Log likelihood at all observations. 
#' Also returns investment probabilities
#' @param BellmanOut Output from Bellman function
#' @param params Estimation parameters
LogLike <- function(BellmanOut,params){
  
  
  if (!exists("LogLikeParams",where = params)){
    params$LogLikeParams <- LogLikeParams(params)
  }
  
  i1 <- params$LogLikeParams$i1
  i2 <- params$LogLikeParams$i2
  w_below <- params$LogLikeParams$w_below
  w_above <- params$LogLikeParams$w_above
  x <- params$LogLikeParams$x
  o1 <- params$LogLikeParams$o1
  li1 <- params$LogLikeParams$li1
  li2 <- params$LogLikeParams$li2
  viostatus <- params$LogLikeParams$viostatus
  violation <- params$LogLikeParams$violation
  DAV <- params$LogLikeParams$DAV

  p1 <- BellmanOut$Investprob[i1]
  p2 <- BellmanOut$Investprob[i2]
  
  prob <- w_below*p1 + w_above*p2
  
  l_i <- log(x*prob + (1-x)*(1-prob))
  ll <- sum(l_i)

  out <- list()
  out$ll <- ll
  
  lldata <- data.frame(o1,li1,li2,viostatus,violation,DAV,x,prob,l_i)
  
  out$lldata <- lldata
  return(out)
}

#' Defines matrices and indexes that allow faster repeated computation of the 
#' log likelihood function
LogLikeParams <- function(params){
  
  DAVgrid <-params$DAVgrid
  data <-params$panel_data
  n_obs <-nrow(data)
  npointsDAV  <- length(DAVgrid)
  maxpointDAV <- DAVgrid[npointsDAV]
  widthDAV <- DAVgrid[2]-DAVgrid[1]
  
  DAV <- data$DAV
  x   <- data$investment
  o1  <- data$omega1
  li1 <- data$lag_investment
  li2 <- data$lag2_investment
  viostatus <- data$ordered_violator
  violation <- data$violation
  
  # Interpolate newDAV using the grid
  # indexes in DAV grid
  i_below <- pmin(1 + floor(DAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  i_above <- pmin(2 + floor(DAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  # weights
  w_below <- pmax((DAVgrid[i_above]-DAV)/widthDAV,0)
  w_above <- 1-w_below
  # states
  state1 <- cbind(o1,li1,violation,viostatus,i_below)
  state2 <- cbind(o1,li1,violation,viostatus,i_above)
  # indexes in probability matrix
  i1 <- apply(state1,1,ValuesToIndex)
  i2 <- apply(state2,1,ValuesToIndex)
  
  llparams <- list()
  llparams$i1 <- i1
  llparams$i2 <- i2
  llparams$w_below <- w_below
  llparams$w_above <- w_above
  llparams$x <- x
  llparams$o1 <- o1
  llparams$li1 <- li1
  llparams$li2 <- li2
  llparams$viostatus <- viostatus
  llparams$violation <- violation
  llparams$DAV <- DAV
  
  
  return(llparams)
}

#' Computes the Log likelihood at all observations. 
#' Faster than Loglike, but does not output investment probabilities
#' @param BellmanOut Output from Bellman function
#' @param params Estimation parameters
LogLike2 <- function(BellmanOut,params){

  i1 <- params$LogLikeParams$i1
  i2 <- params$LogLikeParams$i2
  w_below <- params$LogLikeParams$w_below
  w_above <- params$LogLikeParams$w_above
  x <- params$LogLikeParams$x
  
  p1 <- BellmanOut$Investprob[i1]
  p2 <- BellmanOut$Investprob[i2]
  
  prob <- w_below*p1 + w_above*p2
  
  l_i <- log(x*prob + (1-x)*(1-prob))
  ll <- sum(l_i)
  return(ll)
}

#' Function to be maximized in ML estimation
#' @param theta Coefficient at which the value function will be found and then
#' the loglikelihood will be computed
#' @param params Estimation parameters 
EvalFunction <- function(theta,params){

  print(theta)
  
  # Find the value function
  if (exists("V0")){
    outBellman <- BellmanLoop(theta, params,V0)
  }
  else {
    outBellman <- BellmanLoop(theta, params)
  }
  assign("V0",outBellman$NewV,envir = .GlobalEnv)
  
  # Estimate the log likelihood
  LL <- LogLike2(outBellman,params)
  return(LL)
}

#' Nested fixed point estimation
#' @param theta0 Initial value 
#' @param params Estimation parameters
#' @param solver Solver choice
NestedFixedPoint<-function(theta0,params,solver){
  
  # Construct matrices to compute the Bellman equation
  if (!exists("BellmanParams",where = params)){
    params$BellmanParams <- PackBellmanParams(params)
  }
  
  if (!exists("LogLikeParams",where = params)){
    params$LogLikeParams <- LogLikeParams(params)
  }
  
  # Solver selection
  if (solver== 1){
    library(pracma)
    out<-fminsearch(EvalFunction, theta0, params = params,
                    minimize= FALSE, method = "Hooke-Jeeves", maxiter = 500, tol = 1e-6)
  }
  else if (solver == 2){
    library(pracma)
    out<-fminsearch(EvalFunction, theta0, params = params,
                    minimize= FALSE, method = "Nelder-Mead", maxiter = 500, tol = 1e-6)
  }
  else if (solver == 3){
    library(optimx)
    out <- optimx(theta0, EvalFunction, gr=NULL, hess=NULL, lower=c(0,-Inf,-Inf,-Inf,-Inf),
                  upper=c(Inf,0,0,0,0), 
           method="CG", itnmax=500, hessian=FALSE,
           control=list(maximize=TRUE),
           params = params)
  }
  else if (solver == 4){
    library(optimx)
    out <- optimx(theta0, EvalFunction, gr=NULL, hess=NULL, lower=c(0,-Inf,-Inf,-Inf,-Inf),
                  upper=c(Inf,0,0,0,0), 
                  method="L-BFGS-B", itnmax=500, hessian=FALSE,
                  control=list(maximize=TRUE),
                  params = params)
  }
  else{
    library(maxLik)
    
    A <- matrix(c(0,0,0,0,-1), 1, 5)
    B <- 0
    out<- maxNR(EvalFunction, grad = NULL, hess = NULL, start = theta0, print.level = 1,
                tol = 1e-06, reltol=sqrt(.Machine$double.eps), gradtol = 1e-06,
                steptol = 1e-10, lambdatol = 1e-06, qrtol = 1e-10,
                iterlim = 500,
                constraints=list(eqA=A, eqB=B),
                params = params)
  }
  
  return(out)
}

#' Mapping from indexes (of the vector representing the value function) and the
#' corresponding omega states
IndexToValues <- function(index){
  
  N1 <- 161
  N2 <- 80
  N3 <- 40
  N4 <- 20
  
  omega1 <- 1 + (index-1)%/%N1
  index <- (index-1)%%N1-1 # index goes from -1 to N1-2
  violator <- index>-1
  lagi1 <- (index%/%N2)*violator
  index <- index%%N2  # index goes from 0 to N2-1
  lagi2 <- (index%/%N3)*violator
  index <- index%%N3
  vio   <- (index%/%N4+1)*violator
  index <- index%%N4+1
  dav   <- index*violator + (violator==0)
  return(c(omega1,lagi1,lagi2,vio,dav))
}

#' Mapping from omega states to indexes (of the vector representing the value 
#' function) 
ValuesToIndex <- function(x){
  
  N1 <- 161
  N2 <- 80
  N3 <- 40
  N4 <- 20
  
  i_o1  <- x[1]
  i_l1  <- x[2]
  i_l2  <- x[3]
  i_vio <- x[4]
  i_DAV <- x[5]
  
  i <- (i_o1-1)*N1 + 1 + (i_vio>=1)*((i_l1)*N2 + (i_l2)*N3 + (i_vio-1)*N4 + i_DAV)
  return(i)
}

#' Mapping from omega states to indexes (of the vector representing the value 
#' function) 
ValuesToIndexExclCompliance <- function(x){
  
  N1 <- 160
  N2 <- 80
  N3 <- 40
  N4 <- 20
  
  i_o1  <- x[1]
  i_l1  <- x[2]
  i_l2  <- x[3]
  i_vio <- x[4]
  i_DAV <- x[5]
  
  i <- (i_o1-1)*N1 + (i_vio>=1)*((i_l1)*N2 + (i_l2)*N3 + (i_vio-1)*N4 + i_DAV)
  return(i)
}

# Transition matrices from omega_tilde to omega. There are four possible transitions
# with and without investment, and DAVgrid above and below the actual DAV.
# This transition matrix is the same for each omega1
InvestmentTransitionMatrices <- function(params){
  N1 <- params$N1
  paramsBellman <- params$BellmanParams
  w_below <- paramsBellman$w_below[1:N1]
  w_above <- paramsBellman$w_above[1:N1]
  indexI_below <- paramsBellman$indexI_below[1:N1]
  indexNI_below <- paramsBellman$indexNI_below[1:N1]
  indexI_above <- paramsBellman$indexI_above[1:N1]
  indexNI_above <- paramsBellman$indexNI_above[1:N1]
  
  A <- matrix(0,N1,N1)
  B <- matrix(0,N1,N1)
  
  for (i in 1:N1){
    A[i,indexI_below[i]]  <- A[i,indexI_below[i]] + w_below[i]
    A[i,indexI_above[i]]  <- A[i,indexI_above[i]] + w_above[i]
    B[i,indexNI_below[i]] <- B[i,indexNI_below[i]] + w_below[i]
    B[i,indexNI_above[i]] <- B[i,indexNI_above[i]] + w_above[i]
  }
  
  params$InvTr <- A
  params$NoInvTr <- B
  return(params)
}

# Transition from omega to the next omega_tilde. Transition probabilities
# depend on the regulator CCPs. We compute a transition matrix for each oemga1
CCPTransitionMatrix <- function(params){
  N1 <- params$N1
  paramsBellman <- params$BellmanParams
  n_omega1 <- params$n_omega1
  nstates <- n_omega1*N1   
  
  indexTransition <- paramsBellman$indexTransition
  ProbTransition  <- paramsBellman$ProbTransition
  
  indexTransition  <- RetainN1States(indexTransition,nstates,N1,5,1)
  
  CCPTr <- list()
  
  for (o1 in 1:n_omega1){
    ProbTransition_o1   <- RetainN1States(ProbTransition,nstates,N1,5,o1)
    
    C <- matrix(0,N1,N1)
    for (i in 1:N1){
      for (t in 1:5){
        C[i,indexTransition[i,t]] <- ProbTransition_o1[i,t]
      }
    }
    
    CCPTr[[o1]] <- C
  }
  
  params$CCPTr <- CCPTr
  return(params)
}

# Given a value of omega1, return the corresponding NAICS, region and gravity
GetOmega1Mapping <- function(params){
  N1 <- params$N1
  n_omega1 <- params$n_omega1
  index <- 1 + c(0:(n_omega1-1))*N1
  data <- params$data[index,c("omega1","naics_recode","region","gravity")]
  params$omega1Data <- data
  return(params)
}

# Add transition matrices to the params structure
AddTransitionParams <- function(params){
  params <- InvestmentTransitionMatrices(params)
  params <- CCPTransitionMatrix(params)
  return(params)
}

# Helper function. Returns the N1 first states from vec
RetainN1States <- function(vec,nstates,N1,N2,o1){
  vec <- matrix(vec,nstates,N2)
  vec <- vec[((o1-1)*N1+1):(N1*o1),]
  return(vec)
}

# Given omega1, compute the transition probability fron one omega_tilde
# to the next omega_tilde
GetTransitionProb <- function(params,omega1,InvProb){
  N1 <- params$N1
  InvProb <- InvProb[(N1*(omega1-1)+1):(N1*omega1)]
  C <- params$CCPTr[[omega1]]
  P <- GetTransitionProbOmega1(params,InvProb,C)
  return(P)
}

# Given the transition matrices, compute the transition probability matrix 
# fron one omega_tilde to the next omega_tilde
GetTransitionProbOmega1 <- function(params,InvProb,C){
  N1 <- params$N1
  D <- diag(c(InvProb))
  P  <- (params$NoInvTr + D%*%(params$InvTr-params$NoInvTr))%*%C
  return(P)
}

# Check results with the provided file
CompareTransitionProb <- function(P_omega1_BGL,ExpectedResult1,params){
  expectedResults <- read.csv(ExpectedResult1)
  N1 <- params$N1
  state <- c(1,1,1,2,3)
  index <- ValuesToIndex(state)
  myResult <- c(P_omega1_BGL[index,])
  states <- params$BellmanParams$states[1:N1,]
  results <- data.frame(states)
  colnames(results) <- c("omega1","lag1_inv","current_vio","ord_vio","intDAV")
  results$myResult  <- myResult
  
  results <- results[
    order(results$omega1,results$lag1_inv, results$ord_vio,results$current_vio,
          results$intDAV),]
  results$expectedResult <- expectedResults$transprob
  
  diff <- max(abs(results$myResult-results$expectedResult))
  
  print(paste0("the difference is ",diff))
  
  
  results$index <- c(1:params$N1)
  results %>%
    gather(key,value, myResult, expectedResult) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
}

# Prepare answer for question 1 in problem 3
Question_P3Q1 <- function(P_omega1_BGL,file,params){
  N1 <- params$N1
  state <- c(1,1,1,2,5)
  index <- ValuesToIndex(state)
  myResult <- c(P_omega1_BGL[index,])
  states <- params$BellmanParams$states[1:N1,]
  results <- data.frame(states)
  colnames(results) <- c("omega1","lag1_inv","current_vio","ord_vio","intDAV")
  results$DAV <- (results$intDAV-1)*0.5 
  results$transProb  <- myResult
  
  results <- results[
    order(results$omega1,results$lag1_inv, results$ord_vio,results$current_vio,
          results$intDAV),]
  write.csv(results,file, row.names = FALSE)
}

# Compute the steady state distribution given the transition matrix P
SolveSS <- function(P,params){
  N1 <- params$N1

  A <- rbind(t(P)-diag(N1),matrix(1,1,N1))
  b <- rbind(matrix(0,N1,1),1)

  pi <- lm(b ~ A +0)$coefficients
  
  return(as.matrix(pi))
}

# Check the results with the provided file
CompareSS <- function(piSS,ExpectedResult,params){
  expectedResults <- read.csv(ExpectedResult)
  N1 <- params$N1
  states <- params$BellmanParams$states[1:N1,]
  states <- states[states[,5]==3,]
  
  index <- apply(states,1,ValuesToIndex)
  myResult <- c(piSS[index,])
  
  results <- data.frame(states)
  colnames(results) <- c("omega1","lag1_inv","current_vio","ord_vio","intDAV")
  results$myResult  <- myResult
  
  results <- results[
    order(results$omega1,results$lag1_inv, results$ord_vio), ]
  results$expectedResult <- expectedResults$transprob
  
  diff <- max(abs(results$myResult-results$expectedResult))
  
  print(paste0("the difference is ",diff))
  
  results$index <- c(1:nrow(states))
  results %>%
    gather(key,value, myResult, expectedResult) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
}

# Prepare solution to question 2 in problem 3
Question_P3Q2 <- function(piSS,outFile,params){

  N1 <- params$N1
  states <- params$BellmanParams$states[1:N1,]
  states <- states[states[,5]==5,]
  
  index <- apply(states,1,ValuesToIndex)
  myResult <- c(piSS[index,])
  
  results <- data.frame(states)
  colnames(results) <- c("omega1","lag1_inv","current_vio","ord_vio","intDAV")
  results$DAV <- (results$intDAV-1)*0.5 
  results$ssProbability  <- myResult
  
  results <- results[
    order(results$omega1,results$lag1_inv, results$ord_vio), ]
  
  write.csv(results,outFile, row.names = FALSE)
}

# Compute the sets of model moments given a parameter
ComputeMoments <- function(coeff, params,V0){
  
  n_omega1 <- params$n_omega1
  N1       <- params$N1
  n_states <- n_omega1*N1
  
  # Run the loop to find the fixed point of the Bellman equation
  if(!missing(V0)) {
    bellmanSol <- BellmanLoop(coeff, params,V0)
  }
  else {
    bellmanSol <- BellmanLoop(coeff, params)
  }
  
  invProb <- as.matrix(bellmanSol$Investprob)
  
  moments1 <- matrix(0,n_states,1)
  moments2 <- matrix(0,n_states-n_omega1,1)
  
  for (o1 in 1:n_omega1){
    P <- GetTransitionProb(params,o1,invProb)
    piSS <- SolveSS(P,params)
    
    # moment 1: steady state distribution
    moments1[((o1-1)*N1+1):(o1*N1),1] <- piSS
    
    # moment 2: probabi;ity of investment in steady state distribution
    m2 <- piSS*invProb[((o1-1)*N1+1):(o1*N1),1]
    moments2[((o1-1)*(N1-1)+1):(o1*(N1-1)),1] <- m2[2:N1,1]
  }
  out <-list()
  out$moments <- rbind(moments1,moments2)
  out$V0 <- bellmanSol$NewV
  
  return(out)
}

# Retrieve the omega_tilde state for all the moments
GetStatesForMoments <- function(params){
  
  n_omega1 <- params$n_omega1
  N1       <- params$N1
  n_states <- n_omega1*N1
  states <- params$BellmanParams$states
  
  states1 <- states
  states2 <- matrix(0,n_states-n_omega1,5)
  
  for (o1 in 1:n_omega1){
    m2 <- states[((o1-1)*N1+1):(o1*N1),]
    states2[((o1-1)*(N1-1)+1):(o1*(N1-1)),] <- m2[2:N1,]
  }
  
  mstates <- rbind(states1,states2)
  mstates <- cbind(rbind(matrix(1,n_states,1),matrix(2,n_states-n_omega1,1)),mstates)
  return(mstates)
}

# Check the results with the provided file
CompareMoments <- function(moments,ExpectedResult,params){
  library(tidyr)
  library(ggplot2)
  expectedResults <- read.csv(ExpectedResult)
  N1 <- params$N1
  states <- params$BellmanParams$states[1:N1,]
  
  mstates <- GetStatesForMoments(params)
  results <- data.frame(mstates)
  colnames(results) <- c("momset","omega1","lag1_inv","current_vio","ord_vio","intDAV")
  results$DAV <- (results$intDAV-1)*0.5 
  results$myMoments <- moments
  
  results$naics_recode <- params$omega1Data[results$omega1,"naics_recode"]
  results$region <- params$omega1Data[results$omega1,"region"]
  results$gravity <- params$omega1Data[results$omega1,"gravity"]
  
  # sort columns
  results <- results[,c("momset","naics_recode","region","gravity","intDAV","DAV",
                        "ord_vio","lag1_inv","current_vio","myMoments")]
  
  # sort rows (momset2 follow a different order in the file)
  results1 <- results[results$momset == 1,]
  results2 <- results[results$momset == 2,]
  
  results1 <- results1[
    order(results1$momset,results1$naics_recode,results1$region,results1$gravity,
          results1$lag1_inv, results1$ord_vio,results1$current_vio,
          results1$intDAV), ]
  
  results2 <- results2[
    order(results2$momset,results2$naics_recode,results2$region,results2$gravity,
          results2$lag1_inv, results2$ord_vio,results2$current_vio,
          results2$intDAV), ]
  
  
  # Add momnumber
  results1$momnumber <- c(1:nrow(results1))-1
  results2$momnumber <- c(1:nrow(results2))-1
  
  # join the results
  results <- rbind(results1,results2)
  
  # filter by DAV
  results <- results[results$intDAV==3,]
  
  # Add expected results
  results$expectedMoments <- expectedResults$modelmoment1
  
  write.csv(results,"Output/moments_test.csv")
  
  diff <- max(abs(results$myMoments-results$expectedMoments))
  
  print(paste0("the difference is ",diff))
  
  results$index <- c(1:nrow(results))
  results %>%
    gather(key,value, myMoments, expectedMoments) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
}

# Prepare solution to question 3 in part 3
Question_P3Q3 <- function(parameterGrid,params,file){
  Nparams <- 10
  results <- matrix(0,Nparams,5+3)
  for (i in 1:Nparams){
    theta <- parameterGrid[i,]
    if (i==1){
      out <- ComputeMoments(theta, params)
    }
    else {
      out <- ComputeMoments(theta, params,V0)
    }
    moments <-out$moments
    V0 <- out$V0
    
    results[i,1:5] <- c(theta)
    results[i,5+1] <- mean(moments)
    results[i,5+2] <- sd(moments)
    results[i,5+3] <- nrow(moments)
  }
  results <- data.frame(results)
  colnames(results) <- c("inv_cost","inspect_benefit","violation_benefit",
                         "fine_benefit","hpv_benefit",
                         "moment_mean","moment_sd","n_moments")
  write.csv(results,file,row.names = FALSE)
}

# Do DAV interpolation
GetDAVInterpolation <- function(params,DAV){
  DAVgrid <- params$DAVgrid  # retrieve DAV grid
  npointsDAV  <- length(DAVgrid) # number of points in the grid
  maxpointDAV <- DAVgrid[npointsDAV] # maximum DAV value
  widthDAV <- DAVgrid[2]-DAVgrid[1]  # grid width
  
  # Interpolate newDAV using the grid: find interpolation DAVgrid indexes and 
  # weights
  
  # DAVgrid indexes
  i_below <- pmin(1 + floor(DAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  i_above <- pmin(2 + floor(DAV*(npointsDAV-1)/maxpointDAV),npointsDAV)
  # weights
  w_below <- pmax((DAVgrid[i_above]-DAV)/widthDAV,0)
  w_above <- 1-w_below
  
  out <- c(i_below,i_above,w_below,w_above)
  return(out)
}

# Compute data moments
ComputeDataMoments <- function(params){
  N1 <- params$N1
  n_omega1 <- params$n_omega1
  nstates <- N1*n_omega1
  panel <- params$panel_data
  cols <- c("lag_investment","violation","ordered_violator","DAV","investment")
  
  moment1 <- matrix(0,nstates,1)
  moment2 <- matrix(0,nstates-n_omega1,1)
  
  for (o1 in 1:n_omega1){
    panel_o1 <- panel[panel$omega1 == o1,cols]
    n_obs <- nrow(panel_o1) # number of states in omega1
    n_obs_nc <- nrow(panel_o1[panel_o1$ordered_violator > 0,]) # number of states not in compliance
    for (i in 1:n_obs){
      omega2 <- panel_o1[i,]
      DAV <- omega2[1,4]
      
      # Interpolate DAV
      dav_interpolation <- GetDAVInterpolation(params,DAV)
      
      # Recover the two states related to DAV
      state1 <- c(o1,omega2[1,1],omega2[1,2],omega2[1,3],dav_interpolation[1])
      state2 <- c(o1,omega2[1,1],omega2[1,2],omega2[1,3],dav_interpolation[2])
      
      # Get the indexes of these two states
      index1 <- ValuesToIndex(state1)
      index2 <- ValuesToIndex(state2)
      
      # Add relative frequency of state omega_2
      moment1[index1,1] <- moment1[index1,1] + dav_interpolation[3]/n_obs
      moment1[index2,1] <- moment1[index2,1] + dav_interpolation[4]/n_obs
      
      if (omega2[1,5]>0) { # if the plant invests
        
        # Get the indexes of these two states
        index1 <- ValuesToIndexExclCompliance(state1)
        index2 <- ValuesToIndexExclCompliance(state2)
        
        # Add relative frequency of state with investment
        moment2[index1,1] <- moment2[index1,1] + dav_interpolation[3]/(n_obs)
        moment2[index2,1] <- moment2[index2,1] + dav_interpolation[4]/(n_obs)
      }
    }
  }
  
  moments <- rbind(moment1,moment2)
  return(moments)
}

# Check results with the provided file
CompareDataMoments <- function(moments,ExpectedResult,params){
  library(tidyr)
  library(ggplot2)
  expectedResults <- read.csv(ExpectedResult)
  N1 <- params$N1
  states <- params$BellmanParams$states[1:N1,]
  
  mstates <- GetStatesForMoments(params)
  results <- data.frame(mstates)
  colnames(results) <- c("momset","omega1","lag1_inv","current_vio","ord_vio","intDAV")
  results$DAV <- (results$intDAV-1)*0.5 
  results$myMoments <- moments
  
  results$naics_recode <- params$omega1Data[results$omega1,"naics_recode"]
  results$region <- params$omega1Data[results$omega1,"region"]
  results$gravity <- params$omega1Data[results$omega1,"gravity"]
  
  # sort columns
  results <- results[,c("momset","naics_recode","region","gravity","intDAV","DAV",
                        "ord_vio","lag1_inv","current_vio","myMoments")]
  
  # sort rows (momset2 follow a different order in the file)
  results1 <- results[results$momset == 1,]
  results2 <- results[results$momset == 2,]
  
  results1 <- results1[
    order(results1$momset,results1$naics_recode,results1$region,results1$gravity,
          results1$lag1_inv, results1$ord_vio,results1$current_vio,
          results1$intDAV), ]
  
  results2 <- results2[
    order(results2$momset,results2$naics_recode,results2$region,results2$gravity,
          results2$lag1_inv, results2$ord_vio,results2$current_vio,
          results2$intDAV), ]
  
  
  # Add momnumber
  results1$momnumber <- c(1:nrow(results1))-1
  results2$momnumber <- c(1:nrow(results2))-1
  
  # join the results
  results <- rbind(results1,results2)
  
  # filter by DAV
  results <- results[results$intDAV==3,]

  # Add expected results
  results$expectedMoments <- expectedResults$datamoment
  
  diff <- max(abs(results$myMoments-results$expectedMoments))
  
  print(paste0("the difference is ",diff))
  
  results$index <- c(1:nrow(results))
  results %>%
    gather(key,value, myMoments, expectedMoments) %>%
    ggplot(aes(x=index, y=value, colour=key)) +
    geom_line()
}

# Prepare solution to question 4, part 3
Question_P3Q4 <- function(dataMoments,params,file){
  Nparams <- 10
  results <- matrix(0,1,3)
  
  results[1,1] <- mean(dataMoments)
  results[1,2] <- sd(dataMoments)
  results[1,3] <- nrow(dataMoments)
  results <- data.frame(results)
  colnames(results) <- c("moment_mean","moment_sd","n_moments")
  write.csv(results,file,row.names = FALSE)
}

# Prepare solution to question 5, part 3
GenerateMomentFiles <- function(parameterGrid,dataMoment,params,file1,file2){
  Nparams <- 500
  Nmoments <- 14445
  results <- matrix(0,Nmoments,Nparams+1)
  
  # Add data moments
  results[,1] <- c(dataMoment)
  
  # Add moments for all parameters 
  colname <- c()
  for (i in 1:Nparams){
    theta <- parameterGrid[i,]
    if (i==1){
      out <- ComputeMoments(theta, params)
    }
    else {
      out <- ComputeMoments(theta, params,V0)
    }
    moments <-out$moments
    V0 <- out$V0
    
    results[,i+1] <- c(moments)
    
    colname <- c(colname, paste0("param_",i))
    
    print(paste0("Added moments for parameter number ",i))
  }
  
  results <- data.frame(results)
  colnames(results) <- c("data",colname)
  
  mstates <- GetStatesForMoments(params)
  first_columns <- data.frame(mstates)
  colnames(first_columns) <- c("momset","omega1","lag1_inv","current_vio",
                               "ord_vio","intDAV")
  
  results <- cbind(first_columns,results)
  
  results <- results[
    order(results$momset,results$omega1,
          results$lag1_inv, results$ord_vio,results$current_vio,
          results$intDAV), ]
  
  write.csv(results,file2,row.names = FALSE,)
  write.csv(results[,7:(Nparams+7)],file1,row.names = FALSE,col.names = FALSE)
}

