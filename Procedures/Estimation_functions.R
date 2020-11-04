theta.to.param <- function(theta, nb_states){
  
  daystorecover <- 14
  delta         <- 0.6/(daystorecover * 100)
  gamma         <- 1/daystorecover
  
  Parameters <- list("delta" = delta, # Death rate
                     "gamma" = gamma, # Recovery rate
                     "kappa" = 0.001,
                     "beta"  = plogis(theta[1]),
                     "sigma" = exp(theta[2])*10,
                     "theta_F_high" = rep(1, nb_states),
                     "theta_F_low"  = rep(.1, nb_states),
                     "theta_I_high" = rep(1, nb_states),
                     "theta_I_low"  = rep(.64, nb_states),
                     "start.beta"   = plogis(theta[1]),
                     "Mask_reduc" = 0.58,
                     "rho" = plogis(theta[3]),
                     "tau_com" = (5*8)/(16*7),
                     "tau_trav" = 4,
                     "meas_std_death" = sqrt(.000001),
                     "meas_std_pop" = NA)
  

  Parameters$Premult_Beta <- diag(1, nb_states)
  Parameters$Sigma_mat <- diag(sqrt(Parameters$sigma), nb_states) %*%
    (diag(1-Parameters$rho, nb_states) + Parameters$rho * matrix(1, nb_states, nb_states)) %*% diag(sqrt(Parameters$sigma), nb_states)
  
  return(Parameters)
}


Loglik <- function(theta, 
                   Death.selec,
                   Dummy_theta_F,
                   Dummy_theta_I,
                   Dummy_masks,
                   Weight_commuting,
                   Weight_migration,
                   Pop_vector,
                   for_estim = T){
  
  # Retrieving the number of states
  nb_states <- ncol(Death.selec)
  
  # Backing out the parameters
  Parameters <- theta.to.param(theta, nb_states)
  
  # Initializing the filter
  #------------------------
  # Initialization of factors
  W0 <- c(rep(0, nb_states), Pop_vector, rep(10, nb_states), rep(0, nb_states), 
          rep(Parameters$start.beta, nb_states))
  
  P0 <- diag(0, 5 * nb_states)
  # block of S
  P0[(nb_states + 1):(2 * nb_states), (nb_states + 1):(2 * nb_states)]          <- diag(500, nb_states)#0.01/100 * diag(Pop_vector, nb_states)
  # block of I
  P0[(2 * nb_states + 1):(3 * nb_states), (2 * nb_states + 1):(3 * nb_states)]  <- diag(500, nb_states)#0.01/100 * diag(Pop_vector, nb_states)
  # Covariance terms 
  P0[(nb_states + 1):(2 * nb_states), (2 * nb_states + 1):(3 * nb_states)]      <- - 0.5 * (
    P0[(nb_states + 1):(2 * nb_states), (nb_states + 1):(2 * nb_states)] + 
      P0[(2 * nb_states + 1):(3 * nb_states), (2 * nb_states + 1):(3 * nb_states)])
  P0[(2 * nb_states + 1):(3 * nb_states), (nb_states + 1):(2 * nb_states)]      <- 
    t(P0[(nb_states + 1):(2 * nb_states), (2 * nb_states + 1):(3 * nb_states)] )

  
  P0[(4 * nb_states + 1):(5 * nb_states),(4 * nb_states + 1):(5 * nb_states)]   <- diag(.001^2, nb_states)
    
  
  # Performing filtering
  Filtrage_EKF_estim <- EKF_filter_cpp(Parameters,
                                       W0,
                                       P0,
                                       Death.selec,
                                       Dummy_theta_F,
                                       Dummy_theta_I,
                                       Dummy_masks,
                                       Weight_commuting,
                                       Weight_migration,
                                       Pop_vector)
  
  if(for_estim == T){
    #This Penalty controls for too large values with the parameter constraints
    penalty = 1e5 * sum((theta - 15)^2 * (abs(theta)>15))
    
    return(-sum(Filtrage_EKF_estim$loglik.vector[1:length(Filtrage_EKF_estim$loglik.vector)])
           + penalty)
  } else{
    return(Filtrage_EKF_estim)
  }
}




