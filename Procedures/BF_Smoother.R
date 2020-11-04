#########################################
# SMOOTHER
#########################################

Modified_Branzon_Frazier = function(W.up,
                                    W.pred,
                                    P.up,
                                    P.pred,
                                    Jacobians,
                                    Pop_vector,
                                    Kalman_Gain,
                                    M.pred,
                                    y_residual){
  
  require(MASS)
  require(stringr)
  
  # Objects to retrieve
  nb.states <- length(Pop_vector)
  time      <- ncol(W.up)
  
  
  
  W.smooth = NA * W.up
  P.smooth = list()
  
  
  
  # Initialization
  
  W.smooth[, time] = W.up[, time]
  P.smooth[[time]] = P.up[[time]]
  
  H_matrix <- matrix(0, nrow = nb.states, ncol = nrow(W.up))
  H_matrix[1:nb.states,1:nb.states] <- diag(1, nb.states)
  
  Id <- diag(1, nrow(W.up))
  
  lambda.minus <- 0 * (- t(H_matrix) %*% solve(M.pred[[time]]) %*% y_residual[,time])
  Lambda.minus <- 0 * (t(H_matrix) %*% solve(M.pred[[time]]) %*% H_matrix)
  
  
  # Looop
  for (i in time:1){
    
    invert_matrix <- solve(M.pred[[i]])
    
    # Updates for variance covariance matrix
    C_matrix      <- Id - Kalman_Gain[[i]] %*% H_matrix
    Phi_matrix    <- Jacobians[[i]]
    
    lambda.plus   <- t(Phi_matrix) %*% lambda.minus
    Lambda.plus   <- t(Phi_matrix) %*% Lambda.minus %*% Phi_matrix
    
    W.smooth[, i] <- W.up[,i] - P.up[[i]] %*% lambda.plus
    P.smooth[[i]] <- P.up[[i]] - P.up[[i]] %*% Lambda.plus %*% P.up[[i]]
    
    lambda.minus  <- t(C_matrix) %*% lambda.plus - t(H_matrix) %*% invert_matrix %*% y_residual[,i]
    Lambda.minus  <- t(C_matrix) %*% Lambda.plus %*% C_matrix + t(H_matrix) %*% invert_matrix %*% H_matrix
    
    
    
    # Controlling negatives
    W.smooth[W.smooth[,i] < 0, i] <- 0
    
    
    #print(i)
  }
  
  return(list(
    "W.smooth" = W.smooth,
    "P.smooth" = P.smooth
  ))
  
}









