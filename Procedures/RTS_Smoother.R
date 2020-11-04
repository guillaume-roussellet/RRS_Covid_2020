#########################################
# SMOOTHER
#########################################

Extended_smoother = function(W.up,
                             W.pred,
                             P.up,
                             P.pred,
                             Jacobians,
                             Pop_vector){
  
  require(MASS)
  require(stringr)
  
  # Objects to retrieve
  nb.states <- length(Pop_vector)
  time      <- ncol(W.up)
  
  
  # Initialization
  W.smooth = NA * W.up
  P.smooth = list()
  
  W.smooth[, time] = W.up[, time]
  P.smooth[[time]] = P.up[[time]]
  
  
  
  # loop
  for (i in (time-1):1){
    
    matrix_to_inverse <- P.pred[[i+1]]
    inverse.matrix    <- ginv(matrix_to_inverse)
    
    
    Gain.smooth = P.up[[i]] %*% t(Jacobians[[i+1]]) %*% inverse.matrix
    
    W.smooth[, i] = W.up[, i] + Gain.smooth %*% (W.smooth[, i+1] - W.pred[, i+1])
    P.smooth[[i]] = P.up[[i]] + Gain.smooth %*% (P.smooth[[i+1]] - P.pred[[i+1]]) %*% t(Gain.smooth)

    # Controlling the conditions on the smoothed W
    # Controlling negatives
    W.smooth[W.smooth[,i] < 0, i] <- 0
    
  }
  
  return(list(
    "W.smooth" = W.smooth,
    "P.smooth" = P.smooth
  ))
  
}









