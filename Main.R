####################################
# RENNE ROUSSELLET SCHWENKLER 2020 #
####################################
# Load packages
library(seas)
library(stringr)
library(Rcpp)
library(RcppEigen)
library(lubridate)
library(rstudioapi)    

# Set path
# GETS CURRENT PATH OF THE FILE
main.path <- str_sub(
  rstudioapi::getActiveDocumentContext()$path, 
  end = nchar(rstudioapi::getActiveDocumentContext()$path) - 6)
setwd(main.path)

# Load Data and functions
source("Loading_data/loadJHU.R")

Rcpp::sourceCpp('Procedures/filter_EKF_cpp.cpp')
Rcpp::sourceCpp('Procedures/filter_EKF_cpp_counterfactual.cpp')
source("Procedures/Estimation_functions.R")
source("Procedures/RTS_Smoother.R")
source("Procedures/BF_Smoother.R")

# vector of populations:
Pop_vector  <- population$population2019
nb_states   <- length(Pop_vector)
total_dates <- DEATH[,1]
nb_dates    <- length(total_dates)


########################################################
# Smoothing of Deaths
# matrix of death numbers:
Death_diff_data                       <- as.matrix(DEATH[,2:(nb.states+1)])
Death_diff_data[Death_diff_data < 0]  <- 0
D.smoothed                            <- Death_diff_data
Non.Zero.vector                       <- rep(NA, nb.states)

# De-seasonalizing of deaths:
for(i in 1:nb.states){
  D.i       <- D.smoothed[,i]
  
  # Constraining to two weeks before the first declared death
  Non.Zero.vector[i]  <- min(which(D.i != 0)) - 14
  non.zero            <- Non.Zero.vector[i]
  D.i                 <- D.i[non.zero:length(D.i)]
  
  # Transforming in time series and detrending
  D.i.ts      <- ts(c(D.i), frequency = 7)
  res.sa.D.i  <- stl(D.i.ts, s.window = 7)
  
  weekly.seas <- res.sa.D.i$time.series[,1]
  D.i         <- D.i - weekly.seas
  D.i[D.i<0]  <- 0
  
  # Storing
  D.smoothed[non.zero:nrow(D.smoothed),i] <- D.i
  
}
##############################################################
#=============================================================
#=============================================================
# Restraining to potentially a subset of states
#=============================================================
Selec.state   <- which(vector.of.states %in% vector.of.states) 

Death.selec   <- D.smoothed[, Selec.state]
if(length(Selec.state) > 1){
  Death.selec <- matrix(Death.selec[min(Non.Zero.vector[Selec.state]):nrow(Death.selec),], ncol = length(Selec.state))
} else{
  Death.selec <- matrix(Death.selec[Non.Zero.vector[Selec.state]:length(Death.selec)], ncol = 1)
}

Dates.selec   <- total_dates[min(Non.Zero.vector[Selec.state]):length(total_dates)]

nb_dates      <- nrow(Death.selec)
nb_states     <- length(Selec.state)
Pop_vector    <- population$population2019[Selec.state]

New.non.zero  <- Non.Zero.vector[Selec.state] + 1 - min(Non.Zero.vector[Selec.state])

#==============================================================
#==============================================================
# Migration and transportation matrices
Weight_migration  <- DATA.travel 
Weight_commuting  <- Matrix.commuting 
Weight_migration  <- Weight_migration[Selec.state, Selec.state]
Weight_commuting  <- Weight_commuting[Selec.state, Selec.state]

# Policy matrices
Dummy_theta_F <- as.matrix(1 - DATA.travelban[which(DATA.travelban[,1] %in% Dates.selec), (Selec.state + 1)])   # Dummy = 1 <=> high state
Dummy_theta_I <- as.matrix(1 - DATA.stayathome[which(DATA.stayathome[,1] %in% Dates.selec), (Selec.state + 1)]) # Dummy = 1 <=> high state
Dummy_masks   <- as.matrix(1 - DATA.masks[which(DATA.masks[,1] %in% Dates.selec), (Selec.state + 1)]) # Dummy = 1 <=> high state

# Completing the missing rows if we do not have data until the last date.
if(nrow(Dummy_theta_F) < nb_dates){
  Dummy_theta_F <- rbind(Dummy_theta_F, 
                         matrix(1, nb_dates - nrow(Dummy_theta_F), ncol = 1) %x% matrix(Dummy_theta_F[nrow(Dummy_theta_F),], nrow = 1))  
}
if(nrow(Dummy_theta_I) < nb_dates){
  Dummy_theta_I <- rbind(Dummy_theta_I, 
                         matrix(1, nb_dates - nrow(Dummy_theta_I), ncol = 1) %x% matrix(Dummy_theta_I[nrow(Dummy_theta_I),], nrow = 1))  
}
if(nrow(Dummy_masks) < nb_dates){
  Dummy_masks <- rbind(Dummy_masks, 
                         matrix(1, nb_dates - nrow(Dummy_masks), ncol = 1) %x% matrix(Dummy_masks[nrow(Dummy_masks),], nrow = 1))  
}

#--------------------------------------------
# RESTRICTING THE DATES
#--------------------------------------------
last.date     <- which(Dates.selec == "2020-09-30")

Dates.selec   <- Dates.selec[1:last.date]
Death.selec   <- Death.selec[1:last.date,]
Dummy_theta_F <- Dummy_theta_F[1:last.date,]
Dummy_theta_I <- Dummy_theta_I[1:last.date,]
Dummy_masks   <- Dummy_masks[1:last.date,]


######################################################################################################
######################################################################################################
#
# ESTIMATION
#
######################################################################################################
######################################################################################################
theta.0 <- c(-1.63164831,
             -2.28910817,
             -0.05679141) 



Loglik(theta.0, 
       Death.selec,
       Dummy_theta_F,
       Dummy_theta_I,
       Dummy_masks,
       Weight_commuting,
       Weight_migration,
       Pop_vector)

library(optimx)

Optimisation <- optimx(par = theta.0, fn = Loglik,
                       Death.selec = Death.selec,
                       Dummy_theta_F = Dummy_theta_F,
                       Dummy_theta_I = Dummy_theta_I,
                       Dummy_masks = Dummy_masks,
                       Weight_commuting = Weight_commuting,
                       Weight_migration = Weight_migration,
                       Pop_vector = Pop_vector,
                       method = "nlminb",
                       itnmax = 100,
                       hessian = F, 
                       control = list(trace = 1,
                                      kkt = F,
                                      starttests = F,
                                      maxit = 100))
theta.0 <- as.numeric(coef(Optimisation))
  
# Computation of Estimation Results
Filtrage_EKF <- Loglik(theta.0, 
                      Death.selec,
                      Dummy_theta_F,
                      Dummy_theta_I,
                      Dummy_masks,
                      Weight_commuting,
                      Weight_migration,
                      Pop_vector, for_estim = F)

# Print the total likelihood
cat("Total LogLik: \n",sum(Filtrage_EKF$loglik.vector[1:length(Filtrage_EKF$loglik.vector)]))

# Printing the Parameters estimated through MaxLik
theta.to.param(theta.0,2)
Parameters <- theta.to.param(theta.0,nb_states)

######################################################################################################
######################################################################################################
#
# PLOTS AND SMOOTHING
#
######################################################################################################
######################################################################################################
WINDOWS <- c(6,9) 


Smoothed_data <- Extended_smoother(Filtrage_EKF$W.updated,
                                   Filtrage_EKF$W.predicted,
                                   Filtrage_EKF$P.updated,
                                   Filtrage_EKF$P.predicted,
                                   Filtrage_EKF$Grads, 
                                   Pop_vector)

lower.bounds.KF <- Filtrage_EKF$W.updated - 
  2 * matrix(
    unlist(
      lapply(
        Filtrage_EKF$P.updated, function(x){sqrt(diag(x))})), 
    nrow(Filtrage_EKF$W.updated), ncol(Filtrage_EKF$W.updated))
lower.bounds.KF[lower.bounds.KF < 0] <- 0
upper.bounds.KF <- Filtrage_EKF$W.updated + 2 * matrix(
  unlist(
    lapply(
      Filtrage_EKF$P.updated, function(x){sqrt(diag(x))})), 
  nrow(Filtrage_EKF$W.updated), ncol(Filtrage_EKF$W.updated))

lower.bounds.KS <- Smoothed_data$W.smooth - 2 * matrix(
  unlist(
    lapply(
      Smoothed_data$P.smooth, function(x){sqrt(diag(x))})), 
  nrow(Filtrage_EKF$W.updated), ncol(Filtrage_EKF$W.updated))
lower.bounds.KS[lower.bounds.KS < 0] <- 0
upper.bounds.KS <- Smoothed_data$W.smooth + 2 * matrix(
  unlist(
    lapply(
      Smoothed_data$P.smooth, function(x){sqrt(diag(x))})), 
  nrow(Filtrage_EKF$W.updated), ncol(Filtrage_EKF$W.updated))


# Bounds for Betas
lower.bounds.Betas.KF <- Filtrage_EKF$Betas.updated - 2 * matrix(
  unlist(
    lapply(
      Filtrage_EKF$P.updated, function(x){sqrt(diag(
        Parameters$Premult_Beta %*% x[(4 * nb_states + 1):(5 * nb_states),(4 * nb_states + 1):(5 * nb_states)] %*% Parameters$Premult_Beta
      ))})), nb_states, ncol(Filtrage_EKF$Betas.updated))
lower.bounds.Betas.KF[lower.bounds.Betas.KF<0] <- 0
upper.bounds.Betas.KF <- Filtrage_EKF$Betas.updated + 2 * matrix(
  unlist(
    lapply(
      Filtrage_EKF$P.updated, function(x){sqrt(diag(
        Parameters$Premult_Beta %*% x[(4 * nb_states + 1):(5 * nb_states),(4 * nb_states + 1):(5 * nb_states)] %*% Parameters$Premult_Beta
      ))})), nb_states, ncol(Filtrage_EKF$Betas.updated))

lower.bounds.Betas.KS <- Parameters$Premult_Beta %*% Smoothed_data$W.smooth[(4 * nb_states + 1):(5 * nb_states),] - 2 * matrix(
  unlist(
    lapply(
      Smoothed_data$P.smooth, function(x){sqrt(diag(
        Parameters$Premult_Beta %*% x[(4 * nb_states + 1):(5 * nb_states),(4 * nb_states + 1):(5 * nb_states)] %*% Parameters$Premult_Beta
      ))})), nb_states, ncol(Filtrage_EKF$Betas.updated))
lower.bounds.Betas.KS[lower.bounds.Betas.KS<0] <- 0
upper.bounds.Betas.KS <- Parameters$Premult_Beta %*% Smoothed_data$W.smooth[(4 * nb_states + 1):(5 * nb_states),] + 2 * matrix(
  unlist(
    lapply(
      Smoothed_data$P.smooth, function(x){sqrt(diag(
        Parameters$Premult_Beta %*% x[(4 * nb_states + 1):(5 * nb_states),(4 * nb_states + 1):(5 * nb_states)] %*% Parameters$Premult_Beta
      ))})), nb_states, ncol(Filtrage_EKF$Betas.updated))


########################################################
# GET PLOTS
########################################################
# Death plot
order.to.plot.states <- sort(ordered.state.names, index.return = T)$ix

#--------------------------------------------------------------------
par(mfrow = WINDOWS, mai = c(0.3,0.3,0.3,0.1))
for(i in order.to.plot.states){
  plot(Dates.selec, Filtrage_EKF$W.updated[i,], type='l', 
       main = ordered.state.names[Selec.state[i]], cex.main = 1,
       panel.first = grid(),
       xlab = "", ylab = "per day", 
       col = "dark grey",
       ylim=c(min(lower.bounds.KF[i,],lower.bounds.KS[i,],
                  Smoothed_data$W.smooth[i,],na.rm = TRUE),
              max(upper.bounds.KF[i,],upper.bounds.KS[i,],
                  Smoothed_data$W.smooth[i,],na.rm = TRUE)))
  
  
  # Conf Bands
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KF[i,], rev(upper.bounds.KF[i,])),
          col="#44444433",border = NaN)
  
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KS[i,], rev(upper.bounds.KS[i,])),
          col="#99444455",border = NaN)
  
  lines(Dates.selec, Smoothed_data$W.smooth[i,], col = 'red', lwd = 2)
  lines(Dates.selec, Death.selec[,i], col = "blue")
}
plot(Dates.selec, apply(Filtrage_EKF$W.updated[1:nb.states,], 2,sum), panel.first = grid(),
     xlab = "", ylab = "per day", 
     col = "dark grey", type = 'l', main = "U.S. Total")
# Conf Bands
polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KF[1:nb.states,],2,sum), rev(apply(upper.bounds.KF[1:nb.states,],2,sum))),
        col="#44444433",border = NaN)

polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KS[1:nb.states,],2,sum), rev(apply(upper.bounds.KS[1:nb.states,],2,sum))),
        col="#99444455",border = NaN)
lines(Dates.selec, apply(Smoothed_data$W.smooth[1:nb.states,],2,sum), col = 'red', lwd = 2)
lines(Dates.selec, apply(Death.selec,1,sum), col = "blue")

plot(c(0:1), c(0:1) , type = "n", xaxt = "n", yaxt = "n", main = "", ylab = "", xlab = "", bty = 'n')
legend("center",
  c("data", "filter", "smoother"),
 col = c("blue", "dark grey", "red" ), 
 lty = c(1,1,1), bty = "n")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Susceptible plot
#--------------------------------------------------------------------
par(mfrow = WINDOWS, mai = c(0.3,0.3,0.3,0.1))
for(i in order.to.plot.states){
  
  index.to.plot <- nb_states + i
  
  plot(Dates.selec, Filtrage_EKF$W.updated[index.to.plot,], type='l', 
       main = ordered.state.names[Selec.state[i]], cex.main = 1,
       panel.first = grid(),
       col = "dark grey",
       xlab = "", ylab = "total", 
       ylim=c(min(lower.bounds.KF[index.to.plot,],lower.bounds.KS[index.to.plot,],
                  Smoothed_data$W.smooth[index.to.plot,],na.rm = TRUE),
              max(upper.bounds.KF[index.to.plot,],upper.bounds.KS[index.to.plot,],
                  Smoothed_data$W.smooth[index.to.plot,],na.rm = TRUE)),
       lwd = 2)
  
  polygon(c(Dates.selec[(Dummy_theta_F[,Selec.state[i]]==0)],rev(Dates.selec[(Dummy_theta_F[,Selec.state[i]]==0)])), 
          c(rep(1e9, sum(Dummy_theta_F[,Selec.state[i]]==0)), rep(0, sum(Dummy_theta_F[,Selec.state[i]]==0))),
          col=rgb(44/255, 236/255,17/255, alpha = .33),border = NaN)
  # Conf Bands
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KF[index.to.plot,], rev(upper.bounds.KF[index.to.plot,])),
          col="#44444433",border = NaN)
  
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KS[index.to.plot,], rev(upper.bounds.KS[index.to.plot,])),
          col="#99444455",border = NaN)
  
  lines(Dates.selec, Smoothed_data$W.smooth[index.to.plot,], col = 'red', lwd = 2)
  
  
 
}
plot(Dates.selec, apply(Filtrage_EKF$W.updated[(nb.states + 1):(2 * nb.states),],2,sum), type='l', 
     main = "U.S. Total",
     panel.first = grid(),
     col = "dark grey",
     xlab = "", ylab = "total", 
     ylim=c(min(apply(lower.bounds.KF[(nb.states + 1):(2 * nb.states),],2,sum),
                apply(lower.bounds.KS[(nb.states + 1):(2 * nb.states),],2,sum), na.rm = TRUE),
            max(apply(upper.bounds.KF[(nb.states + 1):(2 * nb.states),],2,sum),
                apply(upper.bounds.KS[(nb.states + 1):(2 * nb.states),],2,sum), na.rm = TRUE)),
     lwd = 2)
# Conf Bands
polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KF[(nb.states + 1):(2 * nb.states),],2,sum), rev(apply(upper.bounds.KF[(nb.states + 1):(2 * nb.states),],2,sum))),
        col="#44444433",border = NaN)

polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KS[(nb.states + 1):(2 * nb.states),],2,sum), rev(apply(upper.bounds.KS[(nb.states + 1):(2 * nb.states),],2,sum))),
        col="#99444455",border = NaN)

lines(Dates.selec, apply(Smoothed_data$W.smooth[(nb.states + 1):(2 * nb.states),],2,sum), col = 'red', lwd = 2)
plot(c(0:1), c(0:1) , type = "n", xaxt = "n", yaxt = "n", main = "", ylab = "", xlab = "", bty = 'n')
legend("center",
       c("filter", "smoother", "travelban"),
       col = c("dark grey", "red" ,rgb(44/255, 236/255,17/255, alpha = .33)), 
       lty = c(1,1,1), bty = "n")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Infected plot
#--------------------------------------------------------------------
par(mfrow = WINDOWS, mai = c(0.3,0.3,0.3,0.1))
for(i in order.to.plot.states){
  
  index.to.plot <- 2 * nb_states + i
  
  plot(Dates.selec, Filtrage_EKF$W.updated[index.to.plot,], type='l', 
       main = ordered.state.names[Selec.state[i]], cex.main = 1,
       panel.first = grid(),
       col = "dark grey",
       xlab = "", ylab = "total", 
       ylim=c(min(lower.bounds.KF[index.to.plot,],lower.bounds.KS[index.to.plot,],
                  Smoothed_data$W.smooth[index.to.plot,],na.rm = TRUE),
              max(upper.bounds.KF[index.to.plot,],upper.bounds.KS[index.to.plot,],
                  Smoothed_data$W.smooth[index.to.plot,],na.rm = TRUE)), lwd = 2)
  
  
  # Conf Bands
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KF[index.to.plot,], rev(upper.bounds.KF[index.to.plot,])),
          col="#44444433",border = NaN)
  
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KS[index.to.plot,], rev(upper.bounds.KS[index.to.plot,])),
          col="#99444455",border = NaN)
  
  lines(Dates.selec, Smoothed_data$W.smooth[index.to.plot,], col = 'red', lwd = 2)
  
}
plot(Dates.selec, apply(Filtrage_EKF$W.updated[(2 * nb_states + 1):(3 * nb_states),],2, sum), type='l', 
     main = "U.S. Total",
     panel.first = grid(),
     col = "dark grey",
     xlab = "", ylab = "total", 
     ylim=c(min(apply(lower.bounds.KF[(2 * nb_states + 1):(3 * nb_states),],2, sum),
                apply(lower.bounds.KS[(2 * nb_states + 1):(3 * nb_states),],2, sum),na.rm = TRUE),
            max(apply(upper.bounds.KF[(2 * nb_states + 1):(3 * nb_states),],2, sum),
                apply(upper.bounds.KS[(2 * nb_states + 1):(3 * nb_states),],2, sum),na.rm = TRUE)), lwd = 2)


# Conf Bands
polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KF[(2 * nb_states + 1):(3 * nb_states),],2, sum), rev(apply(upper.bounds.KF[(2 * nb_states + 1):(3 * nb_states),],2, sum))),
        col="#44444433",border = NaN)

polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KS[(2 * nb_states + 1):(3 * nb_states),],2, sum), rev(apply(upper.bounds.KS[(2 * nb_states + 1):(3 * nb_states),],2, sum))),
        col="#99444455",border = NaN)

lines(Dates.selec, apply(Smoothed_data$W.smooth[(2 * nb_states + 1):(3 * nb_states),],2, sum), col = 'red', lwd = 2)

plot(c(0:1), c(0:1) , type = "n", xaxt = "n", yaxt = "n", main = "", ylab = "", xlab = "", bty = 'n')
legend("center",
       c("filter", "smoother"),
       col = c("dark grey", "red" ), 
       lty = c(1,1), bty = "n")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Recovered plot
#--------------------------------------------------------------------
par(mfrow = WINDOWS, mai = c(0.3,0.3,0.3,0.1))
for(i in order.to.plot.states){
  
  index.to.plot <- 3 * nb_states + i
  
  plot(Dates.selec, Filtrage_EKF$W.updated[index.to.plot,], type='l', 
       main = ordered.state.names[Selec.state[i]], cex.main = 1,
       panel.first = grid(),
       col = "dark grey",
       xlab = "", ylab = "total", 
       ylim=c(min(lower.bounds.KF[index.to.plot,],lower.bounds.KS[index.to.plot,],
                  Smoothed_data$W.smooth[index.to.plot,],na.rm = TRUE),
              max(upper.bounds.KF[index.to.plot,],upper.bounds.KS[index.to.plot,],
                  Smoothed_data$W.smooth[index.to.plot,],na.rm = TRUE)), lwd = 2)
  
  
  # Conf Bands
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KF[index.to.plot,], rev(upper.bounds.KF[index.to.plot,])),
          col="#44444433",border = NaN)
  
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(lower.bounds.KS[index.to.plot,], rev(upper.bounds.KS[index.to.plot,])),
          col="#99444455",border = NaN)
  
  lines(Dates.selec, Smoothed_data$W.smooth[index.to.plot,], col = 'red', lwd = 2)
  lines(total_dates, RECOV[,str_c("recovered.",vector.of.states[Selec.state[i]])], col = "blue", lwd = 2)
}
plot(Dates.selec, apply(Filtrage_EKF$W.updated[(3 * nb_states + 1):(4 * nb_states),], 2, sum), type='l', 
     main = "U.S. Total", cex.main = 1,
     panel.first = grid(),
     col = "dark grey",
     xlab = "", ylab = "total", 
     ylim=c(min(apply(lower.bounds.KF[(3 * nb_states + 1):(4 * nb_states),], 2, sum),
                apply(lower.bounds.KS[(3 * nb_states + 1):(4 * nb_states),], 2, sum),na.rm = TRUE),
            max(apply(upper.bounds.KF[(3 * nb_states + 1):(4 * nb_states),], 2, sum),
                apply(upper.bounds.KS[(3 * nb_states + 1):(4 * nb_states),], 2, sum),na.rm = TRUE)), lwd = 2)


# Conf Bands
polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KF[(3 * nb_states + 1):(4 * nb_states),], 2, sum), rev(apply(upper.bounds.KF[(3 * nb_states + 1):(4 * nb_states),], 2, sum))),
        col="#44444433", border = NaN)

polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(lower.bounds.KS[(3 * nb_states + 1):(4 * nb_states),], 2, sum), rev(apply(upper.bounds.KS[(3 * nb_states + 1):(4 * nb_states),], 2, sum))),
        col="#99444455", border = NaN)

lines(Dates.selec, apply(Smoothed_data$W.smooth[(3 * nb_states + 1):(4 * nb_states),], 2, sum), col = 'red', lwd = 2)
lines(total_dates, apply(RECOV[,str_c("recovered.",vector.of.states[Selec.state])], 1, 
                         function(x){sum(x, na.rm = T)}), col = "blue", lwd = 2)

plot(c(0:1), c(0:1) , type = "n", xaxt = "n", yaxt = "n", main = "", ylab = "", xlab = "", bty = 'n')
legend("center",
       c("filter", "smoother", "data"),
       col = c("dark grey", "red" , "blue"), 
       lty = c(1,1), bty = "n")


#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Betas plot
#--------------------------------------------------------------------
par(mfrow = WINDOWS, mai = c(0.3,0.3,0.3,0.1))
Smoothed_Betas <- Parameters$Premult_Beta %*% Smoothed_data$W.smooth[(4 * nb_states + 1):(5 * nb_states),]
premult_theta_I_total <- NULL
for(i in order.to.plot.states){
  
  index.to.plot   <- 4 * nb_states + i
  premult_theta_I <- (Parameters$Mask_reduc + Dummy_masks[,Selec.state[i]] * (1-Parameters$Mask_reduc)) * (
    Parameters$theta_I_low[Selec.state[i]] + 
      Dummy_theta_I[,Selec.state[i]] *(
        Parameters$theta_I_high[Selec.state[i]] - Parameters$theta_I_low[Selec.state[i]]
        )
  )
  premult_theta_I_total <- cbind(premult_theta_I_total, premult_theta_I)
  
  plot(Dates.selec, premult_theta_I * Filtrage_EKF$Betas.updated[i,], type='l', 
       main = ordered.state.names[Selec.state[i]], cex.main = 1,
       panel.first = grid(),
       col = "dark grey",
       xlab = "", ylab = "total", 
       ylim=c(min(premult_theta_I * lower.bounds.Betas.KF[i,],
                  premult_theta_I * lower.bounds.Betas.KS[i,],
                  premult_theta_I * Smoothed_Betas[i,], na.rm = TRUE),
              min(max(premult_theta_I * upper.bounds.Betas.KF[i,],
                      premult_theta_I * upper.bounds.Betas.KS[i,],
                      premult_theta_I * Smoothed_Betas[i,],na.rm = TRUE), 1)), lwd = 2)
  
  polygon(c(Dates.selec[(Dummy_masks[,Selec.state[i]]==0)],rev(Dates.selec[(Dummy_masks[,Selec.state[i]]==0)])), 
          c(rep(1e5, sum(Dummy_masks[,Selec.state[i]]==0)), rep(0, sum(Dummy_masks[,Selec.state[i]]==0))),
          col=rgb(139/255, 34/255, 139/255, alpha = .33),border = NaN)
  polygon(c(Dates.selec[(Dummy_theta_I[,Selec.state[i]]==0)],rev(Dates.selec[(Dummy_theta_I[,Selec.state[i]]==0)])), 
          c(rep(1e5, sum(Dummy_theta_I[,Selec.state[i]]==0)), rep(0, sum(Dummy_theta_I[,Selec.state[i]]==0))),
          col=rgb(44/255, 236/255,17/255, alpha = .33),border = NaN)
  # Conf Bands
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(premult_theta_I * lower.bounds.Betas.KF[i,], 
            rev(premult_theta_I * upper.bounds.Betas.KF[i,])),
          col="#44444433",border = NaN)
  
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(premult_theta_I * lower.bounds.Betas.KS[i,], 
            rev(premult_theta_I * upper.bounds.Betas.KS[i,])),
          col="#99444455",border = NaN)
  
  lines(Dates.selec, premult_theta_I * Smoothed_Betas[i,], col = 'red', lwd = 2)
  
}

plot(Dates.selec, apply(
  diag(Pop_vector) %*% (t(premult_theta_I_total) * Filtrage_EKF$Betas.updated)/sum(Pop_vector),2,sum), 
     type='l', 
     main = "U.S. Total", cex.main = 1,
     panel.first = grid(),
     col = "dark grey",
     xlab = "", ylab = "total", 
     ylim=c(min(apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KF)/
                        sum(Pop_vector),2,sum),
                apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KS)/
                        sum(Pop_vector),2,sum), 
                na.rm = TRUE),
            max(apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KF)/
                        sum(Pop_vector),2,sum),
                apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KS)/
                        sum(Pop_vector),2,sum), 
                na.rm = TRUE)), lwd = 2)
# Conf Bands
polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(
          diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KF)/sum(Pop_vector),2,sum), 
          rev(apply(
            diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KF)/sum(Pop_vector),2,sum))),
        col="#44444433",border = NaN)

polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(
          diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KS)/sum(Pop_vector),2,sum), 
          rev(apply(
            diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KS)/sum(Pop_vector),2,sum))),
        col="#99444455",border = NaN)

lines(Dates.selec, apply(
  diag(Pop_vector) %*% (t(premult_theta_I_total) * Smoothed_Betas)/sum(Pop_vector),2,sum), 
  col = 'red', lwd = 2)

plot(c(0:1), c(0:1) , type = "n", xaxt = "n", yaxt = "n", main = "", ylab = "", xlab = "", bty = 'n')
legend("center",
       c("filter", "smoother", "stay-at-home", "masks"),
       col = c("dark grey", "red" , rgb(44/255, 236/255,17/255, alpha = .33),
               rgb(139/255, 34/255, 139/255, alpha = .33)), 
       lty = c(1,1,1,1), lwd = c(2,2,2,2), bty = "n")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Plot Infected once
#--------------------------------------------------------------------
par(mfrow = WINDOWS, mai = c(0.3,0.3,0.3,0.1))
for(i in order.to.plot.states){
  
  index.to.plot.I <- 2 * nb_states + i
  index.to.plot.R <- 3 * nb_states + i
  index.to.plot.D <- 0 * nb_states + i
  
  plot(Dates.selec, Filtrage_EKF$W.updated[index.to.plot.I,] + 
         Filtrage_EKF$W.updated[index.to.plot.R,]+ 
         Filtrage_EKF$W.updated[index.to.plot.D,], 
       type='l', main = ordered.state.names[Selec.state[i]], cex.main = 1,
       panel.first = grid(),
       col = "dark grey",
       xlab = "", ylab = "total")

  
  lines(Dates.selec, Smoothed_data$W.smooth[index.to.plot.I,] + 
          Smoothed_data$W.smooth[index.to.plot.R,]+ 
          Smoothed_data$W.smooth[index.to.plot.D,], col = 'red', lwd = 2)
  lines(total_dates, INFEC[,str_c("positiveCasesViral.",vector.of.states[Selec.state[i]])], col = "blue", lwd = 2)
}
plot(Dates.selec, apply(Filtrage_EKF$W.updated[(2 * nb_states + 1):(3 * nb_states),], 2, sum) + 
       apply(Filtrage_EKF$W.updated[(3 * nb_states + 1):(4 * nb_states),], 2, sum)+ 
       apply(Filtrage_EKF$W.updated[(0 * nb_states + 1):(1 * nb_states),], 2, sum), 
     type='l', main = "U.S. Total", cex.main = 1,
     panel.first = grid(),
     col = "dark grey",
     xlab = "", ylab = "total")


lines(Dates.selec, apply(Smoothed_data$W.smooth[(2 * nb_states + 1):(3 * nb_states),], 2, sum) + 
        apply(Smoothed_data$W.smooth[(3 * nb_states + 1):(4 * nb_states),], 2, sum)+ 
        apply(Smoothed_data$W.smooth[(0 * nb_states + 1):(1 * nb_states),], 2, sum), col = 'red', lwd = 2)
lines(total_dates, apply(INFEC[,str_c("positiveCasesViral.",vector.of.states[Selec.state])],1,
                         function(x){sum(x, na.rm=T)}), col = "blue", lwd = 2)

plot(c(0:1), c(0:1) , type = "n", xaxt = "n", yaxt = "n", main = "", ylab = "", xlab = "", bty = "n")
legend("center",
       c("filter", "smoother", 'data'),
       col = c("dark grey", "red" , 'blue'), 
       lty = c(1,1), bty = "n")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# R0 plot
#--------------------------------------------------------------------
par(mfrow = WINDOWS, mai = c(0.3,0.3,0.3,0.1))
Smoothed_Betas <- Parameters$Premult_Beta %*% Smoothed_data$W.smooth[(4 * nb_states + 1):(5 * nb_states),]
for(i in order.to.plot.states){
  
  index.to.plot <- 4 * nb_states + i
  premult_theta_I <- (Parameters$Mask_reduc + Dummy_masks[,Selec.state[i]] * (1-Parameters$Mask_reduc)) * (
    Parameters$theta_I_low[Selec.state[i]] + 
      Dummy_theta_I[,Selec.state[i]] *(
        Parameters$theta_I_high[Selec.state[i]] - Parameters$theta_I_low[Selec.state[i]]
      )
  )
  
  plot(Dates.selec, premult_theta_I * Filtrage_EKF$Betas.updated[i,]/
         (Parameters$delta + Parameters$gamma), type='l', 
       main = ordered.state.names[Selec.state[i]], cex.main = 1,
       panel.first = grid(),
       col = "dark grey",
       xlab = "", ylab = "total", 
       ylim=c(0,
              min(max(premult_theta_I * upper.bounds.Betas.KS[i,],
                      premult_theta_I * Smoothed_Betas[i,],na.rm = TRUE), 1))/(Parameters$delta + Parameters$gamma), 
       lwd = 2)
  
  polygon(c(Dates.selec[(Dummy_masks[,Selec.state[i]]==0)],rev(Dates.selec[(Dummy_masks[,Selec.state[i]]==0)])), 
          c(rep(1e5, sum(Dummy_masks[,Selec.state[i]]==0)), rep(0, sum(Dummy_masks[,Selec.state[i]]==0))),
          col=rgb(139/255, 34/255, 139/255, alpha = .33),border = NaN)
  polygon(c(Dates.selec[(Dummy_theta_I[,Selec.state[i]]==0)],
            rev(Dates.selec[(Dummy_theta_I[,Selec.state[i]]==0)])), 
          c(rep(1e5, sum(Dummy_theta_I[,Selec.state[i]]==0)), 
            rep(0, sum(Dummy_theta_I[,Selec.state[i]]==0)))/(Parameters$delta + Parameters$gamma),
          col=rgb(44/255, 236/255,17/255, alpha = .33),border = NaN)
  # Conf Bands
  polygon(c(Dates.selec, rev(Dates.selec)),
          c(premult_theta_I * lower.bounds.Betas.KF[i,], 
            rev(premult_theta_I * upper.bounds.Betas.KF[i,]))/(Parameters$delta + Parameters$gamma),
          col="#44444433",border = NaN)
  
  polygon(c(Dates.selec,rev(Dates.selec)),
          c(premult_theta_I * lower.bounds.Betas.KS[i,], 
            rev(premult_theta_I * upper.bounds.Betas.KS[i,]))/(Parameters$delta + Parameters$gamma),
          col="#99444455",border = NaN)
  
  lines(Dates.selec, premult_theta_I * Smoothed_Betas[i,]/
          (Parameters$delta + Parameters$gamma), col = 'red', lwd = 2)
  
  abline(h=1, lwd = 2, col = 'blue')
}

plot(Dates.selec, apply(
  diag(Pop_vector) %*% (t(premult_theta_I_total) * Filtrage_EKF$Betas.updated)/sum(Pop_vector),2,sum)/
    (Parameters$delta + Parameters$gamma), 
  type='l', 
  main = "U.S. Total", cex.main = 1,
  panel.first = grid(),
  col = "dark grey",
  xlab = "", ylab = "total", 
  ylim=c(min(apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KF)/sum(Pop_vector),2,sum),
             apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KS)/sum(Pop_vector),2,sum), 
             na.rm = TRUE),
         max(apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KF)/sum(Pop_vector),2,sum),
             apply(diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KS)/sum(Pop_vector),2,sum), 
             na.rm = TRUE))/
    (Parameters$delta + Parameters$gamma), lwd = 2)
# Conf Bands
polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(
          diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KF)/sum(Pop_vector),2,sum), 
          rev(apply(
            diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KF)/sum(Pop_vector),2,sum)))/
          (Parameters$delta + Parameters$gamma),
        col="#44444433",border = NaN)

polygon(c(Dates.selec,rev(Dates.selec)),
        c(apply(
          diag(Pop_vector) %*% (t(premult_theta_I_total) * lower.bounds.Betas.KS)/sum(Pop_vector),2,sum), 
          rev(apply(
            diag(Pop_vector) %*% (t(premult_theta_I_total) * upper.bounds.Betas.KS)/sum(Pop_vector),2,sum)))/
          (Parameters$delta + Parameters$gamma),
        col="#99444455",border = NaN)

lines(Dates.selec, apply(
  diag(Pop_vector) %*% (t(premult_theta_I_total) * Smoothed_Betas)/sum(Pop_vector),2,sum)/
    (Parameters$delta + Parameters$gamma), 
  col = 'red', lwd = 2)
abline(h=1, lwd = 2, col = 'blue')
plot(c(0:1), c(0:1) , type = "n", xaxt = "n", yaxt = "n", main = "", ylab = "", xlab = "", bty = 'n')
legend("center",
       c("filter", "smoother", "stay-at-home", "masks"),
       col = c("dark grey", "red" , rgb(44/255, 236/255,17/255, alpha = .33),
               rgb(139/255, 34/255, 139/255, alpha = .33)), 
       lty = c(1,1,1,1), lwd = c(2,2,2,2), bty = "n")

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# COMPUTING AGGREGATE NUMBER OF INFECTED ONCE
#--------------------------------------------------------------------
index.to.plot.I <- (2 * nb_states + 1):(3 * nb_states)
index.to.plot.R <- (3 * nb_states + 1):(4 * nb_states)
index.to.plot.D <- 1:nb_states

par(mfrow = c(1,1), mai = c(0.5,0.5,0.3,0.1))
Filtered_infec_across_US <- apply(Filtrage_EKF$W.updated[index.to.plot.I,] + 
                                    Filtrage_EKF$W.updated[index.to.plot.R,]+ 
                                    Filtrage_EKF$W.updated[index.to.plot.D,], 2,sum)
Smoothed_infed_across_US <- apply(Smoothed_data$W.smooth[index.to.plot.I,] + 
                                    Smoothed_data$W.smooth[index.to.plot.R,]+ 
                                    Smoothed_data$W.smooth[index.to.plot.D,], 2,sum)

DATA_infec_across_US <- apply(INFEC[,2:52],1,function(x){sum(x, na.rm = T)})

plot(Dates.selec, Filtered_infec_across_US/DATA_infec_across_US[INFEC[,1] %in% Dates.selec], 
     type='l', main = "Ratio Model Infections/Data", panel.first = grid(col = "grey60"),
     col = "dark grey", log = "y", lwd = 2,
     xlab = "", ylab = "total")
lines(Dates.selec, Smoothed_infed_across_US/DATA_infec_across_US[INFEC[,1] %in% Dates.selec], col = 'red', lwd = 2)
abline(h=7, lty = 2, lwd = 2)
legend("topright",
       c("filter", "smoother", "The Economist (7)"),
       col = c("dark grey", "red", "black"), 
       lty = c(1,1,2 ), lwd = c(2,2,2), bty = "n")
#--------------------------------------------------------------------



