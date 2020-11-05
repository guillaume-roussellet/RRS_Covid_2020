#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# FILTERING THE COUNTERFACTUALS
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
###########################################################################################
###########################################################################################
# 
# CAUTION: YOU NEED TO HAVE RUN THE MAIN BEFORE USING THIS PIECE OF CODE. 
# - THE DATA IS LOADED IN THE MAIN.
# - PARAMETERS ARE DEFINED IN THE MAIN.
# IF YOU DON'T, IT WON'T COMPILE.
if(exists("Parameters")==F){
  cat("``Main.R'' needs to be run beforehand")
}
###########################################################################################
###########################################################################################
# Stricter Dummies
#-----------------------------------------------------------------------------
# We hereby identify the dummies for the strict scenario, i.e. policy that starts
# as early as the earliest state as stops as late as the latest state.

# STRICTER TRAVEL BANS
stricter_F            <- c(min(apply(Dummy_theta_F,2,function(x){min(which(diff(x)==-1))})) + 1,
                           min(
                             max(apply(Dummy_theta_F[,apply(Dummy_theta_F,2,sum)!=nrow(Dummy_theta_F)],2,function(x){min(which(diff(x)==1))}))+1,
                             nrow(Dummy_theta_F)))
stricter_F_vec        <- rep(1,nrow(Dummy_theta_F))
stricter_F_vec[stricter_F[1]:stricter_F[2]] <- 0  
stricter_F            <- stricter_F_vec
Strict_Dummy_theta_F  <- matrix(stricter_F_vec , ncol = 1) %*% matrix(1, 1, nb_states)

# STRICTER STAY AT HOME
stricter_I            <- c(min(apply(Dummy_theta_I,2,function(x){min(which(diff(x)==-1))})) + 1,
                           min(
                             max(apply(Dummy_theta_I[,apply(Dummy_theta_I,2,sum)!=nrow(Dummy_theta_I)],2,function(x){min(which(diff(x)==1))}))+1,
                             nrow(Dummy_theta_I)))
stricter_I_vec        <- rep(1,nrow(Dummy_theta_I))
stricter_I_vec[stricter_I[1]:stricter_I[2]] <- 0  
stricter_I            <- stricter_I_vec
Strict_Dummy_theta_I  <- matrix(stricter_I_vec , ncol = 1) %*% matrix(1, 1, nb_states)

# STRICTER MASKS
stricter_masks        <- c(min(apply(Dummy_masks,2,function(x){min(which(diff(x)==-1))})) + 1,
                           min(
                             max(apply(Dummy_masks[,apply(Dummy_masks,2,sum)!=nrow(Dummy_masks)],2,function(x){min(which(diff(x)==1))}))+1,
                             nrow(Dummy_masks)))
stricter_masks_vec    <- rep(1,nrow(Dummy_theta_I))
stricter_masks_vec[stricter_masks[1]:stricter_masks[2]] <- 0  
stricter_masks        <- stricter_masks_vec
Strict_Dummy_masks    <- matrix(stricter_masks_vec , ncol = 1) %*% matrix(1, 1, nb_states)


# Looser Dummies
#-----------------------------------------------------------------------------
Loose_Dummy_theta_F  <- 0 * Dummy_theta_F + 1
Loose_Dummy_theta_I  <- 0 * Dummy_theta_I + 1
Loose_Dummy_masks    <- 0 * Dummy_masks + 1



# Other Counterfactual Dummies
#-----------------------------------------------------------------------------
Super_strict_dummy_masks    <- matrix(stricter_I_vec , ncol = 1) %*% matrix(1, 1, nb_states)
Early_dummy_travel_March15  <- matrix(Dates.selec < "2020-03-15" , ncol = 1) %*% matrix(1, 1, nb_states)
Early_dummy_travel_March1   <- matrix(Dates.selec < "2020-03-1" , ncol = 1) %*% matrix(1, 1, nb_states)
Early_dummy_travel_Feb15    <- matrix(Dates.selec < "2020-02-15" , ncol = 1) %*% matrix(1, 1, nb_states)
Early_dummy_travel_Feb12    <- matrix(Dates.selec < "2020-02-12" , ncol = 1) %*% matrix(1, 1, nb_states)

#============================================================================================
#============================================================================================
# AGGREGATED COUNTERFACTUALS OVER US
#============================================================================================
#============================================================================================
# INITIALIZING THE MEAN AND VCOV FOR THE FILTER
#----------------------------------------------
W0 <- c(rep(0, nb_states), Pop_vector, rep(10, nb_states), rep(0, nb_states), 
        rep(Parameters$start.beta, nb_states),
        rep(0, nb_states), Pop_vector, rep(10, nb_states), rep(0, nb_states))

P0 <- diag(0, 9 * nb_states)

# BLOCKS OF THE BASIC SYSTEM
# block of S
P0[(nb_states + 1):(2 * nb_states), (nb_states + 1):(2 * nb_states)]          <- diag(500, nb_states)
# block of I
P0[(2 * nb_states + 1):(3 * nb_states), (2 * nb_states + 1):(3 * nb_states)]  <- diag(500, nb_states)
# Covariance terms 
P0[(nb_states + 1):(2 * nb_states), (2 * nb_states + 1):(3 * nb_states)]      <- - 0.5 * (
  P0[(nb_states + 1):(2 * nb_states), (nb_states + 1):(2 * nb_states)] + 
    P0[(2 * nb_states + 1):(3 * nb_states), (2 * nb_states + 1):(3 * nb_states)])
P0[(2 * nb_states + 1):(3 * nb_states), (nb_states + 1):(2 * nb_states)]      <- 
  t(P0[(nb_states + 1):(2 * nb_states), (2 * nb_states + 1):(3 * nb_states)] )

# BLOCK OF THE COUNTERFACTUALS
P0[(4 * nb_states + 1):(5 * nb_states),(4 * nb_states + 1):(5 * nb_states)]   <- diag(.001^2, nb_states)
P0[(5 * nb_states + 1):(9 * nb_states),(5 * nb_states + 1):(9 * nb_states)]   <-
  P0[(0 * nb_states + 1):(4 * nb_states),(0 * nb_states + 1):(4 * nb_states)]
P0[(5 * nb_states + 1):(9 * nb_states),(0 * nb_states + 1):(4 * nb_states)]   <-
  P0[(0 * nb_states + 1):(4 * nb_states),(0 * nb_states + 1):(4 * nb_states)]
P0[(0 * nb_states + 1):(4 * nb_states),(5 * nb_states + 1):(9 * nb_states)]   <-
  P0[(0 * nb_states + 1):(4 * nb_states),(0 * nb_states + 1):(4 * nb_states)]


# Strict overall
#---------------
Filter_US_Strict_All <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                               Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                               Weight_commuting, Weight_migration, Pop_vector,
                                               Strict_Dummy_theta_F,
                                               Strict_Dummy_theta_I,
                                               Strict_Dummy_masks)

Smoothed_US_Strict_All <- Modified_Branzon_Frazier(Filter_US_Strict_All$W.updated,
                                                   Filter_US_Strict_All$W.predicted,
                                                   Filter_US_Strict_All$P.updated,
                                                   Filter_US_Strict_All$P.predicted,
                                                   Filter_US_Strict_All$Grads, 
                                                   Pop_vector,
                                                   Filter_US_Strict_All$Gain,
                                                   Filter_US_Strict_All$M.predicted,
                                                   Filter_US_Strict_All$resid)

var_here = Smoothed_US_Strict_All$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Stricter all: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Strict_All)

# Stricter stay-at-home
#----------------------
Filter_US_Strict_StayAtHome <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                      Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                      Weight_commuting, Weight_migration, Pop_vector,
                                                      Dummy_theta_F,
                                                      Strict_Dummy_theta_I,
                                                      Dummy_masks)

Smoothed_US_Strict_StayAtHome <- Modified_Branzon_Frazier(Filter_US_Strict_StayAtHome$W.updated,
                                                          Filter_US_Strict_StayAtHome$W.predicted,
                                                          Filter_US_Strict_StayAtHome$P.updated,
                                                          Filter_US_Strict_StayAtHome$P.predicted,
                                                          Filter_US_Strict_StayAtHome$Grads, 
                                                          Pop_vector,
                                                          Filter_US_Strict_StayAtHome$Gain,
                                                          Filter_US_Strict_StayAtHome$M.predicted,
                                                          Filter_US_Strict_StayAtHome$resid)

var_here = Smoothed_US_Strict_StayAtHome$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Stricter stay-at-home: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Strict_StayAtHome)

# Stricter travel bans
#----------------------
Filter_US_Strict_TravelBan <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                     Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                     Weight_commuting, Weight_migration, Pop_vector,
                                                     Strict_Dummy_theta_F,
                                                     Dummy_theta_I,
                                                     Dummy_masks)

Smoothed_US_Strict_TravelBan <- Modified_Branzon_Frazier(Filter_US_Strict_TravelBan$W.updated,
                                                         Filter_US_Strict_TravelBan$W.predicted,
                                                         Filter_US_Strict_TravelBan$P.updated,
                                                         Filter_US_Strict_TravelBan$P.predicted,
                                                         Filter_US_Strict_TravelBan$Grads, 
                                                         Pop_vector,
                                                         Filter_US_Strict_TravelBan$Gain,
                                                         Filter_US_Strict_TravelBan$M.predicted,
                                                         Filter_US_Strict_TravelBan$resid)

var_here = Smoothed_US_Strict_TravelBan$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Stricter travel ban: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Strict_TravelBan)


# Stricter masks
#----------------------
Filter_US_Strict_Masks <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                 Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                 Weight_commuting, Weight_migration, Pop_vector,
                                                 Dummy_theta_F,
                                                 Dummy_theta_I,
                                                 Strict_Dummy_masks)

Smoothed_US_Strict_Masks <- Modified_Branzon_Frazier(Filter_US_Strict_Masks$W.updated,
                                                     Filter_US_Strict_Masks$W.predicted,
                                                     Filter_US_Strict_Masks$P.updated,
                                                     Filter_US_Strict_Masks$P.predicted,
                                                     Filter_US_Strict_Masks$Grads, 
                                                     Pop_vector,
                                                     Filter_US_Strict_Masks$Gain,
                                                     Filter_US_Strict_Masks$M.predicted,
                                                     Filter_US_Strict_Masks$resid)


var_here = Smoothed_US_Strict_Masks$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Stricter masks: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Strict_Masks)

# No regulation
#--------------
Filter_US_Loose_All <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                              Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                              Weight_commuting, Weight_migration, Pop_vector,
                                              Loose_Dummy_theta_F,
                                              Loose_Dummy_theta_I,
                                              Loose_Dummy_masks)

Smoothed_US_Loose_All <- Modified_Branzon_Frazier(Filter_US_Loose_All$W.updated,
                                                  Filter_US_Loose_All$W.predicted,
                                                  Filter_US_Loose_All$P.updated,
                                                  Filter_US_Loose_All$P.predicted,
                                                  Filter_US_Loose_All$Grads, 
                                                  Pop_vector,
                                                  Filter_US_Loose_All$Gain,
                                                  Filter_US_Loose_All$M.predicted,
                                                  Filter_US_Loose_All$resid)

var_here = Smoothed_US_Loose_All$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Loose all: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Loose_All)

# Loose stay-at-home
#----------------------
Filter_US_Loose_StayAtHome <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                     Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                     Weight_commuting, Weight_migration, Pop_vector,
                                                     Dummy_theta_F,
                                                     Loose_Dummy_theta_I,
                                                     Dummy_masks)

Smoothed_US_Loose_StayAtHome <- Modified_Branzon_Frazier(Filter_US_Loose_StayAtHome$W.updated,
                                                         Filter_US_Loose_StayAtHome$W.predicted,
                                                         Filter_US_Loose_StayAtHome$P.updated,
                                                         Filter_US_Loose_StayAtHome$P.predicted,
                                                         Filter_US_Loose_StayAtHome$Grads, 
                                                         Pop_vector,
                                                         Filter_US_Loose_StayAtHome$Gain,
                                                         Filter_US_Loose_StayAtHome$M.predicted,
                                                         Filter_US_Loose_StayAtHome$resid)

var_here = Smoothed_US_Loose_StayAtHome$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Loose stay-at-home: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Loose_StayAtHome)

# Looser travel bans
#----------------------
Filter_US_Loose_TravelBan <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                    Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                    Weight_commuting, Weight_migration, Pop_vector,
                                                    Loose_Dummy_theta_F,
                                                    Dummy_theta_I,
                                                    Dummy_masks)

Smoothed_US_Loose_TravelBan <- Modified_Branzon_Frazier(Filter_US_Loose_TravelBan$W.updated,
                                                        Filter_US_Loose_TravelBan$W.predicted,
                                                        Filter_US_Loose_TravelBan$P.updated,
                                                        Filter_US_Loose_TravelBan$P.predicted,
                                                        Filter_US_Loose_TravelBan$Grads, 
                                                        Pop_vector,
                                                        Filter_US_Loose_TravelBan$Gain,
                                                        Filter_US_Loose_TravelBan$M.predicted,
                                                        Filter_US_Loose_TravelBan$resid)

var_here = Smoothed_US_Loose_TravelBan$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Loose travel ban: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Loose_TravelBan)

# Looser masks
#----------------------
Filter_US_Loose_Masks <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                Weight_commuting, Weight_migration, Pop_vector,
                                                Dummy_theta_F,
                                                Dummy_theta_I,
                                                Loose_Dummy_masks)

Smoothed_US_Loose_Masks <- Modified_Branzon_Frazier(Filter_US_Loose_Masks$W.updated,
                                                    Filter_US_Loose_Masks$W.predicted,
                                                    Filter_US_Loose_Masks$P.updated,
                                                    Filter_US_Loose_Masks$P.predicted,
                                                    Filter_US_Loose_Masks$Grads, 
                                                    Pop_vector,
                                                    Filter_US_Loose_Masks$Gain,
                                                    Filter_US_Loose_Masks$M.predicted,
                                                    Filter_US_Loose_Masks$resid)

var_here = Smoothed_US_Loose_Masks$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Loose masks: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Loose_Masks)

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#
# Other counterfactuals
#
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
# EARLY MASKS
Filter_US_Super_Strict_Masks <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                       Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                       Weight_commuting, Weight_migration, Pop_vector,
                                                       Dummy_theta_F,
                                                       Dummy_theta_I,
                                                       Strict_Dummy_theta_I)
Smoothed_US_Super_Strict_Masks <- Modified_Branzon_Frazier(Filter_US_Super_Strict_Masks$W.updated,
                                                           Filter_US_Super_Strict_Masks$W.predicted,
                                                           Filter_US_Super_Strict_Masks$P.updated,
                                                           Filter_US_Super_Strict_Masks$P.predicted,
                                                           Filter_US_Super_Strict_Masks$Grads, 
                                                           Pop_vector,
                                                           Filter_US_Super_Strict_Masks$Gain,
                                                           Filter_US_Super_Strict_Masks$M.predicted,
                                                           Filter_US_Super_Strict_Masks$resid)

var_here = Smoothed_US_Super_Strict_Masks$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Masks as early as stay-at-home: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Super_Strict_Masks)

# EARLY MASKS 2
Filter_US_Super_Strict_Masks2 <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                        Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                        Weight_commuting, Weight_migration, Pop_vector,
                                                        Loose_Dummy_theta_F,
                                                        Loose_Dummy_theta_I,
                                                        Strict_Dummy_theta_I)
Smoothed_US_Super_Strict_Masks2 <- Modified_Branzon_Frazier(Filter_US_Super_Strict_Masks2$W.updated,
                                                            Filter_US_Super_Strict_Masks2$W.predicted,
                                                            Filter_US_Super_Strict_Masks2$P.updated,
                                                            Filter_US_Super_Strict_Masks2$P.predicted,
                                                            Filter_US_Super_Strict_Masks2$Grads, 
                                                            Pop_vector,
                                                            Filter_US_Super_Strict_Masks2$Gain,
                                                            Filter_US_Super_Strict_Masks2$M.predicted,
                                                            Filter_US_Super_Strict_Masks2$resid)

var_here = Smoothed_US_Super_Strict_Masks2$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Masks as early as stay-at-home without stay-at-home or travel bans: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Super_Strict_Masks2)

# EARLY TRAVELBAN
Filter_US_Travel_Feb12 <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                 Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                 Weight_commuting, Weight_migration, Pop_vector,
                                                 Early_dummy_travel_Feb12,
                                                 Dummy_theta_I,
                                                 Dummy_masks)
Smoothed_US_Travel_Feb12 <- Modified_Branzon_Frazier(Filter_US_Travel_Feb12$W.updated,
                                                     Filter_US_Travel_Feb12$W.predicted,
                                                     Filter_US_Travel_Feb12$P.updated,
                                                     Filter_US_Travel_Feb12$P.predicted,
                                                     Filter_US_Travel_Feb12$Grads, 
                                                     Pop_vector,
                                                     Filter_US_Travel_Feb12$Gain,
                                                     Filter_US_Travel_Feb12$M.predicted,
                                                     Filter_US_Travel_Feb12$resid)

var_here = Smoothed_US_Travel_Feb12$W.smooth
res_here = sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
print("Travel ban by February 12: Completed")
print(round(res_here[length(res_here)]))

rm(Filter_US_Travel_Feb12)

#==============================================================
#### R naught reduction PLOT
#==============================================================
Early_dummy_masks <- matrix(Dates.selec<"2020-03-17" , ncol = 1) %*% matrix(1, 1, nb_states)

premult_theta_I_baseline  <- NULL
premult_theta_I_cf        <- NULL

for(i in order.to.plot.states){
  
  index.to.plot <- 4 * nb_states + i
  premult_theta_I <- (Parameters$Mask_reduc + Dummy_masks[,Selec.state[i]] * (1-Parameters$Mask_reduc)) * (
    Parameters$theta_I_low[Selec.state[i]] + 
      Dummy_theta_I[,Selec.state[i]] *(
        Parameters$theta_I_high[Selec.state[i]] - Parameters$theta_I_low[Selec.state[i]]
      )
  )
  premult_theta_I_baseline <- cbind(premult_theta_I_baseline, premult_theta_I)
  
  premult_theta_I <- (Parameters$Mask_reduc + Early_dummy_masks[,Selec.state[i]] * (1-Parameters$Mask_reduc)) * (
    Parameters$theta_I_low[Selec.state[i]] + 
      Dummy_theta_I[,Selec.state[i]] *(
        Parameters$theta_I_high[Selec.state[i]] - Parameters$theta_I_low[Selec.state[i]]
      )
  )
  premult_theta_I_cf <- cbind(premult_theta_I_cf, premult_theta_I)
  
}

par(mai = c(.5, 0.5,.3,.1))
MY = max(apply(diag(Pop_vector) %*% (t(premult_theta_I_baseline) * Smoothed_Betas)/sum(Pop_vector),2,sum)/
           (Parameters$delta + Parameters$gamma), apply(diag(Pop_vector) %*% (t(premult_theta_I_cf) * Smoothed_Betas)/sum(Pop_vector),2,sum)/
           (Parameters$delta + Parameters$gamma))

plot(Dates.selec, apply(diag(Pop_vector) %*% (t(premult_theta_I_baseline) * Smoothed_Betas)/sum(Pop_vector),2,sum)/
       (Parameters$delta + Parameters$gamma), 
     type='l', 
     main = "Effective reproduction number in the U.S.", cex.main = 1,
     panel.first = grid(),
     col = "black",
     xlab = "", ylab = "total", 
     ylim=c(0, MY), lwd = 1)

lines(Dates.selec, apply(diag(Pop_vector) %*% (t(premult_theta_I_cf) * Smoothed_Betas)/sum(Pop_vector),2,sum)/
        (Parameters$delta + Parameters$gamma), 
      col = 'red', lwd = 2)

abline(h=1, lwd = 2, col = 'blue')
abline(v = as.Date("2020-03-17"), lty = "dashed")
text(x = Dates.selec[which(Dates.selec=="2020-03-17")], y = 2.5, labels ="March 17", srt = 270*3, pos = 2, cex = .8)

legend("topright",
       c("Data", "Counterfactual", expression(paste(R[0], " = 1", sep = ""))),
       col = c("black", "red", "blue"), 
       lty = c("solid", "solid", "solid"), lwd = c(1,2,2,2,2), bg = "white", cex = 0.9)


#==================================
#### Infected flows PLOT
#==================================
library(zoo)
states_here   <- vector.of.states
q_dates       <- which(Dates.selec >= as.Date("2020-03-17"))

for(i in order.to.plot.states){
  
  state_here          <- states_here[i]
  
  index.to.plot.I     <- 2 * nb_states + i
  index.to.plot.R     <- 3 * nb_states + i
  index.to.plot.D     <- 0 * nb_states + i
  
  index.to.plot.I.cf  <- (5 + 2) * nb_states + i
  index.to.plot.R.cf  <- (5 + 3) * nb_states + i
  index.to.plot.D.cf  <- (5 + 0) * nb_states + i
  

  inf_data  <- cummax(smooth(Smoothed_US_Strict_TravelBan$W.smooth[index.to.plot.I,q_dates] + Smoothed_US_Strict_TravelBan$W.smooth[index.to.plot.R,q_dates] + Smoothed_US_Strict_TravelBan$W.smooth[index.to.plot.D,q_dates]))
  inf_cf    <- cummax(smooth(Smoothed_US_Strict_TravelBan$W.smooth[index.to.plot.I.cf,q_dates] + Smoothed_US_Strict_TravelBan$W.smooth[index.to.plot.R.cf,q_dates] + Smoothed_US_Strict_TravelBan$W.smooth[index.to.plot.D.cf,q_dates]))
  
  plot(Dates.selec[q_dates], inf_data, type='l', 
       cex.main = 1, 
       main = paste("Estimated number of infections (", vector.of.states[Selec.state[i]], ")", sep = ""), 
       panel.first = grid(),
       xlab = "", ylab = "", 
       col = "black", lwd = 1,
       ylim=c(min(inf_data, inf_cf,na.rm = TRUE),
              max(inf_data, inf_cf,na.rm = TRUE)))
  
  lines(Dates.selec[q_dates], inf_cf, col = "red", lwd = 2)
  
  legend(x = "topleft", legend = c("Data", "Counterfactual"), lwd = c(1, 2), col = c("black", "red"))
}







# STATE BREAKDOWN PLOTS #

# Strict overall
#---------------
var_here                          <- Smoothed_US_Strict_All$W.smooth
res_strict_all_cf_perstate        <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_strict_all_baseline_perstate  <- rowSums(var_here[1:nb_states,])
res_strict_all_country            <- round(sum(res_strict_all_cf_perstate) - sum(res_strict_all_baseline_perstate))
res_strict_all_perstate           <- round(res_strict_all_cf_perstate - res_strict_all_baseline_perstate)
names(res_strict_all_perstate)    <- names(res_strict_all_baseline_perstate) <- states_here
res_strict_all_perstate_rel       <- res_strict_all_perstate[states_here]/res_strict_all_baseline_perstate[states_here]

write.csv(res_strict_all_perstate, "cf_joint_abs_All_strict.csv")
write.csv(res_strict_all_perstate_rel, "cf_joint_rel_All_strict.csv")

# Stricter stay-at-home
#----------------------
var_here                                <- Smoothed_US_Strict_StayAtHome$W.smooth
res_strict_stayathome_cf_perstate       <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_strict_stayathome_baseline_perstate <- rowSums(var_here[1:nb_states,])
res_strict_stayathome_country           <- round(sum(res_strict_stayathome_cf_perstate) - sum(res_strict_stayathome_baseline_perstate))
res_strict_stayathome_perstate          <- round(res_strict_stayathome_cf_perstate - res_strict_stayathome_baseline_perstate)
names(res_strict_stayathome_perstate)   <- names(res_strict_stayathome_baseline_perstate) <- states_here
res_strict_stayathome_perstate_rel      <- res_strict_stayathome_perstate[states_here]/res_strict_stayathome_baseline_perstate[states_here]

write.csv(res_strict_stayathome_perstate, "cf_joint_abs_StayAtHome_strict.csv")
write.csv(res_strict_stayathome_perstate_rel, "cf_joint_rel_StayAtHome_strict.csv")

# Stricter travel bans
#----------------------
var_here                                <- Smoothed_US_Strict_TravelBan$W.smooth
res_strict_travel_cf_perstate           <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_strict_travel_baseline_perstate     <- rowSums(var_here[1:nb_states,])
res_strict_travel_country               <- round(sum(res_strict_travel_cf_perstate) - sum(res_strict_travel_baseline_perstate))
res_strict_travel_perstate              <- round(res_strict_travel_cf_perstate - res_strict_travel_baseline_perstate)
names(res_strict_travel_perstate)       <- names(res_strict_travel_baseline_perstate) <- states_here
res_strict_travel_perstate_rel          <- res_strict_travel_perstate[states_here]/res_strict_travel_baseline_perstate[states_here]

write.csv(res_strict_travel_perstate, "cf_joint_abs_TravelBan_strict.csv")
write.csv(res_strict_travel_perstate_rel, "cf_joint_rel_TravelBan_strict.csv")

# Stricter masks
#----------------------
var_here                                <- Smoothed_US_Strict_Masks$W.smooth
res_strict_masks_cf_perstate            <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_strict_masks_baseline_perstate      <- rowSums(var_here[1:nb_states,])
res_strict_masks_country                <- round(sum(res_strict_masks_cf_perstate) - sum(res_strict_masks_baseline_perstate))
res_strict_masks_perstate               <- round(res_strict_masks_cf_perstate - res_strict_masks_baseline_perstate)
res_strict_masks_perstate[res_strict_masks_perstate > 0] <- -0.01
names(res_strict_masks_perstate)        <- names(res_strict_masks_baseline_perstate) <- states_here
res_strict_masks_perstate_rel           <- res_strict_masks_perstate[states_here]/res_strict_masks_baseline_perstate[states_here]

write.csv(res_strict_masks_perstate, "cf_joint_abs_Masks_strict.csv")
write.csv(res_strict_masks_perstate_rel, "cf_joint_rel_Masks_strict.csv")

# Loose overall
#---------------

var_here                          <- Smoothed_US_Loose_All$W.smooth
res_loose_all_cf_perstate         <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_loose_all_baseline_perstate   <- rowSums(var_here[1:nb_states,])
res_loose_all_country             <- round(sum(res_loose_all_cf_perstate) - sum(res_loose_all_baseline_perstate))
res_loose_all_perstate            <- round(res_loose_all_cf_perstate - res_loose_all_baseline_perstate)
names(res_loose_all_perstate)     <- names(res_loose_all_baseline_perstate) <- states_here
res_loose_all_perstate_rel        <- res_loose_all_perstate[states_here]/res_loose_all_baseline_perstate[states_here]

write.csv(res_loose_all_perstate, "cf_joint_abs_All_loose.csv")
write.csv(res_loose_all_perstate_rel, "cf_joint_rel_All_loose.csv")

# Looseer stay-at-home
#----------------------
var_here                                <- Smoothed_US_Loose_StayAtHome$W.smooth
res_loose_stayathome_cf_perstate        <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_loose_stayathome_baseline_perstate  <- rowSums(var_here[1:nb_states,])
res_loose_stayathome_country            <- round(sum(res_loose_stayathome_cf_perstate) - sum(res_loose_stayathome_baseline_perstate))
res_loose_stayathome_perstate           <- round(res_loose_stayathome_cf_perstate - res_loose_stayathome_baseline_perstate)
names(res_loose_stayathome_perstate)    <- names(res_loose_stayathome_baseline_perstate) <- states_here
res_loose_stayathome_perstate_rel       <- res_loose_stayathome_perstate[states_here]/res_loose_stayathome_baseline_perstate[states_here]

write.csv(res_loose_stayathome_perstate, "cf_joint_abs_StayAtHome_loose.csv")
write.csv(res_loose_stayathome_perstate_rel, "cf_joint_rel_StayAtHome_loose.csv")

# Looseer travel bans
#----------------------
var_here                            <- Smoothed_US_Loose_TravelBan$W.smooth
res_loose_travel_cf_perstate        <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_loose_travel_baseline_perstate  <- rowSums(var_here[1:nb_states,])
res_loose_travel_country            <- round(sum(res_loose_travel_cf_perstate) - sum(res_loose_travel_baseline_perstate))
res_loose_travel_perstate           <- round(res_loose_travel_cf_perstate - res_loose_travel_baseline_perstate)
names(res_loose_travel_perstate)    <- names(res_loose_travel_baseline_perstate) <- states_here
res_loose_travel_perstate_rel       <- res_loose_travel_perstate[states_here]/res_loose_travel_baseline_perstate[states_here]

write.csv(res_loose_travel_perstate, "cf_joint_abs_Travelan_loose.csv")
write.csv(res_loose_travel_perstate_rel, "cf_joint_rel_TravelBan_loose.csv")

# Looseer masks
#----------------------
var_here                            <- Smoothed_US_Loose_Masks$W.smooth
res_loose_masks_cf_perstate         <- rowSums(var_here[(5 * nb_states + 1):(6 * nb_states),])
res_loose_masks_baseline_perstate   <- rowSums(var_here[1:nb_states,])
res_loose_masks_country             <- round(sum(res_loose_masks_cf_perstate) - sum(res_loose_masks_baseline_perstate))
res_loose_masks_perstate            <- round(res_loose_masks_cf_perstate - res_loose_masks_baseline_perstate)
names(res_loose_masks_perstate)     <- names(res_loose_masks_baseline_perstate) <- states_here
res_loose_masks_perstate_rel        <- res_loose_masks_perstate[states_here]/res_loose_masks_baseline_perstate[states_here]

write.csv(res_loose_masks_perstate, "cf_joint_abs_Masks_loose.csv")
write.csv(res_loose_masks_perstate_rel, "cf_joint_rel_Masks_loose.csv")


#=======================================================
#=======================================================
### PLOTS FOR FEDERAL MANDATES #### 
#=======================================================
#=======================================================
statenames_here = c("Alaska", "Alabama", "Arkansas", "Arizona", 
                    "California", "Colorado", "Connecticut", "D. C.", 
                    "Delaware", "Florida", "Georgia", "Hawaii", "Iowa", 
                    "Idaho", "Illinois", "Indiana", "Kansas", "Kentucky", 
                    "Louisiana", "Massachusetts", "Maryland", "Maine", "Michigan", 
                    "Minnesota", "Missouri", "Mississippi", "Montana", 
                    "North Carolina", "North Dakota", "Nebraska", "New Hampshire",
                    "New Jersey", "New Mexico", "Nevada", "New York", "Ohio", 
                    "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", 
                    "South Carolina", "South Dakota", "Tennessee", "Texas", 
                    "Utah", "Virginia", "Vermont", "Washington", "Wisconsin", 
                    "West Virginia", "Wyoming")


# all
#-----------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_all_perstate, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states adopt the three policies simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.20 * (Plot1.US.Aggreg$x < 0)),
                    min(abs(range(Plot1.US.Aggreg$x)))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_all_perstate, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -min(range(Plot2.US.Aggreg$x)),
                    max(max(abs(Plot2.US.Aggreg$x))*.28 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)



#-----------------
# Stay-at-home
#-----------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))


# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_stayathome_perstate, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states adopt a stay-at-home order simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.20 * (Plot1.US.Aggreg$x < 0)),
                    min(abs(range(Plot1.US.Aggreg$x)))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_stayathome_perstate, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -min(range(Plot2.US.Aggreg$x)),
                    max(max(abs(Plot2.US.Aggreg$x))*.28 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)

#-----------------
# Masks
#-----------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))


# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_masks_perstate, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states adopt a mask mandate simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.20 * (Plot1.US.Aggreg$x < 0)),
                    min(abs(range(Plot1.US.Aggreg$x)))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_masks_perstate, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -min(range(Plot2.US.Aggreg$x)),
                    max(max(abs(Plot2.US.Aggreg$x))*.28 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)


#-----------------
# Travel ban
#-----------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))


# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_travel_perstate, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states ban border crossings simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.55 * (Plot1.US.Aggreg$x < 0)),
                    max(max(abs(Plot1.US.Aggreg$x))*.35 * (Plot1.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_travel_perstate, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -max(max(abs(Plot2.US.Aggreg$x))*.35 * (Plot2.US.Aggreg$x < 0)),
                    max(max(abs(Plot2.US.Aggreg$x))*.35 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)


#-----------------
# all - relative
#-----------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))


# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_all_perstate_rel, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states adopt the three policies simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.3 * (Plot1.US.Aggreg$x < 0)),
                    min(abs(range(Plot1.US.Aggreg$x)))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_all_perstate_rel, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -min(range(Plot2.US.Aggreg$x)),
                    max(max(abs(Plot2.US.Aggreg$x))*.22 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)


#-------------------------
# Stay at home - relative
#-------------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))


# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_stayathome_perstate_rel, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states adopt a stay-at-home order simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.3 * (Plot1.US.Aggreg$x < 0)),
                    min(abs(range(Plot1.US.Aggreg$x)))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_stayathome_perstate_rel, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -min(range(Plot2.US.Aggreg$x)),
                    max(max(abs(Plot2.US.Aggreg$x))*.22 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)

#-------------------------
# Masks - relative
#-------------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))


# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_masks_perstate_rel, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states adopt a mask mandate simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.3 * (Plot1.US.Aggreg$x < 0)),
                    min(abs(range(Plot1.US.Aggreg$x)))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_masks_perstate_rel, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -min(range(Plot2.US.Aggreg$x)),
                    max(max(abs(Plot2.US.Aggreg$x))*.22 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)

#-------------------------
# Travel ban - relative
#-------------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))


# STRICT SCENARIO
Plot1.US.Aggreg <-sort(res_strict_travel_perstate_rel, index.return = T)
midpts <- barplot(Plot1.US.Aggreg$x, 
                  names.arg = statenames_here[Plot1.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("All states ban border crossings simultaneously")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.US.Aggreg$x) + c(
                    -max(max(abs(Plot1.US.Aggreg$x))*.35 * (Plot1.US.Aggreg$x < 0)),
                    max(max(abs(Plot1.US.Aggreg$x))*.375 * (Plot1.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot1.US.Aggreg$x < 0
index.pos <- Plot1.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg], x = Plot1.US.Aggreg$x[index.neg] - 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos], x = Plot1.US.Aggreg$x[index.pos] + 
         max(abs(Plot1.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot1.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot1.US.Aggreg$x, 
        names.arg = statenames_here[Plot1.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.US.Aggreg <-sort(res_loose_travel_perstate_rel, index.return = T)
midpts <- barplot(Plot2.US.Aggreg$x, 
                  names.arg = statenames_here[Plot2.US.Aggreg$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.US.Aggreg$x) + c(
                    -max(max(abs(Plot2.US.Aggreg$x))*.45 * (Plot2.US.Aggreg$x < 0)),
                    max(max(abs(Plot2.US.Aggreg$x))*.35 * (Plot2.US.Aggreg$x > 0))
                  ))
grid()
index.neg <- Plot2.US.Aggreg$x < 0
index.pos <- Plot2.US.Aggreg$x > 0
if (length(which(index.neg) == T) > 0) {
  text(y = midpts[index.neg]-0.1, x = Plot2.US.Aggreg$x[index.neg] - 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.neg], pos = 2,  cex = .6)
}
if (length(which(index.pos) == T) > 0) {
  text(y = midpts[index.pos]-0.1, x = Plot2.US.Aggreg$x[index.pos] + 
         max(abs(Plot2.US.Aggreg$x))*0 , 
       labels = statenames_here[Plot2.US.Aggreg$ix][index.pos], pos = 4, cex = .6)
}
barplot(Plot2.US.Aggreg$x, 
        names.arg = statenames_here[Plot2.US.Aggreg$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)



#=============================================
# Early mask mandate
#=============================================
seq_dates_here  <- seq(as.Date("2020-03-17"), to = as.Date("2020-04-17"), by = "1 day")
N_dates_here    <- length(seq_dates_here)
mask_here       <- data.frame("date" = seq_dates_here, "prev_death" = NA)

for (ii in N_dates_here:1) {
  day_here          <- seq_dates_here[ii]
  Early_dummy_here  <- matrix(Dates.selec<day_here , ncol = 1) %*% matrix(1, 1, nb_states)
  
  Filter_US_Super_Strict_Masks    <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                         Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                         Weight_commuting, Weight_migration, Pop_vector,
                                                         Dummy_theta_F,
                                                         Dummy_theta_I,
                                                         Early_dummy_here)
  Smoothed_US_Super_Strict_Masks  <- Modified_Branzon_Frazier(Filter_US_Super_Strict_Masks$W.updated,
                                                             Filter_US_Super_Strict_Masks$W.predicted,
                                                             Filter_US_Super_Strict_Masks$P.updated,
                                                             Filter_US_Super_Strict_Masks$P.predicted,
                                                             Filter_US_Super_Strict_Masks$Grads, 
                                                             Pop_vector,
                                                             Filter_US_Super_Strict_Masks$Gain,
                                                             Filter_US_Super_Strict_Masks$M.predicted,
                                                             Filter_US_Super_Strict_Masks$resid)
  
  var_here <- Smoothed_US_Super_Strict_Masks$W.smooth
  res_here <- sum(var_here[(5 * nb_states + 1):(6 * nb_states),]) - sum(var_here[1:nb_states,])
  
  mask_here[ii,"prev_death"] <- res_here
  
  print(c(as.character(day_here), res_here))
}
write.csv(mask_here, "early_mask_mandate.csv")

# Plot 
plot(as.Date(mask_here[,"date"]), abs(mask_here[,"prev_death"]), 
     lwd = 2, type='l', 
     main = "Number of preventable deaths with an early federal mask mandate", 
     cex.main = 1, panel.first = grid(), col = "red", xlab = "Date mandate goes into effect", ylab = "Preventable deaths")
abline(v = as.Date("2020-03-17"), lty = "dashed", lwd = 1, col = "black")  
abline(v = as.Date("2020-03-20"), lty = "dashed", lwd = 1, col = "black")  
abline(v = as.Date("2020-04-17"), lty = "dashed", lwd = 1, col = "black")  

text(x = as.Date("2020-03-17"), y = 140000, labels ="March 17", srt = 90*3, pos = 4, cex  =.8)
text(x = as.Date("2020-03-20"), y = 140000, labels ="March 20", srt = 90*3, pos = 4, cex  =.8)
text(x = as.Date("2020-04-17"), y = 140000, labels ="April 17", srt = 90*3, pos = 4, cex  =.8)




##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#
#     RUN COUNTERFACTUALS WITH STATE-BY-STATE POLICIES 
#
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
library(doParallel)
library(parallel)

# Create cluster to parallelize
stopImplicitCluster()
nw <- detectCores()-1
registerDoParallel(cores = nw)

Smoothed_PerState_Strict_All          <- list()
Smoothed_PerState_Strict_StayAtHome   <- list()
Smoothed_PerState_Strict_TravelBan    <- list()
Smoothed_PerState_Strict_Masks        <- list()

Smoothed_PerState_Loose_All           <- list()
Smoothed_PerState_Loose_StayAtHome    <- list()
Smoothed_PerState_Loose_TravelBan     <- list()
Smoothed_PerState_Loose_Masks         <- list()

stricter_F_vector     <- apply(Dummy_theta_F,1,min)
stricter_I_vector     <- apply(Dummy_theta_I,1,min)
stricter_masks_vector <- apply(Dummy_masks,1,min)

print("Start state-by-state counterfactuals")

tempforeach <- foreach(i=1:nb_states, .errorhandling = 'stop') %do% {
  
  # Define dummies
  #---------------
  Strict_Dummy_theta_I.PerState      <- Dummy_theta_I
  Strict_Dummy_theta_I.PerState[,i]  <- stricter_I
  Strict_Dummy_theta_F.PerState      <- Dummy_theta_F
  Strict_Dummy_theta_F.PerState[,i]  <- stricter_F
  Strict_Dummy_masks.PerState        <- Dummy_masks
  Strict_Dummy_masks.PerState[,i]    <- stricter_masks
  
  Loose_Dummy_theta_F.PerState       <- Dummy_theta_F
  Loose_Dummy_theta_F.PerState[,i]   <- 1
  Loose_Dummy_theta_I.PerState       <- Dummy_theta_I
  Loose_Dummy_theta_I.PerState[,i]   <- 1
  Loose_Dummy_masks.PerState         <- Dummy_masks
  Loose_Dummy_masks.PerState[,i]     <- 1
  
  # Strict overall
  #---------------
  Filter_State_Strict_All   <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                      Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                      Weight_commuting, Weight_migration, Pop_vector,
                                                      Strict_Dummy_theta_F.PerState,
                                                      Strict_Dummy_theta_I.PerState,
                                                      Strict_Dummy_masks.PerState)
  Smoothed_State_Strict_All <- Modified_Branzon_Frazier(Filter_State_Strict_All$W.updated,
                                                        Filter_State_Strict_All$W.predicted,
                                                        Filter_State_Strict_All$P.updated,
                                                        Filter_State_Strict_All$P.predicted,
                                                        Filter_State_Strict_All$Grads, 
                                                        Pop_vector,
                                                        Filter_State_Strict_All$Gain,
                                                        Filter_State_Strict_All$M.predicted,
                                                        Filter_State_Strict_All$resid)
  rm(Filter_State_Strict_All)
  
  # Stricter stay-at-home
  #----------------------
  Filter_State_Strict_StayAtHome    <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                              Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                              Weight_commuting, Weight_migration, Pop_vector,
                                                              Dummy_theta_F,
                                                              Strict_Dummy_theta_I.PerState,
                                                              Dummy_masks)
  Smoothed_State_Strict_StayAtHome  <- Modified_Branzon_Frazier(Filter_State_Strict_StayAtHome$W.updated,
                                                                Filter_State_Strict_StayAtHome$W.predicted,
                                                                Filter_State_Strict_StayAtHome$P.updated,
                                                                Filter_State_Strict_StayAtHome$P.predicted,
                                                                Filter_State_Strict_StayAtHome$Grads, 
                                                                Pop_vector,
                                                                Filter_State_Strict_StayAtHome$Gain,
                                                                Filter_State_Strict_StayAtHome$M.predicted,
                                                                Filter_State_Strict_StayAtHome$resid)
  rm(Filter_State_Strict_StayAtHome)
  
  # Stricter travel bans
  #----------------------
  Filter_State_Strict_TravelBan   <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                            Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                            Weight_commuting, Weight_migration, Pop_vector,
                                                            Strict_Dummy_theta_F.PerState,
                                                            Dummy_theta_I,
                                                            Dummy_masks)
  Smoothed_State_Strict_TravelBan <- Modified_Branzon_Frazier(Filter_State_Strict_TravelBan$W.updated,
                                                              Filter_State_Strict_TravelBan$W.predicted,
                                                              Filter_State_Strict_TravelBan$P.updated,
                                                              Filter_State_Strict_TravelBan$P.predicted,
                                                              Filter_State_Strict_TravelBan$Grads, 
                                                              Pop_vector,
                                                              Filter_State_Strict_TravelBan$Gain,
                                                              Filter_State_Strict_TravelBan$M.predicted,
                                                              Filter_State_Strict_TravelBan$resid)
  rm(Filter_State_Strict_TravelBan)
  
  # Stricter masks
  #----------------------
  Filter_State_Strict_Masks   <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                        Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                        Weight_commuting, Weight_migration, Pop_vector,
                                                        Dummy_theta_F,
                                                        Dummy_theta_I,
                                                        Strict_Dummy_masks.PerState)
  Smoothed_State_Strict_Masks <- Modified_Branzon_Frazier(Filter_State_Strict_Masks$W.updated,
                                                          Filter_State_Strict_Masks$W.predicted,
                                                          Filter_State_Strict_Masks$P.updated,
                                                          Filter_State_Strict_Masks$P.predicted,
                                                          Filter_State_Strict_Masks$Grads, 
                                                          Pop_vector,
                                                          Filter_State_Strict_Masks$Gain,
                                                          Filter_State_Strict_Masks$M.predicted,
                                                          Filter_State_Strict_Masks$resid)
  rm(Filter_State_Strict_Masks)
  
  # No regulation
  #--------------
  Filter_State_Loose_All    <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                      Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                      Weight_commuting, Weight_migration, Pop_vector,
                                                      Loose_Dummy_theta_F.PerState,
                                                      Loose_Dummy_theta_I.PerState,
                                                      Loose_Dummy_masks.PerState)
  Smoothed_State_Loose_All  <- Modified_Branzon_Frazier(Filter_State_Loose_All$W.updated,
                                                        Filter_State_Loose_All$W.predicted,
                                                        Filter_State_Loose_All$P.updated,
                                                        Filter_State_Loose_All$P.predicted,
                                                        Filter_State_Loose_All$Grads, 
                                                        Pop_vector,
                                                        Filter_State_Loose_All$Gain,
                                                        Filter_State_Loose_All$M.predicted,
                                                        Filter_State_Loose_All$resid)
  rm(Filter_State_Loose_All)
  
  # Loose stay-at-home
  #----------------------
  Filter_State_Loose_StayAtHome   <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                            Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                            Weight_commuting, Weight_migration, Pop_vector,
                                                            Dummy_theta_F,
                                                            Loose_Dummy_theta_I.PerState,
                                                            Dummy_masks)
  Smoothed_State_Loose_StayAtHome <- Modified_Branzon_Frazier(Filter_State_Loose_StayAtHome$W.updated,
                                                              Filter_State_Loose_StayAtHome$W.predicted,
                                                              Filter_State_Loose_StayAtHome$P.updated,
                                                              Filter_State_Loose_StayAtHome$P.predicted,
                                                              Filter_State_Loose_StayAtHome$Grads, 
                                                              Pop_vector,
                                                              Filter_State_Loose_StayAtHome$Gain,
                                                              Filter_State_Loose_StayAtHome$M.predicted,
                                                              Filter_State_Loose_StayAtHome$resid)
  rm(Filter_State_Loose_StayAtHome)
  
  # Looser travel bans
  #----------------------
  Filter_State_Loose_TravelBan    <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                            Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                            Weight_commuting, Weight_migration, Pop_vector,
                                                            Loose_Dummy_theta_F.PerState,
                                                            Dummy_theta_I,
                                                            Dummy_masks)
  Smoothed_State_Loose_TravelBan  <- Modified_Branzon_Frazier(Filter_State_Loose_TravelBan$W.updated,
                                                              Filter_State_Loose_TravelBan$W.predicted,
                                                              Filter_State_Loose_TravelBan$P.updated,
                                                              Filter_State_Loose_TravelBan$P.predicted,
                                                              Filter_State_Loose_TravelBan$Grads, 
                                                              Pop_vector,
                                                              Filter_State_Loose_TravelBan$Gain,
                                                              Filter_State_Loose_TravelBan$M.predicted,
                                                              Filter_State_Loose_TravelBan$resid)
  rm(Filter_State_Loose_TravelBan)
  
  # Looser masks
  #----------------------
  Filter_State_Loose_Masks    <- EKF_filter_cpp_counter(Parameters, W0, P0, Death.selec,
                                                        Dummy_theta_F, Dummy_theta_I, Dummy_masks,
                                                        Weight_commuting, Weight_migration, Pop_vector,
                                                        Dummy_theta_F,
                                                        Dummy_theta_I,
                                                        Loose_Dummy_masks.PerState)
  Smoothed_State_Loose_Masks  <- Modified_Branzon_Frazier(Filter_State_Loose_Masks$W.updated,
                                                          Filter_State_Loose_Masks$W.predicted,
                                                          Filter_State_Loose_Masks$P.updated,
                                                          Filter_State_Loose_Masks$P.predicted,
                                                          Filter_State_Loose_Masks$Grads, 
                                                          Pop_vector,
                                                          Filter_State_Loose_Masks$Gain,
                                                          Filter_State_Loose_Masks$M.predicted,
                                                          Filter_State_Loose_Masks$resid)
  rm(Filter_State_Loose_Masks)
  
  
  # PRINT PROGRESS
  #---------------
  print(str_c(vector.of.states[i], " ", 100*round(i/nb_states,4),"%"))
  
  # FINAL OBJECT TO SEND BACK
  #--------------------------
  list("strict_all"         = Smoothed_State_Strict_All$W.smooth, 
       "strict_stayathome"  = Smoothed_State_Strict_StayAtHome$W.smooth, 
       "strict_travelban"   = Smoothed_State_Strict_TravelBan$W.smooth, 
       "strict_masks"       = Smoothed_State_Strict_Masks$W.smooth, 
       "loose_all"          = Smoothed_State_Loose_All$W.smooth, 
       "loose_stayathome"   = Smoothed_State_Loose_StayAtHome$W.smooth, 
       "loose_travelban"    = Smoothed_State_Loose_TravelBan$W.smooth, 
       "loose_masks"        = Smoothed_State_Loose_Masks$W.smooth)
  
}

Smoothed_PerState_Strict_All          <- lapply(tempforeach, "[[", 1)
Smoothed_PerState_Strict_StayAtHome   <- lapply(tempforeach, "[[", 2)
Smoothed_PerState_Strict_TravelBan    <- lapply(tempforeach, "[[", 3)
Smoothed_PerState_Strict_Masks        <- lapply(tempforeach, "[[", 4)

Smoothed_PerState_Loose_All           <- lapply(tempforeach, "[[", 5)
Smoothed_PerState_Loose_StayAtHome    <- lapply(tempforeach, "[[", 6)
Smoothed_PerState_Loose_TravelBan     <- lapply(tempforeach, "[[", 7)
Smoothed_PerState_Loose_Masks         <- lapply(tempforeach, "[[", 8)


names(Smoothed_PerState_Strict_All)           <- 
  names(Smoothed_PerState_Strict_StayAtHome)  <- 
  names(Smoothed_PerState_Strict_TravelBan)   <- 
  names(Smoothed_PerState_Strict_Masks)       <- 
  names(Smoothed_PerState_Loose_All)          <- 
  names(Smoothed_PerState_Loose_StayAtHome)   <- 
  names(Smoothed_PerState_Loose_TravelBan)    <- 
  names(Smoothed_PerState_Loose_Masks)        <- 
  vector.of.states

ordered.state.names <- statenames_here <- 
  c("Alaska", "Alabama", "Arkansas", "Arizona", "California", 
    "Colorado", "Connecticut", "D. C.", "Delaware", "Florida", 
    "Georgia", "Hawaii", "Iowa", "Idaho", "Illinois", "Indiana", 
    "Kansas", "Kentucky", "Louisiana", "Massachusetts", "Maryland", 
    "Maine", "Michigan", "Minnesota", "Missouri", "Mississippi", 
    "Montana", "North Carolina", "North Dakota", "Nebraska", 
    "New Hampshire", "New Jersey", "New Mexico", "Nevada", 
    "New York", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", 
    "Rhode Island", "South Carolina", "South Dakota", "Tennessee", 
    "Texas", "Utah", "Virginia", "Vermont", "Washington", "Wisconsin", 
    "West Virginia", "Wyoming")

#===============================================================
# PLOT RESULTS FOR AGGREGATE BARCHARTS
#===============================================================
Dates.picked  <- "2020-09-30"
last.date     <- which(Dates.selec == Dates.picked)

# STAY AT HOME POLICIES
#----------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_StayAtHome$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot1.perstate) <- ordered.state.names
Plot1.perstate        <- sort(Plot1.perstate, index.return = T)
write.csv(Plot1.perstate$x, "cf_individual_abs_StayAtHome_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("A state adopts a stay-at-home order")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.perstate$x) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x > 0))
                  ))
grid()
index.neg <- Plot1.perstate$x < 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}

barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_StayAtHome$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_abs_StayAtHome_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)


#---------------------
# MASK POLICIES
#---------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_Masks$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot1.perstate)               <- ordered.state.names
Plot1.perstate[Plot1.perstate > 0]  <- 0
Plot1.perstate                      <- sort(Plot1.perstate, index.return = T)

write.csv(Plot1.perstate$x, "cf_individual_abs_Masks_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("A state adopts a mask mandate")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.perstate$x) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x > 0))
                  ))
grid()
index.neg <- Plot1.perstate$x <= 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_Masks$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_abs_Masks_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)



#------------------------
# TRAVELBAN POLICIES
#------------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_TravelBan$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot1.perstate) <- ordered.state.names
Plot1.perstate        <- sort(Plot1.perstate, index.return = T)
write.csv(Plot1.perstate$x, "cf_individual_abs_TravelBan_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("A state bans crossings of its borders")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.perstate$x) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x > 0))
                  ))
grid()
index.neg <- Plot1.perstate$x < 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_TravelBan$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_abs_TravelBan_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.40 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)



#-----------------
# ALL POLICIES
#-----------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_All$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot1.perstate) <- ordered.state.names
Plot1.perstate        <- sort(Plot1.perstate, index.return = T)
write.csv(Plot1.perstate$x, "cf_individual_abs_All_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("A state adopts all three policies")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot1.perstate$x) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x > 0))
                  ))
grid()
index.neg <- Plot1.perstate$x < 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_All$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date] -  Smoothed.data[i, 1:last.date])[last.date] 
                       - (Smoothed.data[(5 * nb_states + i), 1] -  Smoothed.data[i, 1]))
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_abs_All_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Loose scenario")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)




###################################################
###################################################
# AGGREGATE BARCHARTS RELATIVE TO BASELINE
###################################################
###################################################

# STAY AT HOME POLICIES
#----------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_StayAtHome$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot1.perstate) <- ordered.state.names
Plot1.perstate        <- sort(Plot1.perstate, index.return = T)
write.csv(Plot1.perstate$x, "cf_individual_rel_StayAtHome_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("Stay-at-home")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(c(Plot1.perstate$x)) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x > 0))
                  ))
grid()
index.neg <- Plot1.perstate$x < 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_StayAtHome$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_rel_StayAtHome_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("No policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)




#-----------------
# MASKS POLICIES 
#-----------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_Masks$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot1.perstate) <- ordered.state.names
Plot1.perstate[Plot1.perstate > 0] <- 0
Plot1.perstate        <- sort(Plot1.perstate, index.return = T)
write.csv(Plot1.perstate$x, "cf_individual_rel_Masks_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("Stay-at-home")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(c(Plot1.perstate$x)) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x > 0))
                  ))
grid()
index.neg <- Plot1.perstate$x < 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_Masks$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_rel_Masks_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("No policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.40 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)




#---------------------
# TRAVELBAN POLICIES
#---------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_TravelBan$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot1.perstate) <- ordered.state.names
Plot1.perstate        <- sort(Plot1.perstate, index.return = T)
write.csv(Plot1.perstate$x, "cf_individual_rel_TravelBan_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("Stay-at-home")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(c(Plot1.perstate$x)) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x > 0))
                  ))
grid()
index.neg <- Plot1.perstate$x < 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_TravelBan$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_rel_TravelBan_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("No policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.40 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)


#---------------------
# ALL POLICIES
#---------------------
par(mfcol=c(1,2), mai = c(.5,1,.3,.1))

# STRICT SCENARIO
Plot1.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Strict_All$", vector.of.states[Selec.state[i]])))
  Plot1.perstate  <- c(Plot1.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot1.perstate) <- ordered.state.names
Plot1.perstate        <- sort(Plot1.perstate, index.return = T)
write.csv(Plot1.perstate$x, "cf_individual_rel_All_strict.csv")

midpts <- barplot(Plot1.perstate$x, 
                  names.arg = ordered.state.names[Plot1.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("Strict policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("Stay-at-home")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(c(Plot1.perstate$x)) + c(
                    -max(max(abs(Plot1.perstate$x))*.35 * (Plot1.perstate$x < 0)),
                    min(abs(range(c(Plot1.perstate$x))))
                  ))
grid()
index.neg <- Plot1.perstate$x < 0
index.pos <- Plot1.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot1.perstate$x[index.neg]),1, min) - 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot1.perstate$x[index.pos]),1, max) + 
         max(abs(Plot1.perstate$x))*0 , 
       labels = ordered.state.names[Plot1.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot1.perstate$x, 
        names.arg = vector.of.states[Plot1.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#13501354"), border = c("#135013"), add = T)


# LOOSE SCENARIO
Plot2.perstate <- NULL
for(i in 1:nb_states){
  eval(str2expression(paste0("Smoothed.data <- Smoothed_PerState_Loose_All$", vector.of.states[Selec.state[i]])))
  Plot2.perstate  <- c(Plot2.perstate, 
                       (cumsum(Smoothed.data[(5 * nb_states + i), 1:last.date])[last.date]  - Smoothed.data[(5 * nb_states + i), 1])
                       /(cumsum(Smoothed.data[i, 1:last.date])[last.date] -  Smoothed.data[i, 1]) - 1)
}
names(Plot2.perstate) <- ordered.state.names
Plot2.perstate        <- sort(Plot2.perstate, index.return = T)
write.csv(Plot2.perstate$x, "cf_individual_rel_All_loose.csv")

midpts <- barplot(Plot2.perstate$x, 
                  names.arg = ordered.state.names[Plot2.perstate$ix], 
                  cex.names = .5, beside = T, horiz = T,
                  main = expression(bold("No policies")),
                  col = rgb(0,0,0,alpha = 0), border = rgb(0,0,0,alpha = 0),
                  ylab = expression(bold("")), cex.lab = 1.2 , yaxt = "n",
                  xlim = range(Plot2.perstate$x) + c(
                    -max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x < 0)),
                    max(max(abs(Plot2.perstate$x))*.35 * (Plot2.perstate$x > 0))
                  ))
grid()
index.neg <- Plot2.perstate$x < 0
index.pos <- Plot2.perstate$x > 0
if(length(which(index.neg == T)) > 0) {
  text(y = midpts[index.neg], x = apply(cbind(Plot2.perstate$x[index.neg]),1, min) - 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.neg], pos = 2,  cex = .7)
}
if(length(which(index.pos == T)) > 0) {
  text(y = midpts[index.pos], x = apply(cbind(Plot2.perstate$x[index.pos]),1, max) + 
         max(abs(Plot2.perstate$x))*0 , 
       labels = ordered.state.names[Plot2.perstate$ix][index.pos], pos = 4, cex = .7)
}
barplot(Plot2.perstate$x, 
        names.arg = vector.of.states[Plot2.perstate$ix], xaxt = "n", yaxt = "n", 
        cex.names = .5, beside = T, horiz = T,
        col = c("#92000054"), border = c("#920000"), add = T)





