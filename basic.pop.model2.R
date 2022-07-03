# Title: To solve stage-structured ODE model with environmental carrying capacity of eggs
# Date created: June 24, 2022
# Author: Joseph Baafi

# Clear workspace
rm(list = ls())

# set working directory
#setwd("/Users/jbaafi/Desktop/Mosquito population dynamics")

# Load packages
pacman::p_load(pacman, deSolve, tidyverse)

# Define a function for the system of equations
pop.model <- function(t, y, ...){
  # The number of states
  E = y[1]
  L = y[2]
  P = y[3]
  A = y[4]
  # The system of equations (This model is based on Lutambi et al, 2013)
  #dE <- b*phiA*(1-E/k)*A - (FE + muE)*E
  dE <- b*(1-E/k)*A - FE*E - muE*E
  dL <- FE*E - (FL + muL + deltaL*L)*L
  dP <- FL*L - (FP + muP)*P
  dA <- sigma*FP*P - (muA)*A
  return(list(c(dE, dL, dP, dA)))
}

# Model parameters (All parameter values are based on Lutambi et al, 2013)
b  = 200      # number of eggs laid per oviposition
phiA = 3.00   # oviposition rate 
k = 10^4     # environment carrying capacity of eggs laid assuming there is not much waters
FE = 0.50     # egg hatching takes 2 days on average
muE = 0.56    # 
FL = 0.14    # the larval period for the Anopheles mosquitoes is found to be 7 days
muL = 0.44
deltaL = 0.05
FP = 0.50    # the pupal period (1/FP) lasts for 2 days on average
muP = 0.37
sigma = 0.5  # Fraction of larvae that emerge as female mosquitoes
muA = 0.10 # 0.05 - 0.15 (This parameter is not clear with Lutambi so I can vary it)

param <- c(b, phiA, k, FE, muE, FL, muL, deltaL, FP, muP, sigma, muA)

# Vector reproduction numbers
R0 <- b*phiA*sigma*FE*FL*FP/((muE+FE)*(muL+FL)*(muP+FP)*muA)
print(paste0("The value of R0 is ", R0))


# Initial conditions (I am not sure what values to use as initial conditions)
y.initial <- c(0.0, 0.0, 0.0, 10.0)

# Time steps (Period with which to run the simulation)
times <- seq(0, 150, by = 0.01)

# Numerical integration using the ode() function from deSolve package
out <-  ode(y = y.initial, func = pop.model, times = times, parms = param)

# Put the simulated values into a dataframe
out <- data.frame(out)

# Set plot environment 
par(mfrow = c(2, 2), mar=c(2,4,2,2))    # 2 rows, 2 columns

# Create plots
plot(out$time, out$X1, type = "l", xlab = "Time (days)", ylab = "Egg density", 
     pch = 19, col = "blue", main = "Popltn Dynamics(Constant Param)")     # plot no. 1

plot(out$time, out$X2, type = "l", xlab = "Time (days)", ylab = "Larval density", 
     pch = 19, col = "green")     # plot no. 2

plot(out$time, out$X3, type = "l", xlab = "Time (days)", ylab = "Pupal density", 
     pch = 19, col = "yellow")     # plot no. 3

plot(out$time, out$X4, type = "l", xlab = "Time (days)", ylab = "Adult females", 
     pch = 19, col = "red")    # plot no. 2

# This is a model with carrying capacity of eggs. Since we have more than the carrying capacity 
# of eggs laid, I was expecting the time series plot of egg population to kind of oscillate around
# the carrying capacity before eventually reaching a constant equilibrium but I don't 
# see that in the plot. Is there any problem with the model? Was the carrying capacity 
# correctly introduced?


# The next file basic.pop.model3.R solves the model with extra mortality rates added to 
# the states. I also extend the basic model by dividing the Adult stage, A, into 3 states 
# ie. Adult seeking host to bite (Ah), Adult at rest to develop the eggs (Ar) and 
# Adult seeking for oviposition site (Ao). 
