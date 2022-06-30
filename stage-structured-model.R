#created by: Joseph Baafi
#Date: November 10, 2020
#Purpose: To analyze mosquito population dynamics (autonomous model)

#Clear workspace
rm(list = ls())

#Load R package for solving diff equations
library(deSolve)
library(ggplot2)
library(scales)


#Set working directory 
setwd("/Users/jbaafi/Desktop/Mosquito population dynamics")

# Function to return rate of change of the states in the ODE model
model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), { 
    dE<- b*phi*(1-(E/K))*A - (deltaE+muE+pieE)*E
    dL<- deltaE*E - (deltaL+muL+pieL + sigmaL*L)*L
    dP <- deltaL*L - (deltaP+muP+pieP)*P
    dA <- chi*deltaP*P - (muA+pieA)*A
    list(c(dE, dL, dP, dA))
  })
}

# Set of parameter values used in the model
parameters <- c(
  b = 300,
  phi = 3,
  K = 10^3,
  deltaE = 0.55,
  deltaL = 0.20,
  deltaP = 0.3,
  muE =  0.3,#0.0143 # value taken from Hamdan & Kilicman, 2020
  muL = 0.30, #0.0143 # value taken from Hamdan & Kilicman, 2020
  muP = 0.15, #0.0143 # value taken from Hamdan & Kilicman, 2020
  muA = 0.10,   # value taken from Hamdan & Kilicman, 2020
  pieE = 0.011,
  pieL = 0.10,
  pieP =  0.01,
  pieA = 0.07,
  sigmaL = 0.004,
  chi = 0.7     # value taken from Hamdan & Kilicman, 2020
) 

# Initial conditions
state <- c(
  E = 100, 
  L = 0, 
  P = 0,
  A = 0)

# Time steps in running the model
times <- seq(0, 200, by = 0.01)

# solving the system with ode function from the deSolve package
out <- ode(y = state, times = times, func = model, parms = parameters)

# Converting the numerical solution to a Dataframe
out <- data.frame(out)
head(out)

#plotting the numrical results 
plot(out$time, out$A, type = "l", col = "blue", ylim = c(0, 10^3), xlab = "Time (in days)", ylab = "Population")
lines(out$time, out$E, col = "green")
lines(out$time, out$L, col = "yellow")
lines(out$time, out$P, col = "red")
legend( "topright", c("A(t)", "E(t)", "L(t)", "P(t)"), 
        text.col=c("blue", "green", "yellow","red") )


# A few math expressions to produce a good plot
out$Egg <- out$E/10^3
out$Larvae <- out$L/10^3
out$Pupae <- out$P/10^3
out$Adult <- out$A/10^3

plot(out$time, out$Egg, type = "l", col = "green", ylab = "Density (in thousands)")
lines(out$time, out$Larvae, col = "red")
lines(out$time, out$Pupae, col = "blue")
lines(out$time, out$Adult, col = "yellow")
legend( "right", c("Eggs", "Larvae", "Pupae", "Adults"), 
        text.col=c("green", "red", "blue","yellow") )


ggplot(out,aes(x=time, y=E)) + 
  geom_line() +
  scale_y_continuous(labels = label_number(suffix = " T", scale = 1e-3))



##################################################################################################
# This lines of code helps to find the basic reproduction number of the model above
phi = 30
K = 10^6
deltaE = 0.4
deltaL = 0.14
deltaP = 0.3
muE =  0.36#0.0143 # value taken from Hamdan & Kilicman, 2020
muL = 0.30 #0.0143 # value taken from Hamdan & Kilicman, 2020
muP = 0.15 #0.0143 # value taken from Hamdan & Kilicman, 2020
muA = 0.10   # value taken from Hamdan & Kilicman, 2020
pieE = 0.11
pieL = 0.10
pieP =  0.01
pieA = 0.07
sigmaL = 0.004
chi = 0.7    # value taken from Hamdan & Kilicman, 2020
# The basic reproduction number  
Rm <- (chi*phi*deltaP*deltaL*deltaP)/((muA+pieA)*(deltaP+muP+pieP)*(deltaL+muL+pieL + sigmaL)*(deltaE+muE+pieE))

print(Rm)

####################################################################################################3##
# Create a vector of temperature values over time 
#time = c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec")
time <- seq(1, 12)
Temperature = c(28.3, 28.1, 29.3, 28.3, 28.6, 28.4, 28.0, 28.5, 28.0, 27.9, 27.0, 27.3)
phi.values <- c(8.294, 8.294, 8.704, 8.294, 8.704, 8.295, 8.294, 8.704, 7.741, 7.082, 7.084, 7.084)
data <- data.frame(time, Temperature, phi.values)

max <- max(out$E)
df <- out$E/max
df

out$new_col <- df
