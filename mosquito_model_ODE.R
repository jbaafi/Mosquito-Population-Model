# This code solves the Abdelrazac & Gumel 2017 model numerically using PBSddesolve package in R.

# clears the workspace
rm(list = ls())

#Load R package for solving DEs
require(PBSddesolve)

# Set working directory for the script
setwd("/Users/jbaafi/Desktop/R-codes/mosquito_population/modeling_project")

#Parameters for seasonal temperature and rainfall
a <- 20
 b <-  0.3
 c <- 20
 d <- 10
 e <- 1
 t <- seq(0, 360, 1)
 
 #Function for temperature
 T <- function(t){
   temp <- a + b*sin(2*pi * t/360)
   return(temp)
 }

 # t <- seq(0, 365, 1)
 # a <- 0.00015
 # Tmax <- 365/2
 # i <- 25
 # 
 # T <- function(t){
 #   Temp<- i*exp(-a*(t-Tmax)^2)
 #   return(Temp)
 # }
 # 
 # Rmax <- 45
 # R <- function(t){
 #   rain<- i*exp(-a*(t-Rmax)^2)
 #   return(rain)
 # }
 # 
 
 #Function for rainfall
 R <- function(t){
   rain <- c + d*sin(2*pi * t/365)
   return(rain)
 }
 
# I will now define the various temperature dependent parameter functions. 
u_b <- function(t){
  ubT <- exp(-a_b*(T(t)-T_b)^2)
  return(ubT)
}

mu_M <- function(t){
  muMT <- c_M*(T(t)-T_Mstar)^2 + d_M
  return(muMT)
}

# Effect of temperature on the hatching rate
g_E <- function(t){
  gE <- exp(-a_E*(T(t)-T_E)^2)
  return(gE)
}

# Effect of temperature on transition rate F_L
g_L <- function(t){
  gL <- exp(-a_L*(T(t)-T_L)^2)
  return(gL)
}
# Effect of temperature on the transition rate F_P
g_P <- function(t){
  gP <- exp(-a_P*(T(t)-T_P)^2)
  return(gP)
}

# Effect of temperature on mortality rate of eggs
p_E <- function(t){
  pE <- c_E*(T(t)-T_Estar)^2 + d_E
  return(pE)
}

# Effect of temperature on mortality rate of larvae
p_L <- function(t){
  pL <- c_L*(T(t)-T_Lstar)^2 + d_L
  return(pL)
}

# Effect of temperature on mortality rate of pupae
p_P <- function(t){
  pP <- c_P*(T(t)-T_Pstar)^2 + d_P
  return(pP)
}

#  The following functions define the various rainfall dependent parameters 
v_b <- function(t){
  vb <- ((1+s_b)*exp(-r_b*(R(t)-R_b)^2))/(exp(-r_b*(R(t)-R_b)^2) + s_b)
  return(vb)
}

# Effect of rainfall on hatching rate of eggs
h_E <- function(t){
  hE <- ((1+s_E)*exp(-r_E*(R(t)-R_E)^2))/(exp(-r_E*(R(t)-R_E)^2) + s_E)
  return(hE)
}

# Effect of rainfall on transition rate from eggs to larvae
h_L <- function(t){
  hL <- ((1+s_L)*exp(-r_L*(R(t)-R_L)^2))/(exp(-r_L*(R(t)-R_L)^2) + s_L)
  return(hL)
}

# Effect of rainfall on transition rate from larvae to pupae
h_P <- function(t){
  hP <- ((1+s_P)*exp(-r_P*(R(t)-R_P)^2))/(exp(-r_P*(R(t)-R_P)^2) + s_P)
  return(hP)
}

# Effect of rainfall on mortality rate of eggs
q_E <- function(t){
  qE <- 1 + (e_E*R(t))/(1+R(t))
  return(qE)
}

# Effect of rainfall on mortality rate of larvae
q_L <- function(t){
  qL <- 1 + (e_L*R(t))/(1+R(t))
  return(qL)
}

# Effect of rainfall on mortality rate of pupae
q_P <- function(t){
  qP <- 1 + (e_P*R(t))/(1+R(t))
  return(qP)
}

# Natural mortality rate of eggs
mu_E <- function(t){
  muE <- p_E(T(t))*q_E(R(t))
  return(muE)
}
# Natural mortality rate of larvae
mu_L <- function(t){
  muL <- p_L(T(t))*q_L(R(t))
  return(muL)
}

# Natural mortality rate of pupae
mu_P <- function(t){
  muP <- p_P(T(t))*q_P(R(t))
  return(muP)
}

# Hatching rate of eggs
F_E <- function(t){
  FE <- alpha_E*g_E(T(t))*h_E(R(t))
  return(FE)
}

# Development rate of larvae into pupae
F_L <- function(t){
  FL <- alpha_L*g_L(T(t))*h_L(R(t))
  return(FL)
}

# Development rate of pupae into adult mosquitoes
F_P <- function(t){
  FP <- alpha_P*g_P(T(t))*h_P(R(t))
  return(FP)
}

# Dfine the parameter b(T, R), which is the eggs oviposition rate
ovipositor <- function(t){
  bTR <- alpha_b*u_b(T(t))*v_b(R(t))
  return(bTR)
}

# The Verhulst-Pearl logistic oviposition function
B <- function(t, M){
  ovip <- ovipositor(t)*(1-M/K)
  return(ovip)
}

# Here are the parameter values used in the model. 
alpha_b <- 300
alpha_P <- 0.5
T_Pstar <- 20
T_L <- 22
c_E <- 0.001
c_M <- 0.0005
d_P <- 0.15
a_E <- 0.011
e_E <- 1.1
e_M <- 0
s_L <- 1.5
r_E <- 0.05
R_b <- 10
alpha_E <- 0.5
T_Estar <- 20
T_Mstar <- 28
T_P <- 22
c_L <- 0.0025
d_E <- 0.15
d_M <- 0.04
a_L <- 0.013
e_L <- 1.1
s_b <- 1.2
s_P <- 1.5
r_L <- 0.05
R_E <- 15
alpha_L <- 0.35
T_Lstar <- 20
T_E <- 22
T_M <- 27
c_P <- 0.001
d_L <- 0.2
a_b <- 0.015
a_P <- 0.014
e_P <- 1.1
s_E <- 1.5
r_b <- 0.05
r_P <- 0.05
R_L <- 15

T_b <- 22
K <- 10^6
delta_L <- 0.2
R_P <- 15
sigma <- 0.5

# Function for the system of equations
mosquito <- function(t, y){
  # The number of eggs
  E = y[1]
  # The number of Larvae
  L = y[2]
  # The number of Pupa
  P = y[3]
  # The number of Matured (Adult) mosquitoes
  M = y[4]
  
  # The system of equatiions
  dE <- M*B(t, M) - (F_E(t) + mu_E(t))*E
  dL <- F_E(t)*E - (F_L(t) + mu_L(t) + delta_L*L)*L
  dP <- F_L(t)*L - (F_P(t) + mu_P(t))*P
  dM <- sigma*F_P(t)*P - mu_M(t)*M
  return(c(dE, dL, dP, dM))
  
  # Parameter values 
  # parameters = c(alpha_b, alpha_P, T_Pstar, T_L, c_E, c_M, d_P, a_E, e_E, e_M, s_L, 
  #                r_E, R_b, alpha_E, T_Estar, T_Mstar, T_P, c_L, d_E, d_M, a_L, e_L, 
  #                s_b, s_P, r_L, R_E, alpha_L, T_Lstar, T_E, T_M, c_P, d_L, a_b, a_P,
  #                e_P, s_E, r_b, r_P, R_L, T_b)
  

}

#Initial conditions
y0 <- c(10, 0, 0, 20)

# Numerical integration. 
x <-  dde(y = y0, func = mosquito, times = t, hbsize = 0)
x <- data.frame(x)

head(x)

# Plots of variables against time.
plot(t, x$y1, type = "l", col = "orange", xlab = "Time", ylab = "States Dynamics")
lines(t, x$y2, col="red")
lines(t, x$y3, col="green")
lines(t, x$y4, col="blue")
legend( "topright", c("E", "L", "P", "M"), 
        text.col=c("orange", "red", "green","blue") )
