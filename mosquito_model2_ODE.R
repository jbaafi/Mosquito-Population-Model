# This code solves the Abdelrazac & Gumel 2017 model numerically using PBSddesolve package in R.

# clears the workspace
rm(list = ls())

#Load R package for solving DEs
require(PBSddesolve)

# Set working directory for the script
setwd("/Users/jbaafi/Desktop/R-codes/mosquito_population/modeling_project")

a <-  1.5
b <-  1
c <- 1
d <- 2
e <- 1
t <- seq(0, 360, 1)
#Temperature and Rainfall depends on time
T <- function(t){
  temp <- a + b*sin( pi * t/360)
  return(temp)
}

R <- function(t){
  rain <- c + d*sin(pi * t/365) 
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

g_E <- function(t){
  gE <- exp(-a_E*(T(t)-T_E)^2)
  return(gE)
}

g_L <- function(t){
  gL <- exp(-a_L*(T(t)-T_L)^2)
  return(gL)
}

g_P <- function(t){
  gP <- exp(-a_P*(T(t)-T_P)^2)
  return(gP)
}

p_E <- function(t){
  pE <- c_E*(T(t)-T_Estar)^2 + d_E
  return(pE)
}

p_L <- function(t){
  pL <- c_L*(T(t)-T_Lstar)^2 + d_L
  return(pL)
}

p_P <- function(t){
  pP <- c_P*(T(t)-T_Pstar)^2 + d_P
  return(pP)
}

#  The following function define the various rainfall dependent parameter functions. 
v_b <- function(t){
  vb <- ((1+s_b)*exp(-r_b*(R(t)-R_b)^2))/(exp(-r_b*(R(t)-R_b)^2) + s_b)
  return(vb)
}

h_E <- function(t){
  hE <- ((1+s_E)*exp(-r_E*(R(t)-R_E)^2))/(exp(-r_E*(R(t)-R_E)^2) + s_E)
  return(hE)
}

h_L <- function(t){
  hL <- ((1+s_L)*exp(-r_L*(R(t)-R_L)^2))/(exp(-r_L*(R(t)-R_L)^2) + s_L)
  return(hL)
}

h_P <- function(t){
  hP <- ((1+s_P)*exp(-r_P*(R(t)-R_P)^2))/(exp(-r_P*(R(t)-R_P)^2) + s_P)
  return(hP)
}

q_E <- function(t){
  qE <- 1 + (e_E*R(t))/(1+R(t))
  return(qE)
}

q_L <- function(t){
  qL <- 1 + (e_L*R(t))/(1+R(t))
  return(qL)
}

q_P <- function(t){
  qP <- 1 + (e_P*R(t))/(1+R(t))
  return(qP)
}

# These lines of code defines the function for the mortality rate at the various life estages
mu_E <- function(t){
  muE <- p_E(T(t))*q_E(R(t))
  return(muE)
}

mu_L <- function(t){
  muL <- p_L(T(t))*q_L(R(t))
  return(muL)
}

mu_P <- function(t){
  muP <- p_P(T(t))*q_P(R(t))
  return(muP)
}

# Define the transition rates in the model. 
F_E <- function(t){
  FE <- alpha_E*g_E(T(t))*h_E(R(t))
  return(FE)
}

F_L <- function(t){
  FL <- alpha_L*g_L(T(t))*h_L(R(t))
  return(FL)
}

F_P <- function(t){
  FP <- alpha_P*g_P(T(t))*h_P(R(t))
  return(FP)
}

# Define the parameter b(T, R), for the rate of eggs laid per oviposition
b <- function(t){
  bTR <- alpha_b*u_b(T(t))*v_b(R(t))
  return(bTR)
}

# Define the Verhulst-Pearl logistic oviposition function as follows
B <- function(t, M){
  ovip <- b(T(t), R(t))*(1-M/K)
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

# Define the system of ODEs that controls the dynamics occurring in the various states 
mosquito <- function(t, y){
  # The number of eggs
  E = y[1]
  # The number of Larvae
  L = y[2]
  # The number of Pupa
  P = y[3]
  # The number of Matured (Adult) mosquitoes
  M = y[4]
  
  # The system of ODEs
  dE <- M*B(M) - (F_E(T(t), R(t)) + mu_E(T(t), R(t)))*E
  dL <- F_E(T(t), R(t))*E - (F_L(T(t), R(t)) + mu_L(T(t), R(t)) + delta_L*L)*L
  dP <- F_L(T(t), R(t))*L - (F_P(T(t), R(t)) + mu_P(T(t), R(t)))*P
  dM <- sigma*F_P(T(t), R(t))*P - mu_M(T(t))*M
  return(c(dE, dL, dP, dM))
  
  # Parameter values 
  # parameters = c(alpha_b, alpha_P, T_Pstar, T_L, c_E, c_M, d_P, a_E, e_E, e_M, s_L, 
  #                r_E, R_b, alpha_E, T_Estar, T_Mstar, T_P, c_L, d_E, d_M, a_L, e_L, 
  #                s_b, s_P, r_L, R_E, alpha_L, T_Lstar, T_E, T_M, c_P, d_L, a_b, a_P,
  #                e_P, s_E, r_b, r_P, R_L, T_b)
  
  #List of initial values
  y0 <- c(10, 0, 0, 0)
  
  # Define a function to solve the ODEs numerically. 
  x = dde(y = y0, func = mosquito, times = t, hbsize = 0)
  return(x)
}

head(x)

