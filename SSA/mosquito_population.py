# coding=utf-8
# want to use this code to perform SSA of a mosquito population dynamics model

# Importing packages
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import random
from scipy.integrate import odeint

# Define initial states
E = [5]
L = [0]
P = [0]
A = [0]

# Initial time defined
t = [0]

#Time to end the simulation
tend = 15

################################################################
#Parameter values
# Fit scalar in gentrophic cycle rate function
alpha = 250.0257252
beta = 34.0632871
gamma = 5.3371833

#genotrophic cycle as a function of time
def gen_cycle(t):
    func = gamma*np.exp(-(t[-1]-beta)**2/alpha)
    return(func)

# Fit scalars in egg development rate function
alpha_e = 1317.666814
beta_e = 75.098187
gamma_e = 20.049403

# egg development as a function of time
def egg_dev(t):
    egg = gamma_e*np.exp(-(t[-1]-beta_e)**2/alpha_e)
    return(egg)

##############################################################

alpha_b = 100    # avarage number of eggs laid per oviposition (day^-1) (source: CDC)
phi = gen_cycle(t) #1/4 #0.051      # oviposition rate per adult mosquito (source: https://www.pascomosquito.org/resources/mosquito-biology/, https://megacatch.com/mosquito-lifecycle-faqs/)
delta_E = egg_dev(t) #1/2    # hatching rate of eggs into larvae (per day)
delta_L = 1/5    # development rate of larvae into pupae (per day)
delta_P = 1/2    # development rate of pupae into adult (per day)
k = 1/2          # fraction of juviniles becoming adult females
mu_E = 0.01      # egg mortality rate
mu_L = 0.02      # larva mortality rate
mu_P = 0.01      # pupal mortality rate
mu_A = 1/42   # mortality rate of adult mosquitoes (source: https://www.vdci.net/mosquito-biology-101-life-cycle/)

# The Gillespie Algorithm
# This stochastic process is a Continuous Time Markov Chain.

while t[-1] < tend:

    rates = [alpha_b * phi * A[-1], delta_E * E[-1], delta_L * L[-1], k * delta_P * P[-1], \
             mu_E * E[-1], mu_L * L[-1], mu_P * P[-1], mu_A * A[-1]]

    rate_sum = sum(rates)
    if rate_sum == 0:
        break

    tau = np.random.exponential(scale=1 / rate_sum)

    t.append(t[-1] + tau)

    rand = random.uniform(0, 1)

    # Egg oviposition event
    if rand * rate_sum <= rates[0]:

        E.append(E[-1] + 1)
        L.append(L[-1])
        P.append(P[-1])
        A.append(A[-1])

    # Hatching of eggs and development of larvae event
    elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):

        E.append(E[-1] - 1)
        L.append(L[-1] + 1)
        P.append(P[-1])
        A.append(A[-1])

    # Development of larvae into pupae event
    elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):

        E.append(E[-1])
        L.append(L[-1] - 1)
        P.append(P[-1] + 1)
        A.append(A[-1])

    # Development of pupa into adults event
    elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):

        E.append(E[-1])
        L.append(L[-1])
        P.append(P[-1] - 1)
        A.append(A[-1] + 1)

    # Egg mortality event
    elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):

        E.append(E[-1] - 1)
        L.append(L[-1])
        P.append(P[-1])
        A.append(A[-1])

    # Larvae mortality event
    elif rand * rate_sum > sum(rates[:5]) and rand * rate_sum <= sum(rates[:6]):

        E.append(E[-1])
        L.append(L[-1] - 1)
        P.append(P[-1])
        A.append(A[-1])


    # Pupae mortality event
    elif rand * rate_sum > sum(rates[:6]) and rand * rate_sum <= sum(rates[:7]):

        E.append(E[-1])
        L.append(L[-1])
        P.append(P[-1] - 1)
        A.append(A[-1])

    # Adult mortality event
    elif rand * rate_sum > sum(rates[:7]) and rand * rate_sum <= sum(rates[:8]):

        E.append(E[-1])
        L.append(L[-1])
        P.append(P[-1])
        A.append(A[-1] - 1)

E_plot, = plt.plot(t,E, label="E")
L_plot, = plt.plot(t,L, label="L")
P_plot, = plt.plot(t,P, label="P")
A_plot, = plt.plot(t,A, label="A")

plt.legend(handles=[E_plot, L_plot, P_plot, A_plot])
plt.xlabel("Time")
plt.ylabel("Abundance")
plt.show()