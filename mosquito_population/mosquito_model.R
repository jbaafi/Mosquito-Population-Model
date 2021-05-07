# Created by: Joseph Baafi
# Date: April 29, 2021
# Purpose: To predict mosquito population dynamics in seasonal environment

#Clear workspace
rm(list = ls())

#set working directory
setwd("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population")

# Import packages into r
library(tidyverse)
library(dplyr)
library(chron)
library(ggplot2)
library(deSolve)

#Import 2011-2016 climate data into r 
climate_11 <- read.csv("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population/data/en_climate_daily_ON_6153301_2011_P1D.csv", header = TRUE)
climate_12 <- read.csv("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population/data/en_climate_daily_ON_6153301_2012_P1D.csv", header = TRUE)
climate_13 <- read.csv("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population/data/en_climate_daily_ON_6153301_2013_P1D.csv", header = TRUE)
climate_14 <- read.csv("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population/data/en_climate_daily_ON_6153301_2014_P1D.csv", header = TRUE)
climate_15 <- read.csv("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population/data/en_climate_daily_ON_6153301_2015_P1D.csv", header = TRUE)
climate_16 <- read.csv("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population/data/en_climate_daily_ON_6153301_2016_P1D.csv", header = TRUE)

# To subset the data for the needed columns
climate_df_11 <- climate_11 %>% 
  select(Station.Name, Year:Day, Mean.Temp...C., Total.Precip..mm.)

climate_df_12 <- climate_12 %>% 
  select(Station.Name, Year:Day, Mean.Temp...C., Total.Precip..mm.)

climate_df_13 <- climate_13 %>% 
  select(Station.Name, Year:Day, Mean.Temp...C., Total.Precip..mm.)

climate_df_14 <- climate_14 %>% 
  select(Station.Name, Year:Day, Mean.Temp...C., Total.Precip..mm.)

climate_df_15 <- climate_15 %>% 
  select(Station.Name, Year:Day, Mean.Temp...C., Total.Precip..mm.)

climate_df_16 <- climate_16 %>% 
  select(Station.Name, Year:Day, Mean.Temp...C., Total.Precip..mm.)

# Combine the above data-frames into one data-frame.
climate_df <- rbind(climate_df_11, climate_df_12, climate_df_13, climate_df_14,
                    climate_df_15, climate_df_16) 

# Rename some of the columns of the climate_df to make it simple
climate.df <- climate_df %>% 
  rename(Mean.Temp = Mean.Temp...C.,
         Total.Precip = Total.Precip..mm.
  )

# Formatting time data in a form that chron can understand 
daily_dates <- dates(paste(climate.df$Month, 
                           climate.df$Day, 
                           climate.df$Year, sep="/")) 
climate.df$Chron.Date <-chron(dates=daily_dates, 
                              origin. = c(month = 1,day = 1,
                                          year = 2011)) # setting date into chron
climate.df <- arrange(climate.df, Chron.Date) # ordering dates and adding it to data-frame

start.date <- as.Date("2011-01-01") # The starting date
climate.df$Days.Since.Origin <- (as.numeric(as.Date(climate.df$Chron.Date) - start.date)) # Producing a sequence from smallest to largest

# counting NANs column-wise
sapply(climate.df, function(x) sum(is.na(x)))

# Formating data to exclude the NANs 
climate.df<- filter(climate.df, !is.na(Mean.Temp), !is.na(Total.Precip)) 

# Time series plot of temperature 
ggplot(data = climate.df, mapping = aes(x=Chron.Date, y=Mean.Temp))+
  geom_line(colour = "steelblue")+
  theme_classic()+
  xlab("Time (Days)")+ 
  ylab("Mean Temperature")+
  theme(axis.text.x=element_text(angle=45, hjust=1))

write.csv(climate.df, file = "climate.df.csv") 

# Time series plot of total precipitation
ggplot(data = climate.df, mapping = aes(x=Chron.Date, y=Total.Precip))+
  geom_line(colour = "steelblue")+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Total Precipitation") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Fitting function to data of mean temperature
Time <- climate.df$Days.Since.Origin
Temp <- climate.df$Mean.Temp

#Define a periodic function to fit to data
Temp.func <- sin(2*pi*Time/365) + cos(2*pi*Time/365)

#fsin <- sin(2*pi*Time/365)
#fcos <- cos(2*pi*Time/365)

fit.lm <- lm(Temp~Temp.func)

fit <- fitted(fit.lm)  

# find predictions for original time series
pred <- predict(fit.lm, newdata=data.frame(Time=Time))

#Plotting the data with fitted function with the base plot function
plot(Temp ~ Days.Since.Origin, data= climate.df)
lines(fit, col="red")
lines(Time, pred, col="blue")
summary(fit.lm)

#Plotting with same data as above with ggplot2
ggplot(data = climate.df)+
  geom_point(mapping = aes(x=Days.Since.Origin, y=Mean.Temp))+
  geom_line(mapping = aes(x=Days.Since.Origin, y=fit, colour="red"))+
  geom_line(mapping = aes(x=Time, y=pred, colour="blue"))

# Just to be sure of the function
# a <- mean(climate.df$Mean.Temp)
# b1 <- -8.7635
# f <- a + b1*sin(2*pi*Time/365) + b1*cos(2*pi*Time/365)
# plot(Temp~Days.Since.Origin, data = climate.df)
# lines(Time, f, col = "red")

# Fit a periodic function to total precipitation data
Precip <- climate.df$Total.Precip

#Define a periodic function to fit to data
Precip.func <- sin(2*pi*Time/365) + cos(2*pi*Time/365)

fitprecip.lm <- lm(Precip~Precip.func)

fitprecip <- fitted(fitprecip.lm)  

# find predictions for original time series
predprecip <- predict(fitprecip.lm, newdata=data.frame(Time=Time))

#Plotting the data with fitted function with the base plot function
plot(Precip ~ Days.Since.Origin, data= climate.df, ylim=c(0, 5))
lines(fitprecip, col="red")
lines(Time, predprecip, col="blue")
summary(fitprecip.lm)

################################################################################

#Model parameters
a <-  8.9231 #mean(climate.df$Mean.Temp)
b1 <- -8.7635
b2 <- b1

t_end <- 365
#t <- Time
#dtime = 0.5 # time steps 
t <-  seq(0, t_end, #dtime
          )
t0 <-t[1]

# Defining temperature function
 temp <- function(t){
   temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365)
   return(temp)
 }

 df <- data.frame(t, temp(t))
 
 #Plot of temperature as a function of time (365 days)
ggplot(data = df, mapping = aes(x=t, y=temp(t)))+
   geom_line()

# Fit scalar in gentrophic cycle rate function
alpha <- 168.0257252
beta <- 34.0632871
gamma <- 0.4371833 

# This function defines the relationship between temperature and genotrophic cycle
gen <- function(t){
  gen <-  gamma*exp(-(temp(t)-beta)^2/alpha)
  return(gen)
}
df1 <- data.frame(t, gen(t))

#plot of genotrophic cycle as a function of temperature as a function of time
ggplot(df1, aes(x=t, y=gen(t)))+
  geom_line()

# Fit scalars in egg development rate function
alpha_e <- 1337.666814
beta_e <- 75.098187
gamma_e <- 4.049403 

# This function defines the relationship between temperature and genotrophic cycle
egg.dev <- function(t){
  egg <-  gamma_e*exp(-(temp(t)-beta_e)^2/alpha_e)
  return(egg)
}
df2 <- data.frame(t, egg.dev(t))

#plot of egg development as a function of temperature as a function of time
ggplot(df2, aes(x=t, y=egg.dev(t)))+
  geom_line()
 
# Fit scalars in larvae development rate function
alpha_l <- 109.03255
beta_l <- 27.18118
gamma_l <- 0.16231

# This function defines the relationship between temperature and genotrophic cycle
larva.dev <- function(t){
  larva <-  gamma_l*exp(-(temp(t)-beta_l)^2/alpha_l)
  return(larva)
}
df3 <- data.frame(t, larva.dev(t))

#plot of larva development as a function of temperature as a function of time
ggplot(df3, aes(x=t, y=larva.dev(t)))+
  geom_line()

# Fit scalars in pupal development rate function
alpha_p <- 339.6702830
beta_p <- 39.9070020
gamma_p <-  0.5920232

# This function defines pupa development as a function of temperature
pupa.dev <- function(t){
  pupa <-  gamma_p*exp(-(temp(t)-beta_p)^2/alpha_p)
  return(pupa)
}

df4 <- data.frame(t, pupa.dev(t))

#plot of pupa development as a function of temperature as a function of time
ggplot(df4, aes(x=t, y=pupa.dev(t)))+
  geom_line()

# Fit scalars in egg mortality function
c_E <- 1.118e-03
T_E <- 2.218e+01
d_E <- 1.488e-02

# This function defines egg mortality as a function of temperature
egg.mortality <- function(t){
  egg.mortality <-  c_E*(temp(t) - T_E)^2 + d_E
  return(egg.mortality)
}

df5 <- data.frame(t, egg.mortality(t))

#plot of egg mortality as a function of temperature as a function of time
ggplot(df5, aes(x=t, y=egg.mortality(t)))+
  geom_line()


# Fit scalars in larvae mortality function
c_L <- 0.0025 #(Abdelrazek) #1.118e-03
T_L <- 20 #(Abdelrazak) #2.218e+01
d_L <- 0.2 #1.488e-02

# This function defines larva mortality as a function of temperature
larva.mortality <- function(t){
  larva.mortality <-  c_L*(temp(t) - T_L)^2 + d_L
  return(larva.mortality)
}

df6 <- data.frame(t, larva.mortality(t))

#plot of larva mortality as a function of temperature as a function of time
ggplot(df6, aes(x=t, y=larva.mortality(t)))+
  geom_line()

# Fit scalars in pupa mortality function
c_P <- 0.001 #(Abdelrazek) #1.118e-03
T_P <- 20 #(Abdelrazak) #2.218e+01
d_P <- 0.15 #1.488e-02

# This function defines pupa mortality as a function of temperature
pupa.mortality <- function(t){
  pupa.mortality <-  c_P*(temp(t) - T_P)^2 + d_P
  return(pupa.mortality)
}

df7 <- data.frame(t, pupa.mortality(t))

#plot of pupa mortality as a function of temperature as a function of time
ggplot(df7, aes(x=t, y=pupa.mortality(t)))+
  geom_line()

# Fit scalars in adult mortality function
c_A <- 0.08841 
T_A <- 21.24746
d_A <- 14.92552

# This function defines adult mortality as a function of temperature
adult.mortality <- function(t){
  adult.mortality <-  c_A*exp(((temp(t)-T_A)/d_A)^4)
  return(adult.mortality)
}

df8 <- data.frame(t, adult.mortality(t))

#plot of adult mortality as a function of temperature as a function of time
ggplot(df8, aes(x=t, y=adult.mortality(t)))+
  geom_line()

# Other model parameters [Values obtained from Hamdan and Kilicman, 2020]
mu_E <- 0.47
mu_L <- 0.47
mu_P <- 0.47
mu_A <- 0.1
tau <- 0.5 # fraction of mosquitoes that emerge as adult females. 
alpha_b <- 300 # Maximum number of eggs laid per oviposition [value taken from Abdelrazec & Gumel]
# Ignore density-dependent mortality for now.

# Function for the system of equations
model <- function(t, y){
  # The number of eggs
  E = y[1]
  # The number of Larvae
  L = y[2]
  # The number of Pupa
  P = y[3]
  # The number of Matured (Adult) mosquitoes
  M = y[4]
  
  # The system of equations
  dE <- alpha_b*gen(t)*M - (egg.dev(t) + egg.mortality(t) + mu_E)*E
  dL <- egg.dev(t)*E - (larva.dev(t) + larva.mortality(t) + mu_L)*L
  dP <- larva.dev(t)*L - (pupa.dev(t) + pupa.mortality(t) + mu_P)*P
  dM <- tau*pupa.dev(t)*P - (adult.mortality(t) + mu_A)*M
  return(c(dE, dL, dP, dM))
  
}

#Initial conditions
y0 <- c(100, 0, 0, 10)

# Numerical integration. 
out <-  ode(y = y0, func = model, times = t, parms = )
out <- data.frame(out)

head(out)



# Define other juvinile mortality functions (all the tree juvinile stages can assume same function)
#and density dependent mortality function (from Beck-Johnson et al, 2013)
# Find natural mortality rates from Hamdan et al, to run the model. 






