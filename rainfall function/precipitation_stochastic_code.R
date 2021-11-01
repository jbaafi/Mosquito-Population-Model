# Purpose:    Code to fit a probability distribution function (PDF) to precipitation data
# Created on: Sept 15, 2021
# By:         jbaafi

# clear workspace
rm(list=ls())

# set working directory
setwd("/Users/jbaafi/Google Drive/My Drive/rainfall function")

# importing packages needed
packages <- c("tidyverse", "stringr", "dplyr", "base", 
      "ggplot2", "patchwork", "bbmle", "fitdistrplus")
lapply(packages, require, character.only = TRUE)

# reading dataset into r
precip <- read.csv("/Users/jbaafi/Google Drive/My Drive/rainfall function/climate.df.csv")

attach(precip)

# minimum precipitation recorded
min(precip$Total.Precip)

# maximum precipitation recorded
max(precip$Total.Precip)

# drawing a histogram!
hist(precip$Total.Precip, breaks = 50)

#counting the number of zero recordings in the total.precip data
precip[(precip$Total.Precip==0),]
sum(precip$Total.Precip==0)

# eliminating zero recordings 
df <- precip %>% 
  filter(Total.Precip != 0 & !is.na(Total.Precip) & Total.Precip <= 50) %>% 
  arrange(Total.Precip)

#ploting the histogram without zero recordings 
hist_precip <-  hist(df$Total.Precip, breaks = 100, plot = T)

plot(hist_precip$mids, hist_precip$counts, "l")

plot(hist_precip$mids, hist_precip$density, type = "l")

df2 <- data.frame(hist_precip$mids, hist_precip$counts, hist_precip$density)

df2 <- df2 %>% 
  rename(precip = hist_precip.mids, freq = hist_precip.counts, density = hist_precip.density)


################################################################################
x <- df2$precip
y <- df2$density

# Using the nls() function to fit exponential function to data
exp<- nls(y ~ lambda*exp(-lambda*x), data = df2, start = list(lambda = 0.67))

plot(x, y, col = "blue", ylim = c(0, 0.7))
lines(x,predict(exp),col="red")

summary(exp)
coef(exp)


# fitting the log-normal PDF to data by the nls() function.
logn <- nls(y ~ (1/(x*delta*sqrt(2*pi)))*exp(-((log(x)-mu)^2)/(2*delta^2)), data = df2, start = list(mu = 0.499, delta = 1.2))

plot(x, y, col = "blue", ylim = c(0, 0.7))
lines(x,predict(logn),col="red")

summary(logn)
coef(logn)

###############################################################################




lambda <- 0.674043   #0.749641 #1/mean(x)
f <- lambda*exp(-lambda*x)
plot(x, f1, "l", col = "red")
points(x, y, col = "blue")

fit.lm <- lm(y~f)
fit <- fitted(fit.lm) 
pred <- predict(fit.lm, newdata = data.frame(x=x))

#Plotting the data with fitted function with the base plot function
plot(y ~ x)
lines( x, fit, col="blue")
summary(fit.lm)

################################################################################
# Estimating the lambda parameter with the MLE
exponential <- function(lambda){
  d <- lambda*exp(-lambda*x)
}

expon.negLL = function(lambda){
  -sum(dlnorm(x, mean=exponential(lambda), sd=1, log=TRUE))
} 

mle.expon <- mle2(expon.negLL, start = list(lambda=0.67)) 

lambda <- unname(coef(mle.expon))
lambda

plot(y~x) # Probability
lines(x, exponential(lambda), col="blue", pch=20)
lines(x, fit, col="red")

#################################################################################


# Fitting a log-normal distribution
mu <-  0.2951328 #0.4994043   
delta <- 1.8758427 #1.2
logn <- (1/(x*delta*sqrt(2*pi)))*exp(-((log(x)-mu)^2)/(2*delta^2))
plot(x, y, col = "red")
lines(x, logn, col = "blue")


fit.logn <- lm(y~logn)
fit <- fitted(fit.logn) 
pred <- predict(fit.logn, newdata = data.frame(x=x))

#Plotting the data with fitted function with the base plot
plot(y ~ x)
lines( x, fit.logn, col="blue")
lines(x, pred)
summary(fit.lm)

######################################################################################


# The idea is to create a vector for the x-axis that will be similar to the 
#vector in the actual data.
dx <- 0.5 # you may change this.
sub.edge <- seq(1, 50, dx) # min max rainfall mm you wanna consider 
sub.mid <- seq(1+dx/2, 50-dx/2, dx) # finding the class midpoints 
sub.x <- sub.mid

# Subsetting to only take values that are greater or equal to 1
# we subset the actual data to set the x-axis similar to that 
# of the distribution done just above this.
sub.precip <- subset(precip, Total.Precip >= 1 & Total.Precip <=50)
sub.precip <- na.omit(sub.precip$Total.Precip)

df <- as.data.frame(sub.precip)
#write.csv(df, file = "precip_df.csv") 

# Plotting the subset distribution. The histogram also comes with the counts and densities
sub.hist <- hist(sub.precip, sub.edge) # Plotting the histogram when subsetted from 15-40.
sub.freq <- sub.hist$counts # This gives the counts of each interval (ie. frequencies)
sub.freq <- sub.freq/sum(sub.freq) # This gives the densities
sum(sub.freq) # Should sum to 1, ie. the probability density should sum up to 1.

plot(sub.x, sub.hist$density, type = "l")

sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts, sub.hist$density))

names(sub.df) <- c("precipitation", "frequency", "density")

# To-do
# 1. Fit a pdf to this data and estimate the parameters
# 2. candidate pdfs includes, Lognormal, Gamma, Inversegaussian distributions!
###############################################################################

# Fitting function to the precipitation density data.
attach(sub.df)

plot(density~precipitation)

#Define an exponential density function to fit to data
expon <- exp(-precipitation) 

fit.lm <- lm(density~expon)

fit <- fitted(fit.lm) 

#Plotting the data with fitted function with the base plot function
plot(density ~ precipitation, data= sub.df)
lines(fit, col="red")
summary(fit.lm)


# Estimating the lambda parameter with the MLE
exponential <- function(lambda){
  d <- lambda*exp(-lambda*precipitation)
  }

#expon.negLL = function(lambda){
#  -sum(dnorm(density, mean=exponential(lambda), sd=1, log=TRUE))
#} 

expon.negLL = function(lambda){
  -sum(dlnorm(density, mean=exponential(lambda), sd=1, log=TRUE))
} 

mle.expon <- mle2(expon.negLL,start=list(lambda = 1), ) 

lambda <- unname(coef(mle.expon))
lambda

plot(density~precipitation, data = sub.df) # Probability 
lines(precipitation, exponential(lambda), col="blue", pch=20)
lines(fit, col="red")

################################################################################


lognorm <- function(mu, sigma){
  d <- dx*(1/(sigma*sqrt(2*pi)))*exp(-((log(sub.x)-mu)^2)/(2*sigma^2))
} #  pdf with x being changed to sub.x. Sub.x is the x-vector for the subset data

sub.negLL = function(mu, sigma){
  -sum(dnorm(sub.freq, mean=lognorm(mu, sigma), sd=1, log=TRUE))
} # to fit to the two parameters of the model

sub.mle.lognorm <- mle2(sub.negLL, start=list(mu=5, sigma=4), ) 

a.est <- unname(coef(sub.mle.lognorm))[-4]

# fitted lognormal distrubution
print(a.est)

# Initial fitted values 
mu <- a.est[1]
sigma <- a.est[2]

# # The probability distribution is described by
plot(sub.x, lognorm(a.est[1], a.est[2]), type="l") # Probability 
#lines(sub.x, cumsum(lognorm(a.est[1], a.est[2])), col="red") # Cumlative
points(sub.x, sub.freq, col="blue", pch=20)

# Combining the distributions into one distribution
sub.sal.df <- as.data.frame(cbind(sub.x, lognorm(a.est[1], a.est[2]),
                                  cumsum(lognorm(a.est[1], a.est[2]))))
names(sub.sal.df) <- c("precipitation", "PDF", "CDF")
head(sub.sal.df)


#####################################################################################
# fitting to lognormal distribution
z <- df$sub.precip #rgamma(50, 3, 10) + rnorm(50, 0, .02)

#fit our dataset to a lognormal distribution using mle
fit <- fitdist(z, distr = "lnorm", method = "mle")

#view the summary of the fit
summary(fit)

#produce plots to visualize the fit
plot(fit)





#load package
#################
#generate 50 random values that follow a gamma distribution with shape parameter = 3
#and shape parameter = 10 combined with some gaussian noise
z <- df$sub.precip #rgamma(50, 3, 10) + rnorm(50, 0, .02)
plotdist(z, histo = TRUE, demp = TRUE)
descdist(z, boot = 1000)


#fit our dataset to a gamma distribution using mle
fit <- fitdist(z, distr = "gamma", method = "mle")

#view the summary of the fit
summary(fit)

#produce plots to visualize the fit
plot(fit)


#######################
# fitting to lognormal distribution
z <- df$sub.precip #rgamma(50, 3, 10) + rnorm(50, 0, .02)

#fit our dataset to a gamma distribution using mle
fit <- fitdist(z, distr = "lnorm", method = "mle")

#view the summary of the fit
summary(fit)

#produce plots to visualize the fit
plot(fit)

plot(z, dlnorm(z, 1.5518819, 0.9536246), col = "red")
points(z)


####################
# Another way to do this is as follows! This works great!
z <- df$sub.precip #rgamma(50, 3, 10) + rnorm(50, 0, .02)

fw <- fitdist(z, "weibull")
fg <- fitdist(z, "gamma")
fln <- fitdist(z, "lnorm")
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal", "gamma")

denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)
####################
summary(fw)
summary(fg)
summary(fln)



##############################################################################
# Find the m value (m is a parameter in the asymmetric Laplace distribution)
sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts))
sub.df[which.max(sub.df$V2),]  # The maxium value for m therefore is 32.2 
sub.m <- sub.df[which.max(sub.df$V2),][1,1] # Pulls only the m value
sigma=1

# Need to update the asymmetric laplace distribution so that it pulls the right 
# x value

sub.asy.laplace <- function(sub.m, pho, k){
  d <- dx*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
} # Same distrubution but x was changed to sub.x. Sub.x is the x-vector for the subset data

sub.negLL = function(pho, k){
  -sum(dnorm(sub.freq, mean=sub.asy.laplace(sub.m, pho, k), sd=sigma, log=TRUE))
} # Only want to fit to pho and k so only include them within the code

sub.mle.asy.laplace.Ex <- mle2(sub.negLL,
                               start=list(pho=1, k=3), ) # Similary, only want to fit to pho and k

a.est <- unname(coef(sub.mle.asy.laplace.Ex))[-4]

# fitted asymmetric distrubution
print(a.est)

# Initial fitted values 
pho0 <- a.est[1]
k0 <- a.est[2]

# # The probability distribution is described by
plot(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]), type="l", ylim=c(0,1.1)) # Probability 
#lines(sub.x, cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2])), col="red") # Cumlative
points(sub.x, sub.freq, col="blue", pch=20)

# Combining the distributions into one distribution
sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, a.est[1], a.est[2]),
                                  cumsum(sub.asy.laplace(sub.m, a.est[1], a.est[2]))))
names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
head(sub.sal.df)
