# To analyse climate data and fit functions to temperature and precipitation data
#Clear workspace
rm(list = ls())

#set working directory
setwd("/Users/jbaafi/Documents/Data/code")

# Import packages into r
library(tidyverse)
library(dplyr)
library(chron)

#Import 2011-2016 climate data into r 
climate_11 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2011_P1D.csv", header = TRUE)
climate_12 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2012_P1D.csv", header = TRUE)
climate_13 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2013_P1D.csv", header = TRUE)
climate_14 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2014_P1D.csv", header = TRUE)
climate_15 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2015_P1D.csv", header = TRUE)
climate_16 <- read.csv("/Users/jbaafi/Documents/Data/en_climate_daily_ON_6153301_2016_P1D.csv", header = TRUE)

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

# Rename some of the columns of the climate_df data to make it simple
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
  theme_light()+
  xlab("Time (Days)")+
  ylab("Mean Temperature")+
  theme(axis.text.x=element_text(angle=45, hjust=1))
 
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
#a <- mean(climate.df$Mean.Temp)
#b1 <- -8.7635
#f <- a + b1*sin(2*pi*Time/365) + b1*cos(2*pi*Time/365)
#plot(Temp~Days.Since.Origin, data = climate.df)
#lines(Time, f, col = "red")

# Fit a periodic function to total precipitation data
Precip <- climate.df$Total.Precip

#Define a periodic function to fit to data
Precip.func <- sin(2*pi*Time/365) + cos(2*pi*Time/365)

fitprecip.lm <- lm(Precip~Precip.func)

fitprecip <- fitted(fitprecip.lm)  

# find predictions for original time series
predprecip <- predict(fitprecip.lm, newdata=data.frame(Time=Time))

#Plotting the data with fitted function with the base plot function
plot(Precip ~ Days.Since.Origin, data= climate.df)
lines(fitprecip, col="red")
lines(Time, predprecip, col="blue")
summary(fitprecip.lm)



