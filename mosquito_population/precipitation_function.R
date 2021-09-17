# Created by: Joseph Baafi
# Date: June 18, 2021
# Purpose: To fit a rainfall/precipitation function to data

#Clear workspace
rm(list = ls())

#set working directory
setwd("/Users/jbaafi/Desktop/Mosquito population dynamics/mosquito_population")

# Import packages into r
library(tidyverse)
library(scales)
library(dplyr)
library(chron)
library(ggplot2)


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

#write.csv(climate.df, file = "climate.df.csv") 

# Time series plot of total precipitation
ggplot(data = climate.df, mapping = aes(x=Chron.Date, y=Total.Precip))+
  geom_line(colour = "steelblue")+
  theme_light()+
  xlab("Time (Days)")+
  ylab("Total Precipitation") +
  theme(axis.text.x=element_text(angle=45, hjust=1))


# Fit a periodic function to total precipitation data
Precip <- climate.df$Total.Precip
Time <- climate.df$Days.Since.Origin

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
summary(predprecip)

#Plotting with same data as above with ggplot2
ggplot(data = climate.df)+
  geom_point(mapping = aes(x=Days.Since.Origin, y=Precip))+
  geom_line(mapping = aes(x=Days.Since.Origin, y=fitprecip, colour="red"))+
  geom_line(mapping = aes(x=Days.Since.Origin, y=predprecip, colour="blue"))

#################################################################################
#function parameters
p <-   2.16164 #mean(climate.df$Total.Precip) #
b <-  -0.05589
c <- b

t_end <- 365

t <-  seq(0, t_end, #dtime
)
t0 <-t[1]

# Defining temperature function
precip <- function(t){
  precip = p + b*sin(2 * pi * t/365) + c*cos(2 * pi * t/365)
  return(precip)
}

precip.df <- data.frame(t, precip(t))

#Plot of temperature as a function of time (365 days)
ggplot(data = precip.df, mapping = aes(x=t, y=precip(t)))+
  geom_line()

#Want to plot only first year data for precipitation 
year1 <- head(climate.df, 375)

# base plot
plot(year1$Days.Since.Origin, year1$Total.Precip)

#ggplot
ggplot(data = year1, mapping = aes(x=Days.Since.Origin, y=Total.Precip))+
  geom_point()



