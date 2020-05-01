## Katie Murenbeeld
## Draft Code for EEB 622 Final Project
## 14 Apr 2020
## - how to work with NetCDF in R. 
## - testing linear and gamma distribution generalized linear models 

setwd("~/Desktop/EEB622/Final_Project/Final")

# Libraries/packages for Bayes and Multilevel (may not be needed ust yet)

library(rstan)
library(rstanarm)
library(tidybayes)
library(lme4)
library(ggplot2)
library(tibble)
library(tidyr)


# Libraries/packages for netcdf

library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)

# The following data process from http://geog.uoregon.edu/bartlein/courses/geog490/week04-netCDF.html#introduction


# Set a path and filename

ncpath <- "data_in/"
agbname <- "AGB"
gppname <- "GPP" # will get to GPP later
npname <- "NPLANT" # will get to NPLANT later

# Focus on AGB first

# Create names for the nc files
ncname_cntrlagb <- "cntrl_AGB"
ncname_logagb <- "log_AGB"
ncname_logseqagb <- "logseq_AGB"

ncfname_cagb <- paste(ncpath, ncname_cntrlagb, ".nc", sep="" )
ncfname_lagb <- paste(ncpath, ncname_logagb, ".nc", sep="" )
ncfname_lsagb <- paste(ncpath, ncname_logseqagb, ".nc", sep="" )

# Open a netCDF file

ncin_cagb <- nc_open(ncfname_cagb)
print(ncin_cagb)

ncin_lagb <- nc_open(ncfname_lagb)
print(ncin_lagb)

ncin_lsagb <- nc_open(ncfname_lsagb)
print(ncin_lsagb)

time <- ncvar_get(ncin_cagb,"year") # get the time variable from the nc file
time

# Create new arrays for control, log single, and log sequential treatments 

cagb_array <- ncvar_get(ncin_cagb,agbname) # get the agb variable from the nc file
dim(cagb_array) # Check the dimensions
cagb_array

cntrl_agb <- cagb_array[]

lagb_array <- ncvar_get(ncin_lagb,agbname) # get the agb variable from the nc file
dim(lagb_array)
lagb_array

log_agb <- lagb_array[]

lsagb_array <- ncvar_get(ncin_lsagb,agbname) # get the agb variable from the nc file
dim(lsagb_array)
lsagb_array

logseq_agb <- lsagb_array[]
treatments <- c("Control", "Log_Single", "Log_Sequential") # Create a list of names for the treatment types

# Create new dataframes from time, AGB, and treatment type

agb_df01 <- matrix(cbind(cntrl_agb, log_agb, logseq_agb), 
                   nrow = length(time), 
                   ncol = length(treatments),
                   dimnames = list(time, treatments))
dim(agb_df01) # check the dimensions

agb_df02 <- matrix(cbind(time, cntrl_agb, log_agb, logseq_agb), 
                   nrow = length(time), 
                   ncol = 4,
                   dimnames = list(time, c("Year", "Control", "Log_Single", "Log_Sequential")))
dim(agb_df02)

# Write out the new agb dataframe as a .csv file

csvpath <- "data_out/"
csvname <- "AGB_annual.csv" 
csvfile <- paste(csvpath, csvname, sep="") 
write.table(agb_df01, csvfile, row.names = TRUE, sep=",") 
# Something is off, I want the rows to be the actual years
# Columns just the three treatments

csvname2 <- "AGB_annualv02.csv"
csvfile2 <- paste(csvpath, csvname2, sep="")
write.table(agb_df02, csvfile2, row.names = FALSE, sep=",") 

# Read in the new csv files

AGB2 <- read.csv("data_out/AGB_annualv02.csv")

# Check the distribution of the data

hist(AGB2$Control)
hist(AGB2$Log_Single)
hist(AGB2$Log_Sequential)

# Testing with a linear model, normal distribution

# Create a model for the Control treatment
cntrl_agb_lm <- lm(Control ~ Year, data=AGB2) 
summary(cntrl_agb_lm)
coef(cntrl_agb_lm)

## Plot the fits of the linear models:
yearsim <- seq(from=min(AGB2$Year),
               to = max(AGB2$Year), 
               length.out=nrow(AGB2))
agbpreds <- predict(cntrl_agb_lm, 
                    data.frame(Year=yearsim),
                    type="response")
predictions<- data.frame(agbpreds, yearsim)

plot(AGB2$Control~AGB2$Year, xlab = "Year of Simulation", ylab = "Above Ground Biomass (AGB)")
lines(predictions$yearsim, predictions$agbpreds, type="l", col="black", lwd=3) ## This plots the gamma model's prediction in blue

# Create a Model for Logging: Single Time Step
log_agb_lm <- lm(Log_Single ~ Year, data=AGB2) 
summary(log_agb_lm)
coef(log_agb_lm)

log_agbpreds <- predict(log_agb_lm, 
                        data.frame(Year=yearsim),
                        type="response")
log_predictions<- data.frame(log_agbpreds, yearsim)

lines(log_predictions$yearsim, log_predictions$log_agbpreds, type="l", col="orange", lwd=3)


# Create a Model for Logging: Sequential Time Steps

logseq_agb_lm <- lm(Log_Sequential ~ Year, data=AGB2)
summary(logseq_agb_lm)
coef(logseq_agb_lm)

logseq_agbpreds <- predict(logseq_agb_lm, 
                           data.frame(Year=yearsim),
                           type="response")
logseq_predictions<- data.frame(logseq_agbpreds, yearsim)

lines(logseq_predictions$yearsim, logseq_predictions$logseq_agbpreds, type="l", col="blue", lty=2, lwd=3)

# Plot all together
plot(AGB2$Control~AGB2$Year,
     main = "Linear Model, Normal Distribution",
     xlab = "Year of Simulation", 
     ylab = "Above Ground Biomass (AGB)", 
     col="black")
lines(predictions$yearsim, predictions$agbpreds, type="l", col="black", lwd=3) ## This plots the gamma model's control prediction in black
lines(log_predictions$yearsim, log_predictions$log_agbpreds, type="l", col="orange", lwd=3) ## This plots the gamma model's log prediction in orange
lines(logseq_predictions$yearsim, logseq_predictions$logseq_agbpreds, type="l", col="blue", lty=2, lwd=3) ## This plots the gamma model's log seq prediction in dashed blue
legend("bottomright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))

# Testing with a generalized linear model, Gamma distribution
# For the Gamma distribution I will do all model years as well as from the logging year ( to get AGB since disturbance)

# Control
cntrl_agb_gam <- glm(Control ~ Year, data=AGB2, family ="Gamma") 
summary(cntrl_agb_gam)
exp(coef(cntrl_agb_gam))

## Plot the fits of the gamma models:
yearsim_gam <- seq(from=min(AGB2$Year),
               to = max(AGB2$Year), 
               length.out=nrow(AGB2))
agbpreds_gam <- predict(cntrl_agb_gam, 
                    data.frame(Year=yearsim_gam),
                    type="response")
predictions_gam<- data.frame(agbpreds_gam, yearsim_gam)

plot(AGB2$Control~AGB2$Year, xlab = "Year of Simulation", ylab = "Above Ground Biomass (AGB)")
lines(predictions_gam$yearsim_gam, predictions_gam$agbpreds_gam, type="l", col="black", lwd=3) ## This plots the gamma model's prediction in blue

# Logging: Single Time Step
log_agb_gam <- glm(Log_Single ~ Year, data=AGB2, family = "Gamma") 
summary(log_agb_gam)
exp(coef(log_agb_gam))

log_agbpreds_gam <- predict(log_agb_gam, 
                        data.frame(Year=yearsim_gam),
                        type="response")
log_predictions_gam<- data.frame(log_agbpreds_gam, yearsim_gam)

lines(log_predictions_gam$yearsim_gam, log_predictions_gam$log_agbpreds_gam, type="l", col="orange", lwd=3)


# Logging: Sequential Time Steps

logseq_agb_gam <- glm(Log_Sequential ~ Year, data=AGB2, family = "Gamma")
summary(logseq_agb_gam)
exp(coef(logseq_agb_gam))

logseq_agbpreds_gam <- predict(logseq_agb_gam, 
                           data.frame(Year=yearsim_gam),
                           type="response")
logseq_predictions_gam<- data.frame(logseq_agbpreds_gam, yearsim_gam)

lines(logseq_predictions_gam$yearsim_gam, logseq_predictions_gam$logseq_agbpreds_gam, type="l", col="blue", lty=2, lwd=3)

# Plot all together
plot(AGB2$Control~AGB2$Year, 
     main = "GLM: Gamma Distribution",
     xlab = "Year of Simulation", 
     ylab = "Above Ground Biomass (AGB)", 
     col="black")
lines(predictions_gam$yearsim_gam, predictions_gam$agbpreds_gam, type="l", col="black", lwd=3) ## This plots the gamma model's control prediction in black
lines(log_predictions_gam$yearsim_gam, log_predictions_gam$log_agbpreds_gam, type="l", col="orange", lwd=3) ## This plots the gamma model's log prediction in orange
lines(logseq_predictions_gam$yearsim_gam, logseq_predictions_gam$logseq_agbpreds_gam, type="l", col="blue", lty=2, lwd=3) ## This plots the gamma model's log seq prediction in dashed blue
legend("bottomright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))

# New idea - scale the AGB against the control? Specifically for Gamma distribution and looking at time since disturbance.

AGB2$Dev_Single <- abs((AGB2$Log_Single - AGB2$Control)) # absolute because values must be positive with a Gamma distribution
AGB2$Dev_Seq <- abs((AGB2$Log_Sequential - AGB2$Control))
AGB2$Dev_Cntrl <- (AGB2$Control - AGB2$Control)

hist(AGB2$Dev_Single) # From the histograms, this looks like a good candidate for a Gamma distribution
hist(AGB2$Dev_Seq)

# Start from 1st logging year 
hist(AGB2$Dev_Single[101:201])
hist(AGB2$Dev_Seq[101:201])

plot(AGB2$Dev_Seq[101:201] ~ AGB2$Year[101:201], pch=16, ylim=c(-10,150), col="blue",
     ylab = "Absolute Diff From Control",
     xlab = "Year")
points(AGB2$Dev_Single[101:201] ~ AGB2$Year[101:201], pch=16, col="orange")
points(AGB2$Dev_Cntrl[101:201]~AGB2$Year[101:201], type="l", col="black", lwd=2)
legend("topright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))

# Test gamma model with "dev_" values

gam_dev_sing <- glm(AGB2$Dev_Single[101:201]~AGB2$Year[101:201], family="Gamma")
summary(gam_dev_sing)
exp(coef(gam_dev_sing))

yearsim_dev <- seq(from=min(AGB2$Year[101:201]),
                   to = max(AGB2$Year[101:201]), 
                   length.out=101)

log_agbpreds_dev <- predict(gam_dev_sing, 
                            data.frame(Year=yearsim_dev),
                            type="response")

log_predictions_dev <- data.frame(log_agbpreds_dev, yearsim_dev)


gam_dev_seq <- glm(AGB2$Dev_Seq[101:201]~AGB2$Year[101:201], family="Gamma")
summary(gam_dev_seq)
exp(coef(gam_dev_seq))

lseq_agbpreds_dev <- predict(gam_dev_seq, 
                             data.frame(Year=yearsim_dev),
                             type="response")
lseq_predictions_dev <- data.frame(lseq_agbpreds_dev, yearsim_dev)


plot(log_predictions_dev$log_agbpreds_dev~log_predictions_dev$yearsim_dev, type="l", col="orange", lwd=3,
     main = "Gamma Model",
     xlab="Year, change to time since distrubance",
     ylab = "Absolute Difference from Control AGB",
     ylim = c(-10,100))
points(lseq_predictions_dev$lseq_agbpreds_dev~lseq_predictions_dev$yearsim_dev, type="l", col="blue", lwd=3)
points(AGB2$Dev_Single[101:201] ~ AGB2$Year[101:201], pch=16, col="orange")
points(AGB2$Dev_Seq[101:201] ~ AGB2$Year[101:201], pch=16, ylim=c(-10,150), col="blue")
points(AGB2$Dev_Cntrl[101:201]~AGB2$Year[101:201], type="l", col="black", lwd=2)
legend("topright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))

# Look at a metric of model fit, R2

# Create a function to calculate R2
r2 <- function(y_hat,y){ 
  RSS<-sum((((y_hat))-(y))^2) 
  TSS<-sum(((y)-(mean(y)))^2) 
  return(1-RSS/TSS)}

# Find the R2 value for the linear and gamma models over the entire dataset
lm_full_cntrl <- r2(AGB2$Control, predictions$agbpreds)
print(lm_full_cntrl)
lm_full_log <- r2(AGB2$Log_Single, log_predictions$log_agbpreds)
print(lm_full_log)
lm_full_logseq <- r2(AGB2$Log_Sequential, logseq_predictions$logseq_agbpreds)
print(lm_full_logseq)
# All r2 values are negative, models are no better than a horizontal line

gam_full_cntrl <- r2(AGB2$Control, predictions_gam$agbpreds_gam)
print(gam_full_cntrl)
gam_full_log <- r2(AGB2$Log_Single, log_predictions_gam$log_agbpreds_gam)
print(gam_full_log)
gam_full_logseq <- r2(AGB2$Log_Sequential, logseq_predictions_gam$logseq_agbpreds_gam)
print(gam_full_logseq)
# All r2 values are negative, models are no better than a horizontal line

# Find the R2 value for the gamma model over the data subset for time since logging

gam_dev_log <- r2(AGB2$Dev_Single[101:201], log_predictions_dev$log_agbpreds_dev)
print(gam_dev_log)
# A negative r2, the model is no better than a horizontal line
gam_dev_logseq <- r2(AGB2$Dev_Sequential[101:201], lseq_predictions_dev$lseq_agbpreds_dev)
print(gam_dev_logseq)
# r2 is 1! Model fits the data perfectly...which I do not trust.


