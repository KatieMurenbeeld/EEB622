## Testing 2:
## - how to work with NetCDF in R. 
## - testing models 

setwd("~/Desktop/EEB622/Final_Project/Final")

# Libraries/packages for Bayes and Multilevel

library(rstan)
library(rstanarm)
library(tidybayes)
library(lme4)
library(ggplot2)
library(tibble)

# Libraries/packages for netcdf

library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)

# From http://geog.uoregon.edu/bartlein/courses/geog490/week04-netCDF.html#introduction
# Set a path and filename

ncpath <- "data_in/"
agbname <- "AGB"
gppname <- "GPP"
npname <- "NPLANT"

# Focus on AGB first

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

time <- ncvar_get(ncin_cagb,"year")
time

cagb_array <- ncvar_get(ncin_cagb,agbname)
dim(cagb_array)
cagb_array

cntrl_agb <- cagb_array[]

lagb_array <- ncvar_get(ncin_lagb,agbname)
dim(lagb_array)
lagb_array

log_agb <- lagb_array[]

lsagb_array <- ncvar_get(ncin_lsagb,agbname)
dim(lsagb_array)
lsagb_array

logseq_agb <- lsagb_array[]

treatments <- c("Control", "Log_Single", "Log_Sequential")

agb_df01 <- matrix(cbind(cntrl_agb, log_agb, logseq_agb), 
                   nrow = length(time), 
                   ncol = length(treatments),
                   dimnames = list(time, treatments))
dim(agb_df01)

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

AGB <- read.csv("data_out/AGB_annual.csv")
AGB2 <- read.csv("data_out/AGB_annualv02.csv")

hist(AGB$Control)
hist(AGB$Log_Single)
hist(AGB$Log_Sequential)

hist(AGB2$Control)
hist(AGB2$Log_Single)
hist(AGB2$Log_Sequential)

# Testing with a Gamma distributed model

# Control
cntrl_agb_gam <- glm(Control ~ Year, data=AGB2, family="Gamma"(link="log")) # Set link function to log.
summary(cntrl_agb_gam)
exp(coef(cntrl_agb_gam))

## Plot the fits of the gamma models:
yearsim <- seq(from=min(AGB2$Year),
              to = max(AGB2$Year), 
              length.out=nrow(AGB2))
agbpreds <- predict(cntrl_agb_gam, 
                      data.frame(Year=yearsim),
                      type="response")
predictions<- data.frame(agbpreds, yearsim)

plot(AGB2$Control~AGB2$Year, xlab = "Year of Simulation", ylab = "Above Ground Biomass (AGB)")
lines(predictions$yearsim, predictions$agbpreds, type="l", col="black", lwd=3) ## This plots the gamma model's prediction in blue

# Logging: Single Time Step
log_agb_gam <- glm(Log_Single ~ Year, data=AGB2, family="Gamma"(link="log")) # Set link function to log.
summary(log_agb_gam)
exp(coef(log_agb_gam))

log_agbpreds <- predict(log_agb_gam, 
                    data.frame(Year=yearsim),
                    type="response")
log_predictions<- data.frame(log_agbpreds, yearsim)

lines(log_predictions$yearsim, log_predictions$log_agbpreds, type="l", col="orange", lwd=3)


# Logging: Sequential Time Steps

logseq_agb_gam <- glm(Log_Sequential ~ Year, data=AGB2, family="Gamma"(link="log"))
summary(logseq_agb_gam)
exp(coef(logseq_agb_gam))

logseq_agbpreds <- predict(logseq_agb_gam, 
                        data.frame(Year=yearsim),
                        type="response")
logseq_predictions<- data.frame(logseq_agbpreds, yearsim)

lines(logseq_predictions$yearsim, logseq_predictions$logseq_agbpreds, type="l", col="blue", lty=2, lwd=3)

# Plot all together
plot(AGB2$Control~AGB2$Year, xlab = "Year of Simulation", ylab = "Above Ground Biomass (AGB)", col="black")
lines(predictions$yearsim, predictions$agbpreds, type="l", col="black", lwd=3) ## This plots the gamma model's control prediction in black
lines(log_predictions$yearsim, log_predictions$log_agbpreds, type="l", col="orange", lwd=3) ## This plots the gamma model's log prediction in orange
lines(logseq_predictions$yearsim, logseq_predictions$logseq_agbpreds, type="l", col="blue", lty=2, lwd=3) ## This plots the gamma model's log seq prediction in dashed blue
legend("bottomright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))



# New idea - scale the AGB against the control? 

#AGB2$Dev_Single <- (AGB2$Log_Single - AGB2$Control)
AGB2$Dev_Single <- abs((AGB2$Log_Single - AGB2$Control))
#AGB2$Dev_Seq <- (AGB2$Log_Sequential - AGB2$Control)
AGB2$Dev_Seq <- abs((AGB2$Log_Sequential - AGB2$Control))
AGB2$Dev_Cntrl <- (AGB2$Control - AGB2$Control)

hist(AGB2$Dev_Single)
hist(AGB2$Dev_Seq)

#plot(AGB2$Dev_Cntrl ~ AGB2$Year, type="l")

plot(AGB2$Dev_Seq[101:201] ~ AGB2$Year[101:201], pch=16, ylim=c(-10,150), col="blue",
     ylab = "Absolute Diff From Control",
     xlab = "Year")
points(AGB2$Dev_Single[101:201] ~ AGB2$Year[101:201], pch=16, col="orange")
points(AGB2$Dev_Cntrl[101:201]~AGB2$Year[101:201], type="l", col="black", lwd=2)
legend("topright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))


# Start from 1st logging year 

hist(AGB2$Dev_Single[101:201])
hist(AGB2$Dev_Seq[101:201])


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

lseq_agbpreds_dev <- predict(lm_dev_seq, 
                            data.frame(Year=yearsim_dev),
                            type="response")
lseq_predictions_dev <- data.frame(lseq_agbpreds_dev, yearsim_dev)


plot(log_predictions_dev$log_agbpreds_dev~log_predictions_dev$yearsim_dev, type="l", col="orange", lwd=3,
     main = "Gamma Model",
     xlab="Year, change to time since distrubance",
     ylab = "Predicted Absolute Difference from Control AGB",
     ylim = c(-10,100))
points(lseq_predictions_dev$lseq_agbpreds_dev~lseq_predictions_dev$yearsim_dev, type="l", col="blue", lwd=3)
points(AGB2$Dev_Single[101:201] ~ AGB2$Year[101:201], pch=16, col="orange")
points(AGB2$Dev_Seq[101:201] ~ AGB2$Year[101:201], pch=16, ylim=c(-10,150), col="blue")
points(AGB2$Dev_Cntrl[101:201]~AGB2$Year[101:201], type="l", col="black", lwd=2)
legend("topright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))

mu_l <- exp(-1.49 + 0.0007*AGB2$Year)

dgamma <- dgamma(AGB2$Dev_Single, 0.22, rate=1)

hist(dgamma)

# Test Gaussian model with "dev_" values

lm_dev_sing <- glm(AGB2$Dev_Single[101:201]~AGB2$Year[101:201])
summary(lm_dev_sing)
coef(lm_dev_sing)

yearsim_devlm <- seq(from=min(AGB2$Year[101:201]),
                   to = max(AGB2$Year[101:201]), 
                   length.out=101)

log_agbpreds_devlm <- predict(lm_dev_sing, 
                            data.frame(Year=yearsim_devlm),
                            type="response")

log_predictions_devlm <- data.frame(log_agbpreds_devlm, yearsim_devlm)


lm_dev_seq <- glm(AGB2$Dev_Seq[101:201]~AGB2$Year[101:201])
summary(lm_dev_seq)
coef(lm_dev_seq)

lseq_agbpreds_devlm <- predict(lm_dev_seq, 
                             data.frame(Year=yearsim_devlm),
                             type="response")
lseq_predictions_devlm <- data.frame(lseq_agbpreds_devlm, yearsim_devlm)


plot(log_predictions_devlm$log_agbpreds_devlm~log_predictions_devlm$yearsim_devlm, type="l", col="orange", lwd=3,
     main = "Linear Model",
     xlab="Year, change to time since distrubance",
     ylab = "Predicted Absolute Difference from Control AGB",
     ylim = c(-10,100))
points(lseq_predictions_devlm$lseq_agbpreds_devlm~lseq_predictions_devlm$yearsim_devlm, type="l", col="blue", lwd=3)
points(AGB2$Dev_Single[101:201] ~ AGB2$Year[101:201], pch=16, col="orange")
points(AGB2$Dev_Seq[101:201] ~ AGB2$Year[101:201], pch=16, ylim=c(-10,150), col="blue")
points(AGB2$Dev_Cntrl[101:201]~AGB2$Year[101:201], type="l", col="black", lwd=2)
legend("topright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))







