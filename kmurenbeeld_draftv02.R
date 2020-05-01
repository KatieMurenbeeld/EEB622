## Katie Murenbeeld
## Draft Code for EEB 622 Final Project
## 24 Apr 2020
## - how to work with NetCDF in R. 
## - testing linear and gamma distribution generalized linear models 
## - model fit and error

setwd("~/Desktop/EEB622/Final_Project/Final")

# Libraries/packages for Bayes and Multilevel (may not be needed ust yet)

library(tidyverse)
library(rstan)
library(rstanarm)
library(tidybayes)
library(lme4)
library(ggplot2)
library(tibble)
library(tidyr)
library(brms)
library(loo)


# Libraries/packages for netcdf

library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)

# The following data process from http://geog.uoregon.edu/bartlein/courses/geog490/week04-netCDF.html#introduction

# Set a path and filename

ncpath <- "data_in/"
agbname <- "AGB"

# Focus on AGB first

# Create names for the nc files
cntrl <- "cntrl_AGB"
log <- "log_AGB"
logseq <- "logseq_AGB"
cntrl.p <- "cntrl_AGB_polly_v03"
log.p <- "log_AGB_polly_v03"
logseq.p <- "logseq_AGB_polly_v03"
cntrl.pt3 <- "cntrl_AGB_pt3_v01"
log.pt3 <- "log_AGB_pt3_v01"
logseq.pt3 <- "logseq_AGB_pt3_v01"

ncfname_cagb <- paste(ncpath, cntrl, ".nc", sep="" )
ncfname_lagb <- paste(ncpath, log, ".nc", sep="" )
ncfname_lsagb <- paste(ncpath, logseq, ".nc", sep="" )
ncfname_cagb.p <- paste(ncpath, cntrl.p, ".nc", sep="" )
ncfname_lagb.p <- paste(ncpath, log.p, ".nc", sep="" )
ncfname_lsagb.p <- paste(ncpath, logseq.p, ".nc", sep="" )
ncfname_cagb.pt3 <- paste(ncpath, cntrl.pt3, ".nc", sep="" )
ncfname_lagb.pt3 <- paste(ncpath, log.pt3, ".nc", sep="" )
ncfname_lsagb.pt3 <- paste(ncpath, logseq.pt3, ".nc", sep="" )

# Open a netCDF file

ncin_cagb <- nc_open(ncfname_cagb)
print(ncin_cagb)

ncin_lagb <- nc_open(ncfname_lagb)
print(ncin_lagb)

ncin_lsagb <- nc_open(ncfname_lsagb)
print(ncin_lsagb)

ncin_cagb.p <- nc_open(ncfname_cagb.p)
print(ncin_cagb.p)

ncin_lagb.p <- nc_open(ncfname_lagb.p)
print(ncin_lagb.p)

ncin_lsagb.p <- nc_open(ncfname_lsagb.p)
print(ncin_lsagb.p)

ncin_cagb.pt3 <- nc_open(ncfname_cagb.pt3)
print(ncin_cagb.pt3)

ncin_lagb.pt3 <- nc_open(ncfname_lagb.pt3)
print(ncin_lagb.pt3)

ncin_lsagb.pt3 <- nc_open(ncfname_lsagb.pt3)
print(ncin_lsagb.pt3)

time <- ncvar_get(ncin_cagb,"year") # get the time variable from the nc file
time

time.p <- ncvar_get(ncin_cagb.p,"year") # get the time variable from the nc file
time.p

time.pt3 <- ncvar_get(ncin_cagb.pt3,"year") # get the time variable from the nc file
time.pt3

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

cagb_array.p <- ncvar_get(ncin_cagb.p,agbname) # get the agb variable from the nc file
dim(cagb_array.p) # Check the dimensions
cagb_array.p

cntrl_agb.p <- cagb_array.p[]

lagb_array.p <- ncvar_get(ncin_lagb.p,agbname) # get the agb variable from the nc file
dim(lagb_array.p)
lagb_array.p

log_agb.p <- lagb_array.p[]

lsagb_array.p <- ncvar_get(ncin_lsagb.p,agbname) # get the agb variable from the nc file
dim(lsagb_array.p)
lsagb_array.p

logseq_agb.p <- lsagb_array.p[]

cagb_array.pt3 <- ncvar_get(ncin_cagb.pt3,agbname) # get the agb variable from the nc file
dim(cagb_array.pt3) # Check the dimensions
cagb_array.pt3

cntrl_agb.pt3 <- cagb_array.pt3[]

lagb_array.pt3 <- ncvar_get(ncin_lagb.pt3,agbname) # get the agb variable from the nc file
dim(lagb_array.pt3)
lagb_array.pt3

log_agb.pt3 <- lagb_array.pt3[]

lsagb_array.pt3 <- ncvar_get(ncin_lsagb.pt3,agbname) # get the agb variable from the nc file
dim(lsagb_array.pt3)
lsagb_array.pt3

logseq_agb.pt3 <- lsagb_array.pt3[]

treatments <- c("Control", "Log_Single", "Log_Sequential") # Create a list of names for the treatment types

# Create new dataframes from time, AGB, and treatment type

agb_df <- matrix(cbind(time, cntrl_agb, log_agb, logseq_agb), 
                   nrow = length(time), 
                   ncol = 4,
                   dimnames = list(time, c("Year", "Control", "Log_Single", "Log_Sequential")))
dim(agb_df)

agb_df.p <- matrix(cbind(time.p, cntrl_agb.p, log_agb.p, logseq_agb.p), 
                 nrow = length(time.p), 
                 ncol = 4,
                 dimnames = list(time.p, c("Year", "Control", "Log_Single", "Log_Sequential")))
dim(agb_df.p)

agb_df.pt3 <- matrix(cbind(time.pt3, cntrl_agb.pt3, log_agb.pt3, logseq_agb.pt3), 
                   nrow = length(time.pt3), 
                   ncol = 4,
                   dimnames = list(time.pt3, c("Year", "Control", "Log_Single", "Log_Sequential")))
dim(agb_df.pt3)

# Write out the new agb dataframe as a .csv file

csvpath <- "data_out/"

csvname1 <- "AGB_df.csv" 
csvname2 <- "AGB_dfp.csv"
csvname3 <- "AGB_dfpt3.csv"

csvfile1 <- paste(csvpath, csvname1, sep="")
write.table(agb_df, csvfile1, row.names = FALSE, sep=",") 

csvfile2 <- paste(csvpath, csvname2, sep="")
write.table(agb_df.p, csvfile2, row.names = FALSE, sep=",") 

csvfile3 <- paste(csvpath, csvname3, sep="")
write.table(agb_df.pt3, csvfile3, row.names = FALSE, sep=",") 

# Read in the new csv files

AGB1 <- read_csv("data_out/AGB_df.csv")
AGB2 <- read_csv("data_out/AGB_dfp.csv")
AGB3 <- read_csv("data_out/AGB_dfpt3.csv")

AGB1$Dev_Cntrl <- abs(AGB1$Control - AGB1$Control)
AGB1$Dev_Single <- abs((AGB1$Log_Single - AGB1$Control)) # absolute because values must be positive with a Gamma distribution
AGB1$Dev_Seq <- abs((AGB1$Log_Sequential - AGB1$Control))

AGB1$AGBcntrl.c <- AGB1$Control - mean(AGB1$Control)
AGB1$AGBlog.c <- AGB1$Log_Single - mean(AGB1$Log_Single)
AGB1$AGBlogseq.c <- AGB1$Log_Sequential - mean(AGB1$Log_Sequential)

AGB2$Dev_Cntrl <- abs(AGB2$Control - AGB2$Control)
AGB2$Dev_Single <- abs((AGB2$Log_Single - AGB2$Control)) # absolute because values must be positive with a Gamma distribution
AGB2$Dev_Seq <- abs((AGB2$Log_Sequential - AGB2$Control))

AGB2$AGBcntrl.c <- AGB2$Control - mean(AGB2$Control)
AGB2$AGBlog.c <- AGB2$Log_Single - mean(AGB2$Log_Single)
AGB2$AGBlogseq.c <- AGB2$Log_Sequential - mean(AGB2$Log_Sequential)

AGB3$Dev_Cntrl <- abs(AGB3$Control - AGB3$Control)
AGB3$Dev_Single <- abs((AGB3$Log_Single - AGB3$Control)) # absolute because values must be positive with a Gamma distribution
AGB3$Dev_Seq <- abs((AGB3$Log_Sequential - AGB3$Control))

AGB3$AGBcntrl.c <- AGB3$Control - mean(AGB3$Control)
AGB3$AGBlog.c <- AGB3$Log_Single - mean(AGB3$Log_Single)
AGB3$AGBlogseq.c <- AGB3$Log_Sequential - mean(AGB3$Log_Sequential)

head(AGB1)
head(AGB2)
head(AGB3)

tidyAGB1 <- AGB1 %>%
  gather(Treatment, AGB, Control:Log_Sequential)
  #gather(Dev_Treatment, AGB.d, Dev_Cntrl:Dev_Seq)
  #gather(Cen_Treatment, AGB.c, AGBcntrl.c:AGBlogseq.c)
  #select(-Dev_Treatment, -Cen_Treatment)

tidyAGB1 <- tidyAGB1 %>%
  select(-Dev_Cntrl, -Dev_Single, -Dev_Seq, -AGBcntrl.c, -AGBlog.c, -AGBlogseq.c)

tidyAGB1.d <- AGB1 %>%
  gather(Dev_Treatment, AGB.d, Dev_Cntrl:Dev_Seq)

tidyAGB1.d <- tidyAGB1.d %>%
  select(-Control, -Log_Single, -Log_Sequential, -AGBcntrl.c, -AGBlog.c, -AGBlogseq.c)

tidyAGB1.c <- AGB1 %>%
  gather(Cen_Treatment, AGB.c, AGBcntrl.c:AGBlogseq.c)

tidyAGB1.c <- tidyAGB1.c %>%
  select(-Control, -Log_Single, -Log_Sequential, -Dev_Cntrl, -Dev_Single, -Dev_Seq)


tidyAGB1tot <- cbind(tidyAGB1, "AGB.c" = tidyAGB1.c$AGB.c, "AGB.d" = tidyAGB1.d$AGB.d)

head(tidyAGB1tot)

tidyAGB1tot$rep <- 1 # add a column for replicate (so this is test 1)

head(tidyAGB1tot)

tidyAGB2 <- AGB2 %>% 
  gather(Treatment, AGB, Control:Log_Sequential)

tidyAGB2 <- tidyAGB2 %>%
  select(-Dev_Cntrl, -Dev_Single, -Dev_Seq, -AGBcntrl.c, -AGBlog.c, -AGBlogseq.c)

tidyAGB2.d <- AGB2 %>%
  gather(Dev_Treatment, AGB.d, Dev_Cntrl:Dev_Seq)

tidyAGB2.d <- tidyAGB2.d %>%
  select(-Control, -Log_Single, -Log_Sequential, -AGBcntrl.c, -AGBlog.c, -AGBlogseq.c)

tidyAGB2.c <- AGB2 %>%
  gather(Cen_Treatment, AGB.c, AGBcntrl.c:AGBlogseq.c)

tidyAGB2.c <- tidyAGB2.c %>%
  select(-Control, -Log_Single, -Log_Sequential, -Dev_Cntrl, -Dev_Single, -Dev_Seq)


tidyAGB2tot <- cbind(tidyAGB2, "AGB.c" = tidyAGB2.c$AGB.c, "AGB.d" = tidyAGB2.d$AGB.d)

tidyAGB2tot$rep <- 2

head(tidyAGB2tot)

tidyAGB3 <- AGB3 %>% 
  gather(Treatment, AGB, Control:Log_Sequential)

tidyAGB3 <- tidyAGB3 %>%
  select(-Dev_Cntrl, -Dev_Single, -Dev_Seq, -AGBcntrl.c, -AGBlog.c, -AGBlogseq.c)

tidyAGB3.d <- AGB3 %>%
  gather(Dev_Treatment, AGB.d, Dev_Cntrl:Dev_Seq)

tidyAGB3.d <- tidyAGB3.d %>%
  select(-Control, -Log_Single, -Log_Sequential, -AGBcntrl.c, -AGBlog.c, -AGBlogseq.c)

tidyAGB3.c <- AGB3 %>%
  gather(Cen_Treatment, AGB.c, AGBcntrl.c:AGBlogseq.c)

tidyAGB3.c <- tidyAGB3.c %>%
  select(-Control, -Log_Single, -Log_Sequential, -Dev_Cntrl, -Dev_Single, -Dev_Seq)


tidyAGB3tot <- cbind(tidyAGB3, "AGB.c" = tidyAGB3.c$AGB.c, "AGB.d" = tidyAGB3.d$AGB.d)

tidyAGB3tot$rep <- 3

head(tidyAGB3tot)

AGBtot <- rbind(tidyAGB1tot,tidyAGB2tot,tidyAGB3tot)

AGBtot <- mutate(AGBtot,  
                   AGB.sc = scale(AGBtot$AGB))
 
AGBtot$rep <- as.factor(AGBtot$rep)
head(AGBtot)

# Write out AGBtot as a new .csv file

csvname4 <- "AGBtot.csv"
csvfile4 <- paste(csvpath, csvname4, sep="")
write.table(AGBtot, csvfile4, row.names = FALSE, sep=",") 

AGBpost <- filter(AGBtot, Year > 2079) 
AGBpost2 <- filter(AGBpost, Treatment != "Control")

AGBpost$Replicate <- as.factor(AGBpost$rep)

ggplot(data = AGBpost, aes(x=Year, y=AGB.c, color=Treatment)) + 
  labs(x="Year", y = "Above Ground Biomass - centered (gC/m^2)") + 
  geom_line() +
  theme_classic()

ggplot(data = AGBpost, aes(x=Year, y=AGB.sc, color=Treatment)) + 
  labs(x="Year", y = "Above Ground Biomass - scaled (gC/m^2)") + 
  geom_point() +
  theme_classic()


ggplot(data = AGBpost, aes(x = Treatment, y=AGB.c)) + 
  geom_boxplot()  +
  labs(x="Treatment", y = "Above Ground Biomass (centered) [$gC/m^2$]")

ggplot(data = AGBpost, aes(x = Replicate, y=AGB.c)) + 
  geom_boxplot()  +
  labs(x="Replicate", y = "Above Ground Biomass (centered) [$gC/m^2$]")

ggplot(data = AGBpost, aes(x=Year, y=AGB.d, color=Treatment)) + 
  labs(x="Year", y = "Above Ground Biomass - abs. difference from Control (gC/m^2)") + 
  geom_line() +
  theme_classic()

# Testing out Gamma GLM

modgam_dev <- glm(AGB.d~Year, data=AGBpost2, family="Gamma"(link="log"))
summary(modgam_dev)
exp(coef(modgam_dev))

mulimodgam_dev <- glmer(AGB.d~Year + (1|Treatment/rep), data=AGBpost2, family="Gamma"(link="log"))
summary(multimodgam_dev)
exp(coef(multimodgam_dev))
ranef(multimodgam_dev)

# Breakout the treatments, run glm Gamma and compare the models with AIC or something else

single <- filter(AGBpost, Treatment=="Log_Single") %>%
  select(-AGB.c,-AGB.sc, -rep)

sequential <- filter (AGBpost, Treatment=="Log_Sequential") %>%
  select(-AGB.c,-AGB.sc, -rep)

ggplot(data = single, aes(x=Year, y=AGB, color=Replicate)) + 
  labs(x="Year", y = "Above Ground Biomass (gC/m^2)") + 
  geom_point() +
  theme_classic()

ggplot(data = sequential, aes(x=Year, y=AGB, color=Replicate)) + 
  labs(x="Year", y = "Above Ground Biomass (gC/m^2)") + 
  geom_point() +
  theme_classic()

gam_single <- glm(AGB~Year, data=single, family="Gamma"(link="log"))
summary(gam_single)
exp(coef(gam_single))

yearsim_sing <- seq(from=min(single$Year),
                   to = max(single$Year), 
                   length.out=101)

sing_agbpreds <- predict(gam_single, 
                            data.frame(Year=yearsim_sing),
                            type="response")

sing_predictions <- data.frame(sing_agbpreds, yearsim_sing)




gam_seq <- glm(AGB~Year, data=sequential, family="Gamma"(link="log"))
summary(gam_seq)
exp(coef(gam_seq))

yearsim_seq <- seq(from=min(single$Year),
                    to = max(single$Year), 
                    length.out=101)

seq_agbpreds <- predict(gam_seq, 
                         data.frame(Year=yearsim_seq),
                         type="response")

seq_predictions_dev <- data.frame(seq_agbpreds, yearsim_seq)


plot(sing_predictions$sing_agbpreds~sing_predictions$yearsim_sing, type="l", col="orange", lwd=3,
     main = "Gamma Model",
     xlab="Year, change to time since distrubance",
     ylab = "AGB",
     ylim = c(-10,1050))
points(seq_predictions_dev$seq_agbpreds~seq_predictions_dev$yearsim_seq, type="l", col="blue", lwd=3)
points(single$AGB ~ single$Year, pch=16, col="orange")
points(sequential$AGB ~ sequential$Year, pch=16, ylim=c(-10,150), col="blue")
legend("topright", pch=16, col = c("black", "orange","blue"), legend = c("Control", "Single","Sequential"))


gam_sing_dev <- glm(AGB.d~Year, data=single, family="Gamma"(link="log"))
summary(gam_sing_dev)
exp(coef(gam_sing_dev))

yearsim_singdev <- seq(from=min(single$Year),
                    to = max(single$Year), 
                    length.out=101)

singdev_agbpreds <- predict(gam_sing_dev, 
                         data.frame(Year=yearsim_singdev),
                         type="response")

singdev_predictions <- data.frame(singdev_agbpreds, yearsim_singdev)

gam_seq_dev <- glm(AGB.d~Year, data=sequential, family="Gamma"(link="log"))
summary(gam_seq_dev)
exp(coef(gam_seq_dev))

yearsim_seqdev <- seq(from=min(sequential$Year),
                       to = max(sequential$Year), 
                       length.out=101)

seqdev_agbpreds <- predict(gam_seq_dev, 
                            data.frame(Year=yearsim_singdev),
                            type="response")

seqdev_predictions <- data.frame(singdev_agbpreds, yearsim_singdev)

plot(singdev_predictions$singdev_agbpreds~singdev_predictions$yearsim_singdev, type="l", col="orange", lwd=3,
     main = "Gamma Model",
     xlab="Year, change to time since distrubance",
     ylab = "AGB",
     ylim = c(-10,2000))
points(seqdev_predictions$seqdev_agbpreds~seqdev_predictions$yearsim_seqdev, type="l", col="blue", lwd=3)
points(single$AGB.d ~ single$Year, pch=16, col="orange")
points(sequential$AGB.d ~ single$Year, pch=16, col="blue")


multigam_seq_dev <- glmer(AGB.d~Year + (1|Replicate), data=sequential, family="Gamma")
summary(gam_sing_dev)
exp(coef(gam_sing_dev))
