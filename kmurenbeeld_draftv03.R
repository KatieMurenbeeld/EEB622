## Effect of Logging Activity Sequencing on Above Ground Biomass Through Time
## Author: Katie Murenbeeld
## Date: 1 May 2020

# Load packages/libraries
library(tidyverse)
library(rstan)
library(rstanarm)
library(tidybayes)
library("bayesplot")
library(lme4)
library(ggplot2)
library(tibble)
library(tidyr)
library(brms)
library(loo)


# Read in the csv file

AGBtot <- read_csv("data_out/AGBtot.csv")
AGBtot$rep <- as.factor(AGBtot$rep)

## DATA: This data was collected from model simulations using the Community Land Model Functionally Assembled Terrestrial Ecosystem Simulator (CLM(FATES))
## In this study logging treatment activities were performed at different timelines to determine the impact of logging representation on simulated Above Ground Biomass (AGB).
## This dataset was originally in a NetCDF format. Data was reorganized and "tidied" in another set of R scripts. If you are interested in how the author 
## reorganized and "tidied" the data, then please visit the repository on GitHub: <insert link here>.
## YEAR: Is the model year in the simulation.
## Treatment: Logging treatments tested in the CLM(FATES) simulations. Control = no logging occurred, Log_Single = logging occurred at a single time step (model year 2079),
## Log_Sequential = logging that occured at multiple time steps (model years 2079, 2081, 2083).
## AGB: Average annual simulated Above Ground Biomass (gC/m^2)
## AGB.c: Average annual simulated Above Ground Biomass (AGB) centered around each simulation's mean AGB
## AGB.d: The absolute difference in Above Ground Biomass from the control AGB for each replicate
## rep: The replicate of the simulation. In this case there are three replicates with differing parameter or climate settings. 
## Within each replicate was completed a simulation for each of the three Treatment types.
## AGB.sc: Scaled Above Ground Biomass using default base R scaling function.


# Create data table for only the post treatment years
AGBpost <- filter(AGBtot, Year > 2079) 


# Create data table for 20, 50, and 100 years post treatment, excluding the "Control" treatment
AGBpost2 <- filter(AGBpost, Treatment != "Control")
AGBpost20 <- filter(AGBpost2, Year < 2100)
AGBpost50 <- filter(AGBpost2, Year < 2130)

# Create data tables for the "Single" and "Sequential" treatments
single <- filter(AGBpost, Treatment=="Log_Single") %>%
  select(-AGB.c,-AGB.sc)


sequential <- filter (AGBpost, Treatment=="Log_Sequential") %>%
  select(-AGB.c,-AGB.sc)


# Plot boxplots for the treatments and the replicates

ggplot(data = AGBtot, aes(x = Treatment, y=AGB)) + 
  geom_boxplot()  +
  labs(x="Treatment", y = "Above Ground Biomass (gC/m^2)", title = "Above Ground Biomass: Boxplot by Treatment")

ggplot(data = AGBtot, aes(x = rep, y=AGB)) + 
  geom_boxplot()  +
  labs(x="Replicate", y = "Above Ground Biomass (gC/m^2)", title = "Above, Ground Biomass: Boxplot by Replicate")

# Plot time series of the data

ggplot(data = AGBpost, aes(x=Year, y=AGB, color=Treatment)) + 
  labs(x="Year", y = "Above Ground Biomass (gC/m^2)", title = "Above Ground Biomass Post Disturbance") + 
  geom_line() +
  theme_classic()

ggplot(data = AGBpost, aes(x=Year, y=AGB.c, color=Treatment)) + 
  labs(x="Year", y = "Above Ground Biomass - Centered (gC/m^2)") + 
  geom_line() +
  theme_classic()

ggplot(data = AGBpost, aes(x=Year, y=AGB.d, color=Treatment)) + 
  labs(x="Year", y = "Above Ground Biomass - abs. difference from Control (gC/m^2)", title = " Above Ground Biomass - Diff from Control") + 
  geom_line() +
  theme_classic()

# No Pooling
ggplot(AGBpost) + 
  aes(x = Year, y = AGB) + 
  stat_smooth(method = "glm", se = FALSE) + # Plot a trend line for each truck driver
  geom_point() + # Put the points on top of lines
  facet_wrap("Treatment") + # Facet lets us plot each subject's responses separately
  labs(x = "Days of sleep deprivation", y ="Average reaction time (ms)") + 
  scale_x_continuous(breaks = 0:4 * 2)

ggplot(AGBpost) + 
  aes(x = Year, y = AGB.d) + 
  stat_smooth(method = "lm", se = FALSE) + # Plot a trend line for each truck driver
  geom_point() + # Put the points on top of lines
  facet_wrap("rep") + # Facet lets us plot each subject's responses separately
  labs(x = "Years Post Treatment", y ="Above Ground Biomass - abs. difference from Control (gC/m^2)") + 
  scale_x_continuous(breaks = 0:4 * 2)

# Testing out Gamma GLM
modgam_dev <- glm(AGB.d~Year, data=AGBpost2, family="Gamma"(link="log"))
summary(modgam_dev)
exp(coef(modgam_dev))


# Testing a Gamma GLM for each treatment: Log_Single

gam_single <- glm(AGB.d~Year, data=single, family="Gamma"(link="log"))
summary(gam_single)
exp(coef(gam_single))

yearsim_sing <- seq(from=min(single$Year),  
                    to = max(single$Year), 
                    length.out=101)

sing_agbpreds <- predict(gam_single,    
                         data.frame(Year=yearsim_sing),
                         type="response")

sing_predictions <- data.frame(sing_agbpreds, yearsim_sing)


# Repeat for Log_sequential

gam_seq <- glm(AGB.d~Year, data=sequential, family="Gamma"(link="log"))
summary(gam_seq, dispersion=1)
exp(coef(gam_seq))

yearsim_seq <- seq(from=min(single$Year),
                   to = max(single$Year), 
                   length.out=101)

seq_agbpreds <- predict(gam_seq, 
                        data.frame(Year=yearsim_seq),
                        type="response")

seq_predictions_dev <- data.frame(seq_agbpreds, yearsim_seq)

# Plot the AGB.d and the Gamma models for both the Log_single and Log_sequential models.

plot(sing_predictions$sing_agbpreds~sing_predictions$yearsim_sing, type="l", col="orange", lwd=3,
     main = "Above Ground Biomass, Gamma Model",
     xlab="Year",
     ylab = "AGB, absolute difference from control (gC/m^2)",
     ylim = c(-10,800))
points(seq_predictions_dev$seq_agbpreds~seq_predictions_dev$yearsim_seq, type="l", col="blue", lwd=3)
points(single$AGB.d ~ single$Year, pch=16, col="orange")
points(sequential$AGB.d ~ sequential$Year, pch=16, col="blue")
legend("topright", pch=16, col = c("orange","blue"), legend = c("Single","Sequential"))


# Bayesian Approach: Modeling AGB.d 20, 50, and 100 years post treatment, with and without varying effects.

# The "Null" model: Complete pooling, no varying effects. Can use the stan_glm function, Above Ground Biomass (diff from control) regressed against year, 
# using a Gamma distribution (family) and a log link function (link).
m0.1 <- stan_glm(AGB.d~Year,
               data=AGBpost20,
               family ="Gamma"(link="log")) 
m0.2 <- stan_glm(AGB.d~Year,
                 data=AGBpost50,
                 family ="Gamma"(link="log"))
m0.3 <- stan_glm(AGB.d~Year,
                 data=AGBpost2,
                 family ="Gamma"(link="log"))

# Save the model summaries as a data frame
m0.1_sum <- as.data.frame(summary(m0.1))
m0.2_sum <- as.data.frame(summary(m0.2))
m0.3_sum <- as.data.frame(summary(m0.3))


# Time post treatment: Varying Replicate (changes to CLM(FATES) climate and parameters). Here use stan_glmer for multilevel glms. 
# Distribution (family="gamma") and link="log" remain the same as in model 1. These will remain the same for the other two models.
m1.1 <- stan_glmer(AGB.d~Year + (1|rep), #Similar model syntax as mle glm
             data=AGBpost20,
             family ="Gamma"(link="log"))
m1.2 <- stan_glmer(AGB.d~Year + (1|rep), 
                   data=AGBpost50,
                   family ="Gamma"(link="log"))
m1.3 <- stan_glmer(AGB.d~Year + (1|rep), 
                   data=AGBpost2,
                   family ="Gamma"(link="log"))

# Save the model summaries as a data frame
m1.1_sum <- as.data.frame(summary(m1.1))
m1.2_sum <- as.data.frame(summary(m1.2))
m1.3_sum <- as.data.frame(summary(m1.3))

# Time post treatment: Varying Treatment (changes to scenarios)
m2.1 <-stan_glmer(AGB.d~Year + (1|Treatment), #Similar model syntax as mle glm
               data=AGBpost20,
               family ="Gamma"(link="log"))
m2.2 <-stan_glmer(AGB.d~Year + (1|Treatment), #Similar model syntax as mle glm
                  data=AGBpost50,
                  family ="Gamma"(link="log"))
m2.3 <-stan_glmer(AGB.d~Year + (1|Treatment), #Similar model syntax as mle glm
                  data=AGBpost2,
                  family ="Gamma"(link="log"))

# Save the model summaries as a data frame
m2.1_sum <- as.data.frame(summary(m2.1))
m2.2_sum <- as.data.frame(summary(m2.2))
m2.3_sum <- as.data.frame(summary(m2.3))

# Time post treatment: Varying Treatment + Replicate (parameter and scenarios)
m3.1 <-stan_glmer(AGB.d~Year + (1|Treatment/rep), #Similar model syntax as mle glm
                  data=AGBpost20,
                  family ="Gamma"(link="log"))
m3.2 <-stan_glmer(AGB.d~Year + (1|Treatment/rep), #Similar model syntax as mle glm
                  data=AGBpost50,
                  family ="Gamma"(link="log"))
m3.3 <-stan_glmer(AGB.d~Year + (1|Treatment/rep), #Similar model syntax as mle glm
                  data=AGBpost2,
                  family ="Gamma"(link="log"))

# Save the model summaries as a data frame
m3.1_sum <- as.data.frame(summary(m3.1))
m3.2_sum <- as.data.frame(summary(m3.2))
m3.3_sum <- as.data.frame(summary(m3.3))

# Saves the posteriors from our models
post0.1 <- data.frame(m0.1)
post0.2 <- data.frame(m0.2)
post0.3 <- data.frame(m0.3)


post1.1 <- data.frame(m1.1)
post1.2 <- data.frame(m1.2)
post1.3 <- data.frame(m1.3)

post2.1 <- data.frame(m2.1)
post2.2 <- data.frame(m2.2)
post2.3 <- data.frame(m2.3)

post3.1 <- data.frame(m3.1)
post3.2 <- data.frame(m3.2)
post3.3 <- data.frame(m3.3)

# Although we are using the default weakly informative priors, check to see the prior values used for each model.
prior_summary(m0.1)
prior_summary(m0.2)
prior_summary(m0.3)

prior_summary(m1.1)
prior_summary(m1.2)
prior_summary(m1.3)

prior_summary(m2.1)
prior_summary(m2.2)
prior_summary(m2.3)

prior_summary(m3.1)
prior_summary(m3.2)
prior_summary(m3.3)

# Check model performance by using plot(mod, "trace") to visualize the Monte Carlo Chains completed in the Bayesian models
plot(m0.1, "trace")
plot(m0.2, "trace")
plot(m0.3, "trace")

plot(m1.1, "trace")
plot(m1.2, "trace")
plot(m1.3, "trace")

plot(m2.1, "trace")
plot(m2.2, "trace")
plot(m2.3, "trace")

plot(m3.1, "trace")
plot(m3.2, "trace")
plot(m3.3, "trace")

# Plot the parameter estimates for each of the models. Transform to a probability scale and use "areas" to show the distribution of the parameter values.
# In this case each plot will be saved as a object. Ggplot will be used to plot all of the figure objects as one plot.

m0.1_p <- plot(m0.1, transformations = "plogis",  "areas",  prob=0.95)
m0.2_p <- plot(m0.2, transformations = "plogis",  "areas",  prob=0.95)
m0.3_p <- plot(m0.3, transformations = "plogis",  "areas",  prob=0.95)

# Plot all model plots (figures) as one plot.
figure0 <- ggarrange(m0.1_p, m0.2_p, m0.3_p,
                     labels = c("20 yr post-treatment", "50 yr post-treatment", "100 yr post-treatment"),
                     ncol = 1, nrow=3) 

# View your new figure.
figure0

# Repeat for each of the varying effects models (m1, m2, m3)
# Model 1: Varying Replicate

m1.1_p <- plot(m1.1, transformations = "plogis", pars = "varying", "areas", prob=0.95)
# Can also view as a whisker plot, here the darker inner bar is set to 0.5 probability
plot(m1.1, pars="varying", prob=0.5) 
m1.2_p <- plot(m1.2, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)
#plot(m1.2, pars="varying", prob=0.5)
m1.3_p <- plot(m1.3, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)
#plot(m1.3, pars="varying", prob=0.5)

figure1 <- ggarrange(m1.1_p, m1.2_p, m1.3_p,
                     ncol = 1, nrow=3)
figure1

# Model 2: Varying Treatments
m2.1_p <- plot(m2.1, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)
m2.2_p <- plot(m2.2, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)
m2.3_p <- plot(m2.3, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)

figure2 <- ggarrange(m2.1_p, m2.2_p, m2.3_p,
                     ncol = 1, nrow=3)
figure2

# Model 3: Varying Treatment by Replicate
m3.1_p <- plot(m3.1, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)
m3.2_p <- plot(m3.2, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)
m3.3_p <- plot(m3.3, transformations = "plogis", pars = "varying",  "areas",  prob=0.95)

figure3 <- ggarrange(m3.1_p, m3.2_p, m3.3_p,
                     ncol = 1, nrow=3)

figure3


# Use a Leave One Out (loo) Cross Validation to compare each of the model's predictive performance through time post treatment.


# 100 years post treatment
nullmod.3 <- loo(m0.3)
mod1.3 <- loo(m1.3)
mod2.3 <- loo(m2.3)
fullmod.3 <- loo(m3.3)

loo.3 <- loo_compare(mod1.3, mod2.3, fullmod.3, nullmod.3)

# 50 years post treatment
nullmod.2 <- loo(m0.2)
mod1.2 <- loo(m1.2)
mod2.2 <- loo(m2.2)
fullmod.2 <- loo(m3.2)

loo.2 <- loo_compare(mod1.2, mod2.2, fullmod.2, nullmod.2)

# 20 years post treatment
nullmod.1 <- loo(m0.1)
mod1.1 <- loo(m1.1)
mod2.1 <- loo(m2.1)
fullmod.1 <- loo(m3.1)

loo.1 <- loo_compare(mod1.1, mod2.1, fullmod.1, nullmod.1)

