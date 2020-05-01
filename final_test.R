## Testing:
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
ncname <- "cntrl_AGBv01"
ncfname <- paste(ncpath, ncname, ".nc", sep="" )
dname <- "AGB" 

# Open a netCDF file

ncin <- nc_open(ncfname)
print(ncin)

# get time

time <- ncvar_get(ncin,"time")
time

tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt
tunits

agb_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin, dname,"long name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(agb_array)

# convert time -- split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
chron(time,origin=c(tmonth,tday,tyear))

# Replace netCDF fill values with NAs

agb_array[agb_array==fillvalue$value] <- NA

length(na.omit(as.vector(agb_array[,1]))) # 12 months in a year, AGB is a monthly averaged output


# create dataframe -- reshape data 

# Make a slice of the data, here I am taking from year 100

y <- 100
agb_slice <- agb_array[,y]

# I feel like I need to do something with time...
timemat <- as.matrix(time)
length(timemat)

dim(agb_slice)

plot(agb_slice)

# Now, back to reshaping, create a vector from your slice

agb_vec <- as.vector(agb_slice)
length(agb_vec)

agb_df01 <- data.frame(agb_vec)
names(agb_df01) <- c(paste(dname, as.character(y), sep = "_"))
head(na.omit(agb_df01), 12)

# Write out the data
# set path and filename
csvpath <- "data_out/"
csvname <- "cntrl_agb_test_yr100.csv"
csvfile <- paste(csvpath, csvname, sep="")
write.table(na.omit(agb_df01), csvfile, row.names=FALSE, sep=",")


# Convert the whole array to a data frame, calculate annual mean AGB

# Reshape the entire array
# Reshape the array into vector

agb_vec_long <- as.vector(agb_array)
length(agb_vec_long)

agb_mat <- matrix(agb_vec_long, nrow=nt)
dim(agb_mat)

head(na.omit(agb_mat))

agb_df02 <- data.frame(cbind(time,agb_mat))
head(na.omit(agb_df02))

test_mat <- na.omit(agb_df02)

### Something itsn't working... hopefully not because of how I subset from python.
# I'll work through this again but with the data that has already averaged for the year

# Set path and filename

ncname2 <- "cntrl_AGB"
ncfname2 <- paste(ncpath, ncname2, ".nc", sep="")
#dname2 <- "AGB"

# open a netCDF file

ncin2 <- nc_open(ncfname2)
print(ncin2)
print(ncin)

time2 <- ncvar_get(ncin2,"year")
time2

tunits2 <- ncatt_get(ncin2, "year","units")
nt2 <- dim(time2)
nt2
tunits2
agb2_array <- ncvar_get(ncin2,dname)
dim(agb2_array)
agb2_array

plot(agb2_array)

year <- time2
agb2 <- agb2_array[]

plot(year,agb2)




