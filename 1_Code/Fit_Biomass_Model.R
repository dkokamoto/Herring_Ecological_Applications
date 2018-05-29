###########################################
### Basic Code to Fit the Biomass Model ###
### Author: DK Okamoto                  ###
###########################################


library(ggplot2);library(reshape2);
library(gridExtra);library(scales);library(grid)
library(zoo);library(dplyr);library(gdata);library(rstan)
library(lubridate)

### read in data
biomass <- read.csv("2_Data/CC_BIOMASS_Data.csv")

### turn data into matrix
data_mat <- acast(biomass, SECTION ~ YEAR ~ variable)

### create a numeric scaler so biomass is on a manageable scale 
### (mean sd in observed total biomass)
scaler <- mean(apply(data_mat[,,"TOT_BIOM"],2,sd,na.rm=T))

### create a numeric scaler so biomass is on a manageable scale
data_test <- t(data_mat[,,"SP_BIOMASS"]/scaler)
data_test2 <- t(data_mat[,,"TOTAL_CATCH"]/scaler)

### set observed NAs to zero (for catch) or the observed minimum
data_test[is.na(data_test)] <- min(data_test[data_test>0&!is.na(data_test)])
data_test2[is.na(data_test2)] <- 0

### inputs (p and q)
p=3 ### reproductive lag
q= 1 ### order of moving average in errors 
NS= length(unique(biomass$SECTION)) ### number of sites
T= length(unique(biomass$YEAR)) ### number of time periods

### create a list of data for Stan
stan_data <- list(
  T=T,
  NS=NS,
  p=3,
  S_obs=data_test, 
  H= data_test2,
  q=1
)

### fit the model ###
fit_q <- stan(file= "1_Code/Varma_q.stan", 
                data=stan_data, iter=2000,warmup =1000,chains= 5,cores= 5)
