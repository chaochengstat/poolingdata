### Set the working directory
rm(list=ls())
setwd("/Users/chaocheng/Dropbox/Pooling Paper 2/for JASA/pseudo dataset")

### Load the R packages and functions that contain the ACM and ECM estimation algorithms
library("lme4");
library("nlme");
library("pracma");
library("statmod");
library("mvtnorm");
library("fGarch");
source("functions_CCS.R")

### Load the pooling 25(OH)D dataset
pooled_data=read.csv("pooled_data.csv")[,-1]
head(pooled_data)

### Because the 25(OH)D presents seasonal flucations, we use a periodic function of the
### week of blood draw to represent this seasonal trend
pooled_data[,"T1"]=sin(2*pi*pooled_data[,"draww"]/52)
pooled_data[,"T2"]=cos(2*pi*pooled_data[,"draww"]/52)
pooled_data[,"T3"]=sin(4*pi*pooled_data[,"draww"]/52)
pooled_data[,"T4"]=cos(4*pi*pooled_data[,"draww"]/52)
pooled_data=as.matrix(pooled_data)

#################################################################################
#          Model I
#################################################################################
# 1, Covariates of the model for the 25(OH)D:
# Week of blood draw, Physical activity, Age of blood draw, Smoking, BMI.
# 2, Covariates of the model for the disease outcome (colorectal cancer):
# physical activity, family history of colorectal cancer, smoking, BMI.
#################################################################################

result =main_func(data=pooled_data,             
                  Y="Y",                      ## disease outcome
                  HL="HL",                    ## local lab measurements
                  HC="HC",                    ## reference lab measurements 
                  Wname=c("T1","T2","T3","T4","physical","ageblood","smoke","bmi0"),
                  Zname=c("physical","family","smoke","bmi0"),
                  ID="matchID",               ## ID for the matched set or strata
                  group="group",              ## group or study number
                  calibration="calibration"   ## in the calib. subset?
)

#### (1) calibration parameters
## fixed effects
result$fixed_effect
## gamma estimates
result$gamma
## xi estimates
result$xi
## sigmas
result$sigma
#### (2) log Odds Ratios
## point estimates of log(OR)
result$beta_point
## SE of the log(OR) estimates 
result$beta_se


#################################################################################
#          Model II
#################################################################################
# 1, Covariates of the model for the 25(OH)D:
# Week of blood draw.
# 2, Covariates of the model for the disease outcome (colorectal cancer):
# physical activity, family history of colorectal cancer, smoking, BMI.
#################################################################################

result =main_func(data=pooled_data,             
                  Y="Y",                      ## disease outcome
                  HL="HL",                    ## local lab measurements
                  HC="HC",                    ## reference lab measurements 
                  Wname=c("T1","T2","T3","T4"),
                  Zname=c("physical","family","smoke","bmi0"),
                  ID="matchID",               ## ID for the matched set or strata
                  group="group",              ## group or study number
                  calibration="calibration"   ## in the calib. subset?
)

#### (1) calibration parameters
## fixed effects
result$fixed_effect
## gamma estimates
result$gamma
## xi estimates
result$xi
## sigmas
result$sigma
#### (2) log Odds Ratios
## point estimates of log(OR)
result$beta_point
## SE of the log(OR) estimates 
result$beta_se


#################################################################################
#          Model III
#################################################################################
# 1, Covariates of the model for the 25(OH)D:
# None covariates are used.
# 2, Covariates of the model for the disease outcome (colorectal cancer):
# physical activity, family history of colorectal cancer, smoking, BMI.
#################################################################################

result =main_func(data=pooled_data,             
                  Y="Y",                      ## disease outcome
                  HL="HL",                    ## local lab measurements
                  HC="HC",                    ## reference lab measurements 
                  Wname=NULL,
                  Zname=c("physical","family","smoke","bmi0"),
                  ID="matchID",               ## ID for the matched set or strata
                  group="group",              ## group or study number
                  calibration="calibration"   ## in the calib. subset?
)

#### (1) calibration parameters
## fixed effects
result$fixed_effect
## gamma estimates
result$gamma
## xi estimates
result$xi
## sigmas
result$sigma
#### (2) log Odds Ratios
## point estimates of log(OR)
result$beta_point
## SE of the log(OR) estimates 
result$beta_se
