
EE_analysis= function(data){

## @data is in the long format
## "time" is the time variable for analysis. For example time in the heat data
## for the purpose of simple application, we need to change the time variable in the data into "time" if it was not. 
  
require(stats)
require(lattice)
require(magrittr)
require(Matrix)
require(MASS)
require(reshape)
require(tidyverse)
require(mi) ## for imputation
require(lmerTest) ## for summary() with p-value
require(cplm)
require(optimx) ## to chose different optimer methods than the default ones in lmer package
require(dplyr)

  
####--------------------------------------------------------------------------------------------####
#### linear mixed effect models (raw heat)
####--------------------------------------------------------------------------------------------####

### linear mixed effect model #####
mod.lm = lmer(Heat~ trtgrp + time +(1|animal), data = data)
  
### linear quadratic mixed effect model ####
mod.qlm = lmer(Heat~ trtgrp + time +time2+(1|animal), data = data)


### linear cubic mixed effect model ####
mod.clm = lmer(Heat~ trtgrp + time+time2+time3 +(1|animal), data = data)

### put all the linea-mixed effect model in a list ##
mod.pool = list(mod.lm, mod.qlm, mod.clm)

## get aic for each of them 
aic.pool = lapply(mod.pool, function(s){
  AIC(logLik(s))
})

### choose the one with smallest aic ###
mod.lmm = mod.pool[[which.min(aic.pool)]]


####--------------------------------------------------------------------------------------------####
#### linear mixed effect models (log heat)
####--------------------------------------------------------------------------------------------####

### linear mixed effect model #####
mod.lm.log = lmer(log(Heat)~ trtgrp + time +(1|animal), data = data)
  
### linear quadratic mixed effect model ####
mod.qlm.log = lmer(log(Heat)~ trtgrp + time +time2+(1|animal), data = data)

### linear cubic mixed effect model ####
mod.clm.log = lmer(log(Heat)~ trtgrp + time+time2+time3 +(1|animal), data = data)

### put all the linea-mixed effect model in a list ##
mod.pool = list(mod.lm.log, mod.qlm.log, mod.clm.log)

## get aic for each of them 
aic.pool = lapply(mod.pool, function(s){
  AIC(logLik(s))
})

### choose the one with smallest aic ###
mod.lmm.log = mod.pool[[which.min(aic.pool)]]




####--------------------------------------------------------------------------------------------####
#### B-splines model (raw)
####--------------------------------------------------------------------------------------------####

##### ##### ##### 
##### linear ######
##### ##### ##### 

### get the AIC for different number (5:8) of b-splines
aic.bs = sapply(5:8, function(s) {
  K_AIC.bs(s, data, data$time, model="linear")
})

### chose the k (the number of basis) that has the smallest AIC
k = (5:8)[which.min(aic.bs)]

### generate the linear (one-degree) basis 
bs3.time = bs(data$time, df = k, degree = 1, intercept = F) ## 432*(k)
colnames(bs3.time) = paste0("bs", 1:k)
data.bs = cbind(data, bs3.time) ## combine the orginal data with the bsplines basis for data analysis 

### create random effect part of the formula
## starting with 1, and additional bs term will be added each time until k times
## for example, for the first time, bs1 will be added to 1 and bs2 will be added to 1+bs1 for the second time. 
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

## create the full formula
## add the fixed effect to the formula for random effect
formula = paste0("Heat~ trtgrp + time + (", rf, "|animal)")

### semi-parametric linear mixed effect model with linear time
mod.bs.lm = lmer(as.formula(formula), data = data.bs, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
 

##### ##### ##### 
##### quadratic ######
##### ##### ##### 

aic.bs = sapply(5:8, function(s) {
  K_AIC.bs(s, data, data$time, model="quadratic")
})

k = (5:8)[which.min(aic.bs)]

## quadraticc bspline for time(mins)
bs3.time = bs(data$time, df = k, degree = 2, intercept = F) ## 432*(k)
colnames(bs3.time) = paste0("bs", 1:k)
data.bs = cbind(data, bs3.time)

### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
formula = paste0("Heat~ trtgrp + time + time2+(", rf, "|animal)")

mod.bs.qm = lmer(as.formula(formula), data = data.bs, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
 
##### ##### ##### 
##### cubic ######
##### ##### ##### 

aic.bs = sapply(5:8, function(s) {
  K_AIC.bs(s, data, data$time, model="cubic")
})

k = (5:8)[which.min(aic.bs)]

## cubic bspline for time(mins)
bs3.time = bs(data$time, df = k, degree = 3, intercept = F) ## 432*(k)
colnames(bs3.time) = paste0("bs", 1:k)
data.bs = cbind(data, bs3.time)

### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
formula = paste0("Heat~ trtgrp + time +time2 +time3 +(", rf, "|animal)")

mod.bs.cm = lmer(as.formula(formula), data = data.bs, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
 

##### ##### ##### 
##### linear ######
##### ##### ##### 
### select the number of basis based on AIC of mixed effect model 
aic.tp = sapply(5:8, function(s) {
  K_AIC.tp(s, data,data$time, model="linear")
}) 

k = (5:8)[which.min(aic.tp)]

## linear bspline for time
bs3.time = tp(data$time, k = k+1, degree = 1, allPen = F)$Z ## cubic basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
colnames(bs3.time) = paste0("bs", 1:k)
data.tp = cbind(data, bs3.time)
  
### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time + (", rf, "|animal)")

mod.tp.lm = lmer(as.formula(formula), data = data.tp, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
  ### use REML to avoid sigularity issues 
 ### see https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
  

##### ##### ##### 
##### quadratic ######
##### ##### ##### 
### select the number of basis based on AIC of mixed effect model 
aic.tp = sapply(5:8, function(s) {
  K_AIC.tp(s, data,data$time, model="quadratic")
}) 

k = (5:8)[which.min(aic.tp)]

## quadratic bspline for time

bs3.time = tp(data$time, k = k+2, degree = 2, allPen = F)$Z ## cubic basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
colnames(bs3.time) = paste0("bs", 1:k)
data.tp = cbind(data, bs3.time)
  
### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time + time2+(", rf, "|animal)")

mod.tp.qm = lmer(as.formula(formula), data = data.tp, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
  ### use REML to avoid sigularity issues 
 ### see https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
  

##### ##### ##### 
##### cubic ######
##### ##### ##### 
### select the number of basis based on AIC of mixed effect model 
aic.tp = sapply(5:8, function(s) {
  K_AIC.tp(s, data,data$time, model="cubic")
}) 

k = (5:8)[which.min(aic.tp)]

## cubic bspline for time
bs3.time = tp(data$time, k = k+3, degree = 3, allPen = F)$Z ## cubic basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
colnames(bs3.time) = paste0("bs", 1:k)
data.tp = cbind(data, bs3.time)
  
### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time + time2+time3+(", rf, "|animal)")

mod.tp.cm = lmer(as.formula(formula), data = data.tp, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
 ### use REML to avoid sigularity issues 
 ### see https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
  

 return(list(mod.lm=mod.lm, mod.lm.log=mod.lm.log,
             mod.qlm=mod.qlm, mod.qlm.log=mod.qlm.log, 
             mod.clm=mod.clm,mod.clm.log=mod.clm.log,
             mod.bs.lm=mod.bs.lm,  mod.tp.lm=mod.tp.lm, 
             mod.bs.qm=mod.bs.qm, mod.tp.qm=mod.tp.qm,
             mod.bs.cm=mod.bs.cm, mod.tp.cm=mod.tp.cm))


}

# dat = read.table("../../rats_heat.txt", header=T)
# EE_analysis.min(dat)
