### select k based on AIC of the model  for b-splines
### function for calulating AIC for model with different k values
K_AIC.bs = function(k, data, X, model){
  ### k is the number of basis and should be at least 4 and at most 9
  ### data is the data used for regression 
  ### X is the data for basis 
  ### model: linear, quadratic, or cubic model 
  if(model=="linear"){
    bs3 = bs(X, df = k, degree = 1, intercept = F) ## linear basis
    colnames(bs3) = paste0("bs", 1:k)
    data = cbind(data, bs3) ## create dataframe for regression model 
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

  }else if(model=="quadratic"){
    bs3 = bs(X, df = k, degree = 2, intercept = F) ## quadratic basis
    colnames(bs3) = paste0("bs", 1:k)
    data= cbind(data, bs3) ## create dataframe for regression model 
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time +time2 + (", rf, "|animal)")

  }else{
    bs3 = bs(X, df = k, degree = 3, intercept = F) ## cubic basis
    colnames(bs3) = paste0("bs", 1:k)
    data= cbind(data, bs3) ## create dataframe for regression model 
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time +time2+time3 +(", rf, "|animal)")

  }

 mod = lmer(as.formula(formula), data = data, REML = T,control = lmerControl(
 optimizer ='Nelder_Mead'))
  # mod = lmer(as.formula(formula), data = data.bs, REML = T,control = lmerControl(
  #                          optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

  #mod = lmer(as.formula(formula), data = data, REML = F)
           
 ### use REML to avoid sigularity issues 
 ### see https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
  
  aic = AIC(logLik(mod))
  return(aic)
}



#### select k for log(Heat)
K_AIC.bs.log = function(k, data, X, model){
  ### k is the number of basis and should be at least 4 and at most 9
  ### data is the data used for regression 
  ### X is the data for basis 
  ### model: linear, quadratic, or cubic model 
  if(model=="linear"){
    bs3 = bs(X, df = k, degree = 1, intercept = F) ## linear basis
    colnames(bs3) = paste0("bs", 1:k)
    data = cbind(data, bs3) ## create dataframe for regression model 
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("log(Heat)~ trtgrp + time + (", rf, "|animal)")

  }else if(model=="quadratic"){
    bs3 = bs(X, df = k, degree = 2, intercept = F) ## quadratic basis
    colnames(bs3) = paste0("bs", 1:k)
    data = cbind(data, bs3) ## create dataframe for regression model 
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("log(Heat)~ trtgrp + time +time2 + (", rf, "|animal)")

  }else{
    bs3 = bs(X, df = k, degree = 3, intercept = F) ## cubic basis
    colnames(bs3) = paste0("bs", 1:k)
    data = cbind(data, bs3) ## create dataframe for regression model 
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<k+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("log(Heat)~ trtgrp + time +time2+time3 +(", rf, "|animal)")

  }

 mod = lmer(as.formula(formula), data = data, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
 ### use REML to avoid sigularity issues 
 ### see https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
  
  aic = AIC(logLik(mod))
  return(aic)
}

### select k based on AIC of the model  for truncated power basis
K_AIC.tp = function(knots, data, X, model){
  ### knots is the number of basis and should be at least 4 and at most 9
  ### data is the data for mixed effect regression 
  ### X is he covariate for the smooth function 
   ### model: linear, quadratic, or cubic model 
  
  if(model=="linear"){
    bs3 = tp(X, k=knots+1, degree = 1, allPen = F)$Z ## linear basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
  colnames(bs3) = paste0("bs", 1:knots)
  ## Z is an n * (k-degree) design matrix for penalized part
  data = cbind(data, bs3) ## create dataframe for regression model 
  
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<knots+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time + (", rf, "|animal)")
  }else if(model=="quadratic"){
    bs3 = tp(X, k=knots+2, degree = 2, allPen = F)$Z ## quadratic basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
  colnames(bs3) = paste0("bs", 1:knots)
  ## Z is an n * (k-degree) design matrix for penalized part
  data = cbind(data, bs3) ## create dataframe for regression model 
  
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<knots+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time+ time2+ (", rf, "|animal)")

  }else{
    bs3 = tp(X, k=knots+3, degree = 3, allPen = F)$Z ## cubic basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
  colnames(bs3) = paste0("bs", 1:knots)
  ## Z is an n * (k-degree) design matrix for penalized part
  data = cbind(data, bs3) ## create dataframe for regression model 
  
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<knots+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("Heat~ trtgrp + time + time2+time3+(", rf, "|animal)")

  }
 mod = lmer(as.formula(formula), data = data, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
  ### use REML to avoid sigularity issues 
 ### see https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
 
  aic = AIC(logLik(mod))
  return(aic)
}

K_AIC.tp.log = function(knots, data, X, model){
  ### knots is the number of basis and should be at least 4 and at most 9
  ### data is the data for mixed effect regression 
  ### X is he covariate for the smooth function 
   ### model: linear, quadratic, or cubic model 
  
  if(model=="linear"){
    bs3 = tp(X, k=knots+1, degree = 1, allPen = F)$Z ## linear basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
  colnames(bs3) = paste0("bs", 1:knots)
  ## Z is an n * (k-degree) design matrix for penalized part
  data = cbind(data, bs3) ## create dataframe for regression model 
  
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<knots+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("log(Heat)~ trtgrp + time + (", rf, "|animal)")
  }else if(model=="quadratic"){
    bs3 = tp(X, k=knots+2, degree = 2, allPen = F)$Z ## quadratic basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
  colnames(bs3) = paste0("bs", 1:knots)
  ## Z is an n * (k-degree) design matrix for penalized part
  data = cbind(data, bs3) ## create dataframe for regression model 
  
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<knots+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("log(Heat)~ trtgrp + time+ time2+ (", rf, "|animal)")

  }else{
    bs3 = tp(X, k=knots+3, degree = 3, allPen = F)$Z ## cubic basis 
  ## k in tp() is the dimensionality of the basis k = number of knots + degree
  ## we set k = knots+3 so that k will be the "knots" number of basis
  colnames(bs3) = paste0("bs", 1:knots)
  ## Z is an n * (k-degree) design matrix for penalized part
  data = cbind(data, bs3) ## create dataframe for regression model 
  
  ### create random effect part of the formula
  i = 1
  temp = 1
  while (i<knots+1){
  rf = paste0(temp,"+" ,paste0("bs", i))
  temp = rf
  i = i+1
} ### end of while loop

  ## create the full formula
  formula = paste0("log(Heat)~ trtgrp + time + time2+time3+(", rf, "|animal)")

  }
 mod = lmer(as.formula(formula), data = data, REML = T,control = lmerControl(  optimizer ='Nelder_Mead'))
  ### use REML to avoid sigularity issues 
 ### see https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
 
  aic = AIC(logLik(mod))
  return(aic)
}
