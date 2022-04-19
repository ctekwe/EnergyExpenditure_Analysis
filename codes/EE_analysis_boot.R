EE_analysis_boot = function(dat,boot_iter = 500, computing="parallel"){
################ set up parameters for bootstrap #########
 ## @boot_iter is number of iterations (replications) of bootstrap
 ## @dat is the sample data in the long format (the same one used for data analysis)
 ## for the purpose of simple application, we need to name the variable for subject id "animal" if it was not. 
 ## @computing is to select wether parallel computing or not. Defualt is paralleling computing
  
ids = unique(dat$animal) ### id of subjects 
n = length(ids) ## sample size

#### We can either use parallel computing or non-parallel computing 

if (computing=="parallel"){
#################### parallel computing ############
require(parallel)
cores = detectCores() -1

mclapply(1:boot_iter, function(s){
    #print(s)
  set.seed(s+seeds)
  #### select the id of subject randomly with replacement
 id.boot= sample(ids,n,replace = T) %>% as.data.frame()
 colnames(id.boot) ="animal" 
 dat.boot = merge(dat,id.boot, by ="animal") %>%  #analysis at time level
   ### current data with repeated animal ids which will be considered as one subject
   ### during the analysis and thus we need to rename the subject
   dplyr::group_by(time) %>%  ## for the same time, rename the repeated animal (subject)
  mutate(animal = ave(as.character(animal), animal, FUN= function(x)
    if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])) %>%
  ungroup() %>%
  as.data.frame()

  EE_analysis(dat.boot) ### analysis
  
},mc.cores = cores)} else{
#################### nonparallel computing ############
lapply(1:boot_iter, function(s){
    #print(s)
  set.seed(s+seeds)
  #### select the id of subject randomly with replacement
 id.boot= sample(ids,n,replace = T) %>% as.data.frame()
 colnames(id.boot) ="animal" 
 dat.boot = merge(dat,id.boot, by ="animal") %>%  #analysis at time level
   ### current data with repeated animal ids which will be considered as one subject
   ### during the analysis and thus we need to rename the subject
   dplyr::group_by(time) %>%  ## for the same time, rename the repeated animal (subject)
  mutate(animal = ave(as.character(animal), animal, FUN= function(x)
    if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])) %>%
  ungroup() %>%
  as.data.frame()

  EE_analysis(dat.boot) ### analysis
})
}
}