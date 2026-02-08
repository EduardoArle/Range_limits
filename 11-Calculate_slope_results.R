#load libraries
library(data.table)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables'
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Slopes'
  
#read results table
setwd(wd_tables)
results <- read.csv('20260126_Results_all_sps.csv')

#transform species names in factor
results$species <- as.factor(results$sps)

#select only presences (for now)
results <- results[results$Occurrence == 1,]

#make a species list
sps_list <- unique(results$sps)

#explanatory variables
rangeSize <- numeric()
rangeLoc <- numeric()
roundness <- numeric()
bodyMass <- numeric()
nOcc <- numeric()
order <- character()
elevMedian <- numeric()
elevAmplitude <- numeric()
latAmplitude <- numeric()

#response variables and stats
slope_minT_relPol <- numeric()
MRMSE_minT_relPol <- numeric()
slope_meanT_relPol <- numeric()
MRMSE_meanT_relPol <- numeric()
slope_maxT_relPol <- numeric()
MRMSE_maxT_relPol <- numeric()

slope_minT_distEdge <- numeric()
MRMSE_minT_distEdge <- numeric()
slope_meanT_distEdge <- numeric()
MRMSE_meanT_distEdge <- numeric()
slope_maxT_distEdge <- numeric()
MRMSE_maxT_distEdge <- numeric()

slope_minPPT_relPol <- numeric()
MRMSE_minPPT_relPol <- numeric()
slope_meanPPT_relPol <- numeric()
MRMSE_meanPPT_relPol <- numeric()
slope_maxPPT_relPol <- numeric()
MRMSE_maxPPT_relPol <- numeric()

slope_minPPT_distEdge <- numeric()
MRMSE_minPPT_distEdge <- numeric()
slope_meanPPT_distEdge <- numeric()
MRMSE_meanPPT_distEdge <- numeric()
slope_maxPPT_distEdge <- numeric()
MRMSE_maxPPT_distEdge <- numeric()

#correlation (to exclude highly correlated models from the analysis)
Cor_vars_minT <- numeric()        
Cor_vars_meanT <- numeric()
Cor_vars_maxT <- numeric()
Cor_vars_minPPT <- numeric()
Cor_vars_meanPPT <- numeric()
Cor_vars_maxPPT <- numeric()

#for loop to produce a table with the variables we need to check for species
for(i in 1:length(sps_list))
{
  #select each sps
  res_sps <- results[results$sps == sps_list[i],]
  
  #make table keeping only points up to 250km away from range edges
  res_sps_250 <- res_sps[res_sps$distEdge <= 250,]
  
  #populate the cols with the available explanatory variables
  rangeSize[i] <- unique(res_sps$rangeSize)
  rangeLoc[i] <- unique(res_sps$rangeLoc)
  roundness[i] <- unique(res_sps$roundness)
  bodyMass[i] <- unique(res_sps$bodyMass)
  nOcc[i] <- unique(res_sps$nOcc)
  order[i] <- unique(res_sps$order)
  latAmplitude[i] <- unique(res_sps$latAmpl)
  
  #calculate median elevation
  elevMedian[i] <- median(res_sps$elevation, na.rm = T)
  
  #calculate elevation amplitide (95 % quantile)
  elev_95 <- quantile(res_sps$elevation, probs = c(0.025, 0.975), na.rm = T)
  elevAmplitude[i] <- elev_95[2] - elev_95[1]
 
  #run lms for shap values against relPol and distEdge
  lm_minT_relPol <- try(lm(res_sps$avg_Min_T_SHAP ~ res_sps$relPolewardness,
                           data = res_sps), silent = T)
  lm_meanT_relPol <- try(lm(res_sps$avg_Mean_T_SHAP ~ res_sps$relPolewardness,
                       data = res_sps), silent = T)
  lm_maxT_relPol <- try(lm(res_sps$avg_Max_T_SHAP ~ res_sps$relPolewardness,
                       data = res_sps), silent = T)
  
  lm_minT_distEdge <- try(lm(res_sps_250$avg_Min_T_SHAP ~ res_sps_250$distEdge,
                       data = res_sps_250), silent = T)
  lm_meanT_distEdge <- try(lm(res_sps_250$avg_Mean_T_SHAP ~ res_sps_250$distEdge,
                        data = res_sps_250), silent = T)
  lm_maxT_distEdge <- try(lm(res_sps_250$avg_Max_T_SHAP ~ res_sps_250$distEdge,
                       data = res_sps_250), silent = T)
  
  lm_minPPT_relPol <- try(lm(res_sps$avg_Min_PPT_SHAP ~ res_sps$relPolewardness,
                       data = res_sps), silent = T)
  lm_meanPPT_relPol <- try(lm(res_sps$avg_Mean_PPT_SHAP ~ res_sps$relPolewardness,
                        data = res_sps), silent = T)
  lm_maxPPT_relPol <- try(lm(res_sps$avg_Max_PPT_SHAP ~ res_sps$relPolewardness,
                       data = res_sps), silent = T)
  
  lm_minPPT_distEdge <- try(lm(res_sps_250$avg_Min_PPT_SHAP ~
                               res_sps_250$distEdge,
                         data = res_sps_250), silent = T)
  lm_meanPPT_distEdge <- try(lm(res_sps_250$avg_Mean_PPT_SHAP ~
                                res_sps_250$distEdge,
                          data = res_sps_250), silent = T)
  lm_maxPPT_distEdge <- try(lm(res_sps_250$avg_Max_PPT_SHAP ~
                               res_sps_250$distEdge,
                         data = res_sps_250), silent = T)
  
  #calculate slope and MRMSE of lms
  
  #T relPol
  if(class(lm_minT_relPol) == 'lm'){
    slope_minT_relPol[i] <- coef(lm_minT_relPol)[2]
    
    pred <- predict(lm_minT_relPol)
    rmse_minT_relPol <- sqrt(mean((res_sps$avg_Min_T_SHAP - pred) ^ 2))
    MRMSE_minT_relPol[i] <-  rmse_minT_relPol / mean(res_sps$avg_Min_T_SHAP)
  }else{
    slope_minT_relPol[i] <- NA
    MRMSE_minT_relPol[i] <- NA
  }

  if(class(lm_meanT_relPol) == 'lm'){
    slope_meanT_relPol[i] <- coef(lm_meanT_relPol)[2]
    
    pred <- predict(lm_meanT_relPol)
    rmse_meanT_relPol <- sqrt(mean((res_sps$avg_Mean_T_SHAP - pred) ^ 2))
    MRMSE_meanT_relPol[i] <- rmse_meanT_relPol / mean(res_sps$avg_Mean_T_SHAP)
  }else{
    slope_meanT_relPol[i] <- NA
    MRMSE_meanT_relPol[i] <- NA
  }

  if(class(lm_maxT_relPol) == 'lm'){
    slope_maxT_relPol[i] <- coef(lm_maxT_relPol)[2]
    
    pred <- predict(lm_maxT_relPol)
    rmse_maxT_relPol <- sqrt(mean((res_sps$avg_Max_T_SHAP - pred) ^ 2))
    MRMSE_maxT_relPol[i] <- rmse_maxT_relPol / mean(res_sps$avg_Max_T_SHAP)
  }else{
    slope_maxT_relPol[i] <- NA
    MRMSE_maxT_relPol[i] <- NA
  }
  
  
  #T distEdge
  if(class(lm_minT_distEdge) == 'lm'){
    slope_minT_distEdge[i] <- coef(lm_minT_distEdge)[2]
    
    pred <- predict(lm_minT_distEdge)
    rmse_minT_distEdge <- sqrt(mean((res_sps$avg_Min_T_SHAP - pred) ^ 2))
    MRMSE_minT_distEdge[i] <-  rmse_minT_distEdge / mean(res_sps$avg_Min_T_SHAP)
  }else{
    slope_minT_distEdge[i] <- NA
    slope_minT_distEdge[i] <- NA
  }
 
  if(class(lm_meanT_distEdge) == 'lm'){
    slope_meanT_distEdge[i] <- coef(lm_meanT_distEdge)[2]
      
    pred <- predict(lm_meanT_distEdge)
    rmse_meanT_distEdge <- sqrt(mean((res_sps$avg_Mean_T_SHAP - pred) ^ 2))
    MRMSE_meanT_distEdge[i] <- rmse_meanT_distEdge / mean(res_sps$avg_Mean_T_SHAP)
  }else{
    slope_meanT_distEdge[i] <- NA
    MRMSE_meanT_distEdge[i] <- NA
  }
  
  if(class(lm_maxT_distEdge) == 'lm'){
    slope_maxT_distEdge[i] <- coef(lm_maxT_distEdge)[2]
    
    pred <- predict(lm_maxT_distEdge)
    rmse_maxT_distEdge <- sqrt(mean((res_sps$avg_Max_T_SHAP - pred) ^ 2))
    MRMSE_maxT_distEdge[i] <- rmse_maxT_distEdge / mean(res_sps$avg_Max_T_SHAP)
  }else{
    slope_maxT_distEdge[i] <- NA
    MRMSE_maxT_distEdge[i] <- NA
  }
 
  
  #PPT relPol
  if(class(lm_minPPT_relPol) == 'lm'){
    slope_minPPT_relPol[i] <- coef(lm_minPPT_relPol)[2]
    
    pred <- predict(lm_minPPT_relPol)
    rmse_minPPT_relPol <- sqrt(mean((res_sps$avg_Min_PPT_SHAP - pred) ^ 2))
    MRMSE_minPPT_relPol[i] <-  rmse_minPPT_relPol / mean(res_sps$avg_Min_PPT_SHAP)
  }else{
    slope_minPPT_relPol[i] <- NA
    MRMSE_minPPT_relPol[i] <- NA
  }

  if(class(lm_meanPPT_relPol) == 'lm'){
    slope_meanPPT_relPol[i] <- coef(lm_meanPPT_relPol)[2]
    
    pred <- predict(lm_meanPPT_relPol)
    rmse_meanPPT_relPol <- sqrt(mean((res_sps$avg_Mean_PPT_SHAP - pred) ^ 2))
    MRMSE_meanPPT_relPol[i] <- rmse_meanPPT_relPol / mean(res_sps$avg_Mean_PPT_SHAP)
  }else{
    slope_meanPPT_relPol[i] <- NA
    MRMSE_meanPPT_relPol[i] <- NA
  }
 
  if(class(lm_maxPPT_relPol) == 'lm'){
    slope_maxPPT_relPol[i] <- coef(lm_maxPPT_relPol)[2]
    
    pred <- predict(lm_maxPPT_relPol)
    rmse_maxPPT_relPol <- sqrt(mean((res_sps$avg_Max_PPT_SHAP - pred) ^ 2))
    MRMSE_maxPPT_relPol[i] <- rmse_maxPPT_relPol / mean(res_sps$avg_Max_PPT_SHAP)
  }else{
    slope_maxPPT_relPol[i] <- NA
    MRMSE_maxPPT_relPol[i] <- NA
  }
 
  
  #PPT distEdge
  if(class(lm_minPPT_distEdge) == 'lm'){
    slope_minPPT_distEdge[i] <- coef(lm_minPPT_distEdge)[2]
    
    pred <- predict(lm_minPPT_distEdge)
    rmse_minPPT_distEdge <- sqrt(mean((res_sps$avg_Min_PPT_SHAP - pred) ^ 2))
    MRMSE_minPPT_distEdge[i] <- rmse_minPPT_distEdge /mean(res_sps$avg_Min_PPT_SHAP)
  }else{
    slope_minPPT_distEdge[i] <- NA
    MRMSE_minPPT_distEdge[i] <- NA
  }
  
  if(class(lm_meanPPT_distEdge) == 'lm'){
    slope_meanPPT_distEdge[i] <- coef(lm_meanPPT_distEdge)[2]
    
    pred <- predict(lm_meanPPT_distEdge)
    rmse_meanPPT_distEdge <- sqrt(mean((res_sps$avg_Mean_PPT_SHAP - pred) ^ 2))
    MRMSE_meanPPT_distEdge[i] <- rmse_meanPPT_distEdge / mean(res_sps$avg_Mean_PPT_SHAP)
  }else{
    slope_meanPPT_distEdge[i] <- NA
    MRMSE_meanPPT_distEdge[i] <- NA
  }
  
  if(class(lm_maxPPT_distEdge) == 'lm'){
    slope_maxPPT_distEdge[i] <- coef(lm_maxPPT_distEdge)[2]
    
    pred <- predict(lm_maxPPT_distEdge)
    rmse_maxPPT_distEdge <- sqrt(mean((res_sps$avg_Max_PPT_SHAP - pred) ^ 2))
    MRMSE_maxPPT_distEdge[i] <- rmse_maxPPT_distEdge /mean(res_sps$avg_Max_PPT_SHAP)
  }else{
    slope_maxPPT_distEdge[i] <- NA
    MRMSE_maxPPT_distEdge[i] <- NA
  }
 
  
  #populate vectors with correlation values
  Cor_vars_minT[i] <- unique(res_sps$Cor_vars_minT)
  Cor_vars_meanT[i] <- unique(res_sps$Cor_vars_meanT)
  Cor_vars_maxT[i] <- unique(res_sps$Cor_vars_maxT)
  Cor_vars_minPPT[i] <- unique(res_sps$Cor_vars_minPPT)
  Cor_vars_meanPPT[i] <- unique(res_sps$Cor_vars_meanPPT)
  Cor_vars_maxPPT[i] <- unique(res_sps$Cor_vars_maxPPT)

  print(i)
}


#create a data frame with results
slopes <- data.frame(species = sps_list, rangeSize = rangeSize,
                     rangeLoc = rangeLoc, roundness = roundness, 
                     bodyMass = bodyMass, nOcc = nOcc, order = order,
                     elevMedian = elevMedian, elevAmplitude = elevAmplitude,
                     latAmplitude = latAmplitude,
                     slope_minT_relPol = slope_minT_relPol,
                     MRMSE_minT_relPol = MRMSE_minT_relPol,
                     slope_meanT_relPol = slope_meanT_relPol,
                     MRMSE_meanT_relPol = MRMSE_meanT_relPol,
                     slope_maxT_relPol = slope_maxT_relPol,
                     MRMSE_maxT_relPol = MRMSE_maxT_relPol,
                     slope_minT_distEdge = slope_minT_distEdge,
                     MRMSE_minT_distEdge = MRMSE_minT_distEdge, 
                     slope_meanT_distEdge = slope_meanT_distEdge,
                     MRMSE_meanT_distEdge = MRMSE_meanT_distEdge,
                     slope_maxT_distEdge = slope_maxT_distEdge,
                     MRMSE_maxT_distEdge = MRMSE_maxT_distEdge,
                     slope_minPPT_relPol = slope_minPPT_relPol,
                     MRMSE_minPPT_relPol = MRMSE_minPPT_relPol,
                     slope_meanPPT_relPol = slope_meanPPT_relPol,
                     MRMSE_meanPPT_relPol = MRMSE_meanPPT_relPol,
                     slope_maxPPT_relPol = slope_maxPPT_relPol,
                     MRMSE_maxPPT_relPol = MRMSE_maxPPT_relPol,
                     slope_minPPT_distEdge = slope_minPPT_distEdge,
                     MRMSE_minPPT_distEdge = MRMSE_minPPT_distEdge,
                     slope_meanPPT_distEdge = slope_meanPPT_distEdge,
                     MRMSE_meanPPT_distEdge = MRMSE_meanPPT_distEdge,
                     slope_maxPPT_distEdge = slope_maxPPT_distEdge,
                     MRMSE_maxPPT_distEdge = MRMSE_maxPPT_distEdge,
                     Cor_vars_minT = Cor_vars_minT,       
                     Cor_vars_meanT = Cor_vars_meanT,
                     Cor_vars_maxT = Cor_vars_maxT,
                     Cor_vars_minPPT = Cor_vars_minPPT,
                     Cor_vars_meanPPT = Cor_vars_meanPPT,
                     Cor_vars_maxPPT = Cor_vars_maxPPT)

#save table with slopes
setwd(wd_slopes)
write.csv(slopes, '20260208_Slopes.csv', row.names = F)
