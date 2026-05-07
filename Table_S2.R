#load libraries
library(mgcv); library(itsadug)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Models'
wd_result <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SI/Tables'

#empty vectors
AIC_Model_G <- numeric()
AIC_Model_GS <- numeric()
DE_Model_G <- numeric()
DE_Model_GS <- numeric()

#get models names
setwd(wd_models)
model_files <- list.files()

#get unique model pair names
names_models <- unique(gsub('_G|_GS', '', model_files))

#loop through model pairs
for(i in 1:length(names_models))
{
  #load G model
  mod_G <- readRDS(paste0(names_models[i], '_G'))
  
  AIC_Model_G[i] <- AIC(mod_G)
  DE_Model_G[i] <- summary(mod_G)$dev.expl
  
  rm(mod_G)
  gc()
  
  #load GS model
  mod_GS <- readRDS(paste0(names_models[i], '_GS'))
  
  AIC_Model_GS[i] <- AIC(mod_GS)
  DE_Model_GS[i] <- summary(mod_GS)$dev.expl
  
  rm(mod_GS)
  gc()
  
  print(i)
}

#build final table
GAM_evaluation <- data.frame(
  model = names_models,
  AIC_G = round(AIC_Model_G, 2),
  AIC_GS = round(AIC_Model_GS, 2),
  Delta_AIC = round(AIC_Model_GS - AIC_Model_G, 2),
  DE_G = round(DE_Model_G * 100, 2),
  DE_GS = round(DE_Model_GS * 100, 2),
  Delta_DE = round((DE_Model_GS - DE_Model_G) * 100, 2)
)

#rename labels
GAM_evaluation$model <- gsub('distEdge_250', 'distEdge',
                             GAM_evaluation$model)

#save table
setwd(wd_result)

write.csv(GAM_evaluation, 'Table_S2.csv', row.names = FALSE)
