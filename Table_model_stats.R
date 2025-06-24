#load libraries
library(mgcv); library(itsadug)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20250504_GAMs/Models'
wd_result <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/SI/Tables'

#load all models
setwd(wd_models)
models <- lapply(list.files(), readRDS)
names(models) <- list.files()

#get AIC for all models
AIC <- sapply(models, function(x){AIC(x)})

#get dev explanined for all models
DE <- sapply(models, function(x){summary(x)$dev.expl})

#combine results in a table

#get names of models without _G or _GS extensions
models_gen <- gsub('_GS', '', names(models))
models_gen <- gsub('_G', '', models_gen)
models_gen <- unique(models_gen)

#create a vector ordering the models as I want in the table
names_models <- c('minT_distEdge', 'meanT_distEdge', 'maxT_distEdge',
                  'minT_distEdge_250', 'meanT_distEdge_250', 'maxT_distEdge_250',
                  'minT_relPol', 'meanT_relPol', 'maxT_relPol',
                  'minT_absPol', 'meanT_absPol', 'maxT_absPol',
                  'minT_elev', 'meanT_elev', 'maxT_elev')


#create empty vectors to populate with results
AIC_Model_G <- numeric()
AIC_Model_GS <- numeric()
DE_Model_G <- numeric()
DE_Model_GS <- numeric()

#populate the vectors in a for loop
for(i in 1:length(names_models))
{
  AIC_Model_G[i] <- AIC[which(names(AIC) == paste0(names_models[i], '_G'))]
  AIC_Model_GS[i] <- AIC[which(names(AIC) == paste0(names_models[i], '_GS'))]
  DE_Model_G[i] <- DE[which(names(DE) == paste0(names_models[i], '_G'))] * 100
  DE_Model_GS[i] <- DE[which(names(DE) == paste0(names_models[i], '_GS'))] * 100
  
  print(i)
}

#join table
table <- data.frame(Models = names_models,
                    AIC_G = round(AIC_Model_G), AIC_GS = round(AIC_Model_GS),
                    Δ_AIC = round(AIC_Model_G - AIC_Model_GS),
                    DE_G = round(DE_Model_G, 2), DE_GS = round(DE_Model_GS, 2),
                    Δ_DE = round(DE_Model_GS - DE_Model_G, 2))

#calculate diff in deviance explained between distEdge models using all vars and those using only points up to 250km from the margins

diff_DE <- mean(table$DE_GS)

#calculate mean and SD of differences between G and GS
#(excluding elevation for now)

table_noElev <- table[-c(13:15),]

mean_Δ_AIC <- mean(table_noElev$Δ_AIC)
sd_Δ_AIC <- sd(table_noElev$Δ_AIC)

mean_D_DE <- mean(table_noElev$D_DE)
sd_D_DE <- sd(table_noElev$D_DE)

#save table
setwd(wd_result)
write.csv(table, 'GLM_evaluation.csv', row.names = F)


