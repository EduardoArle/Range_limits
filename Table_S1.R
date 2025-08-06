#load libraries
library(data.table)

#list WDs
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Tables'
wd_all_sps_res <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_All_species_analysis'
wd_tab_SI <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/SI/Tables'

#load file with all results
setwd(wd_all_sps_res)
all_res <- read.csv('20241210_Results_all_sps.csv')

#select only species used in any model
all_res_models <- all_res[which(!is.na(all_res$avg_Min_T_SHAP) |
                                  !is.na(all_res$avg_Mean_T_SHAP) |
                                  !is.na(all_res$avg_Max_T_SHAP)),]

n_sps_model <- length(unique(all_res_models$species))

#get the total number of points that were used in any models 
n_pts_models <- nrow(all_res_models)

#get the total number of species for which any model was made
n_species_all <- length(unique(all_res$species))
n_species_all

#create a version of the file including only entries that were not too correlated
all_res_cor <- all_res_models[which(abs(all_res_models$Cor_vars_minT) <= 0.7 |
                                      abs(all_res_models$Cor_vars_meanT) <= 0.7 |
                                      abs(all_res_models$Cor_vars_maxT) <= 0.7),]

#n points that were used valid models (by var correl)
n_pts_val_models <- nrow(all_res_cor)
n_pts_val_models

#get the total number of species for which any model was made
n_species_model <- length(unique(all_res_cor$species))
n_species_model

#get info on biomes represented
length(unique(all_res_cor$biome))

#make a table containing species used and number of records
used_species_nOcc <- unique(as.data.table(all_res_cor),
                            by = c('species', 'nOcc'))

#select columns
species_info <- used_species_nOcc[,c(1,6,7)]

#save
setwd(wd_tab_SI)
write.csv(species_info, 'Table_S1.csv', row.names = F)
