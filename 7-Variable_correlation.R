#load packages
library(raster); library(sf); library(data.table)

#list wds
# wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_variables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Species_occ'
wd_res <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables/Variable_correlation'
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables'

#list species 
setwd(wd_occ)
sps_list <- gsub('_occ.csv', '', list.files())

#load all six BioCLim variables being used
setwd(wd_variables)
AnnualMeanTemperature <- raster('wc2.1_2.5m_bio_1.tif')
MaxTemperatureOfWarmestMonth <- raster('wc2.1_2.5m_bio_5.tif')
MinTemperatureOfColdestMonth <- raster('wc2.1_2.5m_bio_6.tif')
AnnualPrecipitation <- raster('wc2.1_2.5m_bio_12.tif')
PrecipitationOfWettestMonth <- raster('wc2.1_2.5m_bio_13.tif')
PrecipitationOfDriestMonth <- raster('wc2.1_2.5m_bio_14.tif')

#stack all variables 
var_stack <- stack(AnnualMeanTemperature, MaxTemperatureOfWarmestMonth,
                   MinTemperatureOfColdestMonth, AnnualPrecipitation,
                   PrecipitationOfWettestMonth, PrecipitationOfDriestMonth)


######## Correlation analysis for each species ######

#create empty object to save species names and n_occ
sps_name <- character()
n_occ <- numeric()

for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #save name in list
  sps_name[i] <- sps
  
  #load species occurrences
  setwd(wd_occ)
  sps_occ <- read.csv(paste0(sps, '_occ.csv'))
  
  #select only presence occurrences
  pr_sps <- sps_occ[which(sps_occ$Type == 'Presence'),]
  
  #save n_occ in list
  n_occ[i] <- nrow(pr_sps)
  
  #create sf object of the species occurrences
  pr_sps_sp <- st_as_sf(pr_sps,
                        coords = c('decimalLongitude', 'decimalLatitude'))
  
  #get values of variables at each poin
  values <- extract(var_stack, pr_sps_sp)
  
  #calculate correlation
  correl <- cor(values)
  
  #save results per species
  setwd(wd_res)
  write.csv(correl, paste0(sps,'.csv'), row.names = T)
  
  print(i)
}


######## Make table with correl values per each pair of variables ######

#make objects to save all info for the table
species_name <- character()

preds_min_T <- numeric()
preds_mean_T <- numeric()
preds_max_T <- numeric()

preds_min_PPT <- numeric()
preds_mean_PPT <- numeric()
preds_max_PPT <- numeric()

#for loop through results getting correlation of pairs of variable of interest
setwd(wd_res)

for(i in 1:length(list.files())){
  
  species_name[i] <- gsub('.csv', '', list.files()[i])
  
  #read correlation table for the species
  correl_table <- read.csv(list.files()[i])
  
  #fix problem of row.names being read as a column
  row.names(correl_table) <- correl_table[,1]
  correl_table <- correl_table[,-1]
  
  preds_min_T[i] <- correl_table['wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_6']
  preds_mean_T[i] <- correl_table['wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_1']
  preds_max_T[i] <- correl_table['wc2.1_2.5m_bio_12', 'wc2.1_2.5m_bio_5']

  preds_min_PPT[i] <- correl_table['wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_14']
  preds_mean_PPT[i] <- correl_table['wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_12']
  preds_max_PPT[i] <- correl_table['wc2.1_2.5m_bio_1', 'wc2.1_2.5m_bio_13']
  
  print(i)
}


#create data.frame with results
results_correl <- data.frame(Species = species_name,
                             n_occ = n_occ,
                             Cor_vars_minT = preds_min_T,
                             Cor_vars_meanT = preds_mean_T,
                             Cor_vars_maxT = preds_max_T,
                             Cor_vars_minPPT = preds_min_PPT,
                             Cor_vars_meanPPT = preds_mean_PPT,
                             Cor_vars_maxPPT = preds_max_PPT)

#save dataframe with results
setwd(wd_tables)
write.csv(results_correl, 'Correlation_variables.csv', row.names = F)


#include correlation analysis in final table 

#load dataframe with SHAP results
setwd(wd_tables)
info_species <- read.csv('Info_species.csv', stringsAsFactors = FALSE)

#join correlation analysis (including n_occ) to species info table
mistress_file <- merge(info_species, results_correl, by = 'Species',
                       all.x = TRUE)

#save mistress file
setwd(wd_tables)
write.csv(mistress_file, 'Mistress_file.csv', row.names = F)
