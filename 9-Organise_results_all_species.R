#load libraries
library(data.table)

#list wds
wd_pts_metrics <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables/Occurrence_metrics'
wd_res_species <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SHAP_results'
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables'


#list species
setwd(wd_pts_metrics)
sps_list <- gsub('_point_range_metrics.csv', '', list.files())

#load files
setwd(wd_pts_metrics)
sps_metrics <- lapply(list.files(), read.csv)
names(sps_metrics) <- sps_list #name objects

#load shap results
setwd(wd_res_species)
shap_results <- lapply(list.files(pattern = '.csv$'), read.csv)
names(shap_results) <- gsub('.csv', '', list.files(pattern = '.csv$')) 

#load table with range metrics
setwd(wd_tables)
range_metrics <- read.csv('Selected_species.csv')

#make copy of sps tables
sps_tables <- sps_metrics

#include shap results in each sps_table
for(i in 1:length(sps_tables))
{
  #select each file with SHAP results for the species
  minT <- try(shap_results[names(shap_results) == 
                             paste0(names(sps_tables)[i], '_minT')][[1]],
              silent = TRUE)
  meanT <- try(shap_results[names(shap_results) == 
                              paste0(names(sps_tables)[i], '_meanT')][[1]],
               silent = TRUE)
  maxT <- try(shap_results[names(shap_results) == 
                             paste0(names(sps_tables)[i], '_maxT')][[1]],
              silent = TRUE)
  
  minPPT <- try(shap_results[names(shap_results) == 
                               paste0(names(sps_tables)[i], '_minPPT')][[1]],
                silent = TRUE)
  meanPPT <- try(shap_results[names(shap_results) == 
                                paste0(names(sps_tables)[i], '_meanPPT')][[1]],
                 silent = TRUE)
  maxPPT <- try(shap_results[names(shap_results) == 
                               paste0(names(sps_tables)[i], '_maxPPT')][[1]],
                silent = TRUE)
  
  #rename cols to indicate the 'control variables' 
  
  ##NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE##
  
  #control_ is mean temperature for all precipitation analyses 
  #and mean precipitation for all temperature analyses
  
  ##NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE##
  
  names(minT)[names(minT) == 'avg_Mean_PPT_SHAP'] <- 'control_Min_T_SHAP'
  names(meanT)[names(meanT) == 'avg_Mean_PPT_SHAP'] <- 'control_Mean_T_SHAP'
  names(maxT)[names(maxT) == 'avg_Mean_PPT_SHAP'] <- 'control_Max_T_SHAP'
  
  names(minPPT)[names(minPPT) == 'avg_Mean_T_SHAP'] <- 'control_Min_PPT_SHAP'
  names(meanPPT)[names(meanPPT) == 'avg_Mean_T_SHAP'] <- 'control_Mean_PPT_SHAP'
  names(maxPPT)[names(maxPPT) == 'avg_Mean_T_SHAP'] <- 'control_Max_PPT_SHAP'
  
  #select only columns of interest in each table
  cols_minT  <- c('decimalLongitude','decimalLatitude',
                  'avg_Min_T_SHAP','control_Min_T_SHAP')
  cols_meanT <- c('decimalLongitude','decimalLatitude',
                  'avg_Mean_T_SHAP','control_Mean_T_SHAP')
  cols_maxT  <- c('decimalLongitude','decimalLatitude',
                  'avg_Max_T_SHAP','control_Max_T_SHAP')
  
  cols_minPPT  <- c('decimalLongitude','decimalLatitude',
                    'avg_Min_PPT_SHAP','control_Min_PPT_SHAP')
  cols_meanPPT <- c('decimalLongitude','decimalLatitude',
                    'avg_Mean_PPT_SHAP','control_Mean_PPT_SHAP')
  cols_maxPPT  <- c('decimalLongitude','decimalLatitude',
                    'avg_Max_PPT_SHAP','control_Max_PPT_SHAP')
  
  minT2 <- try(minT[, intersect(cols_minT, names(minT)), drop = FALSE],
               silent = TRUE)
  if(inherits(minT2, "try-error")) minT2 <- NULL
  
  meanT2 <- try(meanT[, intersect(cols_meanT, names(meanT)), drop = FALSE],
                silent = TRUE)
  if(inherits(meanT2, "try-error")) meanT2 <- NULL
  
  maxT2 <- try(maxT[, intersect(cols_maxT, names(maxT)), drop = FALSE],
               silent = TRUE)
  if(inherits(maxT2, "try-error")) maxT2 <- NULL
  
  
  minPPT2 <- try(minPPT[, intersect(cols_minPPT, names(minPPT)), drop = FALSE],
                 silent = TRUE)
  if(inherits(minPPT2, "try-error")) minPPT2 <- NULL
  
  meanPPT2 <- try(meanPPT[, intersect(cols_meanPPT, names(meanPPT)), drop = FALSE],
                  silent = TRUE)
  if(inherits(meanPPT2, "try-error")) meanPPT2 <- NULL
  
  maxPPT2 <- try(maxPPT[, intersect(cols_maxPPT, names(maxPPT)), drop = FALSE],
                 silent = TRUE)
  if(inherits(maxPPT2, "try-error")) maxPPT2 <- NULL
  
  #include values for each point into the species tables
  
  if(is.data.frame(minT2) &&
     all(c('decimalLongitude','decimalLatitude') %in% names(minT2)) &&
     ncol(minT2) > 2){
    sps_tables[[i]] <- merge(sps_tables[[i]], minT2,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  if(is.data.frame(meanT2) &&
     all(c('decimalLongitude','decimalLatitude') %in% names(meanT2)) &&
     ncol(meanT2) > 2){
    sps_tables[[i]] <- merge(sps_tables[[i]], meanT2,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  if(is.data.frame(maxT2) &&
     all(c('decimalLongitude','decimalLatitude') %in% names(maxT2)) &&
     ncol(maxT2) > 2){
    sps_tables[[i]] <- merge(sps_tables[[i]], maxT2,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  if(is.data.frame(minPPT2) &&
     all(c('decimalLongitude','decimalLatitude') %in% names(minPPT2)) &&
     ncol(minPPT2) > 2){
    sps_tables[[i]] <- merge(sps_tables[[i]], minPPT2,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  if(is.data.frame(meanPPT2) &&
     all(c('decimalLongitude','decimalLatitude') %in% names(meanPPT2)) &&
     ncol(meanPPT2) > 2){
    sps_tables[[i]] <- merge(sps_tables[[i]], meanPPT2,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  if(is.data.frame(maxPPT2) &&
     all(c('decimalLongitude','decimalLatitude') %in% names(maxPPT2)) &&
     ncol(maxPPT2) > 2){
    sps_tables[[i]] <- merge(sps_tables[[i]], maxPPT2,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  ## Include range geographical details for each species
  
  #species name in this table
  sps <- sps_tables[[i]]$sps[1]
  
  #row in range_metrics
  idx <- match(sps, range_metrics$species)
  
  #attach species-level info
  sps_tables[[i]]$dom_ratio <- range_metrics$dom_ratio[idx]
  sps_tables[[i]]$perc_ocean <- range_metrics$perc_ocean[idx]
  sps_tables[[i]]$N_ocean <- range_metrics$N_ocean[idx]
  sps_tables[[i]]$S_ocean <- range_metrics$S_ocean[idx]
  sps_tables[[i]]$E_ocean <- range_metrics$E_ocean[idx]
  sps_tables[[i]]$W_ocean <- range_metrics$W_ocean[idx]
  
  print(i)
}


#exclude empty items on list that had no models
no_models <- sapply(sps_tables, ncol)

#eliminate tables that produced no models
sps_tables2 <- sps_tables[no_models > min(no_models)]

#check which species do not have info on order
orders <- sapply(sps_tables2, function(x){unique(x$order)})
missing <- which(is.na(orders))
missing 

#make a table with all species values
all_sps_table <- rbindlist(sps_tables2, fill = T)

#load correlation table
setwd(wd_tables)
correl <- read.csv('Correlation_variables.csv')

#delete nOcc col
correl <- correl[,-which(names(correl) == 'n_occ')]

#harmonise col names in both tables
names(correl)[1] <- 'sps'

#eliminate species that are not in the all_sps_table
correl2 <- correl[which(correl$sps %in% unique(all_sps_table$sps)),]

all_sps_table2 <- merge(all_sps_table, correl2,
                        by = 'sps', all.x = T)

all_sps_table2 <- as.data.frame(all_sps_table2)

#save all species table
setwd(wd_tables)
write.csv(all_sps_table2, '20260126_Results_all_sps.csv', row.names = F)

