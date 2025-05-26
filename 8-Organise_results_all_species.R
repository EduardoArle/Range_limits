#load libraries
library(data.table)

#list wds
wd_pts_measure <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_Point_and_range_measurements'
wd_res_shap <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20250421_Comparison'
wd_orders <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Species_lists'
wd_out <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20250504_All_species_analysis'
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Tables'

#list species
setwd(wd_pts_measure)
sps_list <- gsub('_point_range_metrics.csv', '', list.files())

#load files
setwd(wd_pts_measure)
sps_tables <- lapply(list.files(), read.csv)

names(sps_tables) <- sps_list #name objects

#fix name of _Dasypus yepesi_ that for some reason is listed with the synonym 
# _Dasypus mazzai_ in the table
sps_tables[[60]]$species <- names(sps_tables)[60]

#same for _Artibeus anderseni_ listed with the synonym _Dermanura anderseni_ 
sps_tables[[61]]$species <- names(sps_tables)[61]

#same for _Eudorcas albonotata_ listed with the synonym _Eudorcas rufifrons_ 
sps_tables[[87]]$species <- names(sps_tables)[87]

#same for _Gerbillus stigmonyx_ listed with the synonym _Dipodillus stigmonyx_ 
sps_tables[[105]]$species <- names(sps_tables)[105]

#same for _Handleyomys chapmani_ listed with the synonym _Oryzomys chapmani_ 
sps_tables[[110]]$species <- names(sps_tables)[110]

#same for _Micoureus constantiae_ listed with the synonym _Marmosa constantiae_
sps_tables[[170]]$species <- names(sps_tables)[170]

#same for _Micoureus phaeus_ listed with the synonym _Marmosa phaeus_
sps_tables[[171]]$species <- names(sps_tables)[171]

#same for _Neoromicia matroka_ listed with the synonym _Laephotis matroka_
sps_tables[[234]]$species <- names(sps_tables)[234]

#same for _Neotamias alpinus_ listed with the synonym _Tamias alpinus_
sps_tables[[235]]$species <- names(sps_tables)[235]

#same for _Neotamias bulleri_ listed with the synonym _Tamias bulleri_
sps_tables[[236]]$species <- names(sps_tables)[236]

#same for _Neotamias cinereicollis_ listed with the synonym _Tamias cinereicollis_
sps_tables[[237]]$species <- names(sps_tables)[237]

#same for _Neotamias quadrivittatus_ listed with the synonym _Tamias quadrivittatus_
sps_tables[[238]]$species <- names(sps_tables)[238]

#same for _Otomys sloggetti_ listed with the synonym _Myotomys sloggetti_
sps_tables[[285]]$species <- names(sps_tables)[285]

#same for _Otomys sloggetti_ listed with the synonym _Myotomys sloggetti_
sps_tables[[285]]$species <- names(sps_tables)[285]

#same for _Phaiomys leucurus_ listed with the synonym _Neodon leucurus_
sps_tables[[304]]$species <- names(sps_tables)[304]

#same for _Rhabdomys intermedius_ listed with the synonym _Rhabdomys pumilio_
sps_tables[[382]]$species <- names(sps_tables)[382]

#same for _Toromys rhipidurus_ listed with the synonym _Makalata rhipidura_
sps_tables[[490]]$species <- names(sps_tables)[490]

# sps_tables[[i]]$species
# names(sps_tables)[i]

#create a column with n occurrences for each species
# for(i in 1:length(sps_tables))
# {
#   sps_tables[[i]]$n_occ <- nrow(sps_tables[[i]])
# }

#load shap results
setwd(wd_res_shap)
shap_results <- lapply(list.files(pattern = '.csv$'), read.csv)

#name objects
names(shap_results) <- gsub('.csv', '', list.files(pattern = '.csv$')) 

#load table informing each species order
setwd(wd_orders)
sps_orders <- read.csv('Species_order.csv')








#include shap results in each sps_table
for(i in 1:length(sps_tables))
{
  #select each file with SHAP results for the species
  minT <- try(shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_minT')][[1]],
              silent = T)
  meanT <- try(shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_meanT')][[1]],
               silent = T)
  maxT <- try(shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_maxT')][[1]],
              silent = T)
  
  minPPT <- try(shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_minPPT')][[1]],
                silent = T)
  meanPPT <- try(shap_results[names(shap_results) == 
                          paste0(names(sps_tables)[i], '_meanPPT')][[1]],
                 silent = T)
  maxPPT <- try(shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_maxPPT')][[1]],
                silent = T)
  
  #select only presences
  minT2 <- try(minT[minT$Type == 'Presence',], silent = T)
  meanT2 <- try(meanT[meanT$Type == 'Presence',], silent = T)
  maxT2 <- try(maxT[maxT$Type == 'Presence',], silent = T)
  
  minPPT2 <- try(minPPT[minPPT$Type == 'Presence',], silent = T)
  meanPPT2 <- try(meanPPT[meanPPT$Type == 'Presence',], silent = T)
  maxPPT2 <- try(maxPPT[maxPPT$Type == 'Presence',], silent = T)

  #rename cols to indicate the 'control variables' 
  
  ## NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE ##
  
  # control_ is mean temperature for all precipitation analyses 
  # and mean precipitation for all temperature analyses
  
  # control_elev_ is elevation for all analyses
  
  ## NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE ##
  
  names(minT2)[names(minT2) == 'avg_Mean_PPT_SHAP'] <- 'control_Min_T_SHAP'
  names(meanT2)[names(meanT2) == 'avg_Mean_PPT_SHAP'] <- 'control_Mean_T_SHAP'
  names(maxT2)[names(maxT2) == 'avg_Mean_PPT_SHAP'] <- 'control_Max_T_SHAP'

  names(minPPT2)[names(minPPT2) == 'avg_Mean_T_SHAP'] <- 'control_Min_PPT_SHAP'
  names(meanPPT2)[names(meanPPT2) == 'avg_Mean_T_SHAP'] <- 'control_Mean_PPT_SHAP'
  names(maxPPT2)[names(maxPPT2) == 'avg_Mean_T_SHAP'] <- 'control_Max_PPT_SHAP'

  #select only columns of interest in each table
  minT3 <- minT2[c('decimalLongitude', 'decimalLatitude',
                   'avg_Min_T_SHAP', 'control_Min_T_SHAP')]
  meanT3 <- meanT2[c('decimalLongitude', 'decimalLatitude',
                     'avg_Mean_T_SHAP', 'control_Mean_T_SHAP')]
  maxT3 <- maxT2[c('decimalLongitude', 'decimalLatitude',
                   'avg_Max_T_SHAP', 'control_Max_T_SHAP')]
  
  minPPT3 <- minPPT2[c('decimalLongitude', 'decimalLatitude',
                       'avg_Min_PPT_SHAP', 'control_Min_PPT_SHAP')]
  meanPPT3 <- meanPPT2[c('decimalLongitude', 'decimalLatitude',
                         'avg_Mean_PPT_SHAP', 'control_Mean_PPT_SHAP')]
  maxPPT3 <- maxPPT2[c('decimalLongitude','decimalLatitude',
                       'avg_Max_PPT_SHAP','control_Max_PPT_SHAP')]
  
  #include order info into the species tables
  sps_tables[[i]]$order <- sps_orders$order[
    sps_orders$species == unique(sps_tables[[i]]$species)]

  #include values for each point into the species tables
  if(class(minT3) == 'data.frame'){
    sps_tables[[i]] <- merge(sps_tables[[i]], minT3,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  if(class(meanT3) == 'data.frame'){
    sps_tables[[i]] <- merge(sps_tables[[i]], meanT3,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  if(class(maxT3) == 'data.frame'){
    sps_tables[[i]] <- merge(sps_tables[[i]], maxT3,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  if(class(minPPT3) == 'data.frame'){
    sps_tables[[i]] <- merge(sps_tables[[i]], minPPT3,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  if(class(meanPPT3) == 'data.frame'){
    sps_tables[[i]] <- merge(sps_tables[[i]], meanPPT3,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  if(class(maxPPT3) == 'data.frame'){
    sps_tables[[i]] <- merge(sps_tables[[i]], maxPPT3,
                             by = c('decimalLongitude','decimalLatitude'))
  }
  
  
  print(i)
}

#exclude empty items on list that had no models
no_models <- sapply(sps_tables, ncol)
sps_tables2 <- sps_tables[no_models > 17] #17 is the ncol of tables that produced no models

#exclude Isothrix orinoci (synonym)
sps_tables2 <- sps_tables2[-which(names(sps_tables2) == 'Isothrix orinoci')]

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
names(correl)[1] <- 'species'

#eliminate species that are not in the all_sps_table
correl2 <- correl[which(correl$species %in% unique(all_sps_table$species)),]

all_sps_table2 <- merge(all_sps_table, correl2,
                        by = 'species', all.x = T)

all_sps_table2 <- as.data.frame(all_sps_table2)

#save all species table
setwd(wd_out)
write.csv(all_sps_table2, '20250504_Results_all_sps.csv', row.names = F)

