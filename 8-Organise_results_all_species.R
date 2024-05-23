#load libraries
library(data.table)

#list wds
wd_pts_measure <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Point_and_range_measurements'
wd_res_shap <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Comparison'
wd_orders <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Species_lists'
wd_out <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All_species_analysis'

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
for(i in 1:length(sps_tables))
{
  sps_tables[[i]]$n_occ <- nrow(sps_tables[[i]])
}

#load shap results
setwd(wd_res_shap)
shap_results <- lapply(list.files(), read.csv)

names(shap_results) <- gsub('_.csv', '', list.files())  #name objects

#load table informing each species order
setwd(wd_orders)
sps_orders <- read.csv('Species_order.csv')

#include shap results in each sps_table
for(i in 1:length(sps_tables))
{
  #select each file with SHAP results for the species
  minT <- shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_minT')][[1]]
  meanT <- shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_meanT')][[1]]
  maxT <- shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_maxT')][[1]]
  
  minPPT <- shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_minPPT')][[1]]
  meanPPT <- shap_results[names(shap_results) == 
                          paste0(names(sps_tables)[i], '_meanPPT')][[1]]
  maxPPT <- shap_results[names(shap_results) == 
                         paste0(names(sps_tables)[i], '_maxPPT')][[1]]
  
  #select only presences
  minT2 <- minT[minT$Type == 'Presence',]
  meanT2 <- meanT[meanT$Type == 'Presence',]
  maxT2 <- maxT[maxT$Type == 'Presence',]
  
  minPPT2 <- minPPT[minPPT$Type == 'Presence',]
  meanPPT2 <- meanPPT[meanPPT$Type == 'Presence',]
  maxPPT2 <- maxPPT[maxPPT$Type == 'Presence',]

  #rename cols to indicate the 'control variables' 
  
  ## NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE ##
  
  # control_ is mean temperature for all precipitation analyses 
  # and mean precipitation for all temperature analyses
  
  # control_elev_ is elevation for all analyses
  
  ## NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE ##
  
  names(minT2)[names(minT2) == 'Mean_PPT_SHAP'] <- 'control_Min_T_SHAP'
  names(minT2)[names(minT2) == 'elev_SHAP'] <- 'control_elev_Min_T_SHAP'
  names(meanT2)[names(meanT2) == 'Mean_PPT_SHAP'] <- 'control_Mean_T_SHAP'
  names(meanT2)[names(meanT2) == 'elev_SHAP'] <- 'control_elev_Mean_T_SHAP'
  names(maxT2)[names(maxT2) == 'Mean_PPT_SHAP'] <- 'control_Max_T_SHAP'
  names(maxT2)[names(maxT2) == 'elev_SHAP'] <- 'control_elev_Max_T_SHAP'
  
  names(minPPT2)[names(minPPT2) == 'Mean_T_SHAP'] <- 'control_Min_PPT_SHAP'
  names(minPPT2)[names(minPPT2) == 'elev_SHAP'] <- 'control_elev_Min_PPT_SHAP'
  names(meanPPT2)[names(meanPPT2) == 'Mean_T_SHAP'] <- 'control_Mean_PPT_SHAP'
  names(meanPPT2)[names(meanPPT2) == 'elev_SHAP'] <- 'control_elev_Mean_PPT_SHAP'
  names(maxPPT2)[names(maxPPT2) == 'Mean_T_SHAP'] <- 'control_Max_PPT_SHAP'
  names(maxPPT2)[names(maxPPT2) == 'elev_SHAP'] <- 'control_elev_Max_PPT_SHAP'
  
  #select only columns of interest in each table
  minT3 <- minT2[c('decimalLongitude','decimalLatitude',
                   'Min_T_SHAP',
                   'control_Min_T_SHAP', 'control_elev_Min_T_SHAP')]
  meanT3 <- meanT2[c('decimalLongitude','decimalLatitude',
                     'Mean_T_SHAP',
                     'control_Mean_T_SHAP', 'control_elev_Mean_T_SHAP')]
  maxT3 <- maxT2[c('decimalLongitude','decimalLatitude',
                   'Max_T_SHAP',
                   'control_Max_T_SHAP', 'control_elev_Max_T_SHAP')]
  
  minPPT3 <- minPPT2[c('decimalLongitude','decimalLatitude',
                       'Min_PPT_SHAP',
                       'control_Min_PPT_SHAP', 'control_elev_Min_PPT_SHAP')]
  meanPPT3 <- meanPPT2[c('decimalLongitude','decimalLatitude',
                         'Mean_PPT_SHAP',
                         'control_Mean_PPT_SHAP', 'control_elev_Mean_PPT_SHAP')]
  maxPPT3 <- maxPPT2[c('decimalLongitude','decimalLatitude',
                       'Max_PPT_SHAP',
                       'control_Max_PPT_SHAP', 'control_elev_Max_PPT_SHAP')]
  
  #include order info into the species tables
  sps_tables[[i]]$order <- sps_orders$order[
    sps_orders$species == unique(sps_tables[[i]]$species)]

  #include values for each point into the species tables
  sps_tables[[i]] <- merge(sps_tables[[i]], minT3,
                           by = c('decimalLongitude','decimalLatitude'))
  sps_tables[[i]] <- merge(sps_tables[[i]], meanT3,
                           by = c('decimalLongitude','decimalLatitude'))
  sps_tables[[i]] <- merge(sps_tables[[i]], maxT3,
                           by = c('decimalLongitude','decimalLatitude'))
  
  sps_tables[[i]] <- merge(sps_tables[[i]], minPPT3,
                           by = c('decimalLongitude','decimalLatitude'))
  sps_tables[[i]] <- merge(sps_tables[[i]], meanPPT3,
                           by = c('decimalLongitude','decimalLatitude'))
  sps_tables[[i]] <- merge(sps_tables[[i]], maxPPT3,
                           by = c('decimalLongitude','decimalLatitude'))
  
  print(i)
}

#check which species do not have info on order
orders <- sapply(sps_tables, function(x){unique(x$order)})
missing <- which(is.na(orders))
missing 

#make a table with all species values
all_sps_table <- rbindlist(sps_tables, fill = T)

#save all species table
setwd(wd_out)
write.csv(all_sps_table, 'Results_all_sps.csv', row.names = F)

