####  This script combines the SHAP results for each species with the info
####  about range characteristics and position of points within the ranges

#load libraries

#list wds
wd_res_shap <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Comparison'
wd_orders <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Species_lists'
wd_pts_measure <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Point_and_range_measurements'
wd_res <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240603/Each_species'
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Results_20240603'


#list species 
setwd(wd_pts_measure)
sps_list <- gsub('_point_range_metrics.csv', '', list.files())

#load all species measurements
setwd(wd_pts_measure)
sps_measures <- lapply(list.files(), read.csv)
names(sps_measures) <- sps_list

#load table informing each species order
setwd(wd_orders)
sps_orders <- read.csv('Species_order.csv')


##############################################
############## EACH SPECIES TEMP #############
##############################################


## Load results for T variables considering all species
setwd(wd_res_shap)

#minT
sps_minT <- lapply(sps_list, function(x){read.csv(paste0(x,'_minT_.csv'))})
names(sps_minT) <- sps_list

#meanT
sps_meanT <- lapply(sps_list, function(x){read.csv(paste0(x,'_meanT_.csv'))})
names(sps_meanT) <- sps_list

#maxT
sps_maxT <- lapply(sps_list, function(x){read.csv(paste0(x,'_maxT_.csv'))})
names(sps_maxT) <- sps_list


## Select only rows representing presence
sps_minT_pr <- lapply(sps_minT, function(x){x[x$Occurrence == 1,]}) #minT
sps_meanT_pr <- lapply(sps_meanT, function(x){x[x$Occurrence == 1,]}) #meanT
sps_maxT_pr <- lapply(sps_maxT, function(x){x[x$Occurrence == 1,]}) #maxT


## Concatenate point and range info to SHAP results

#minT
sps_minT_SHAP_info <- list()
for(j in 1:length(sps_minT_pr))
{
  sps_minT_SHAP_info[[j]] <- cbind(sps_minT_pr[[j]], 
                                   sps_measures[[j]][,c('key','datasetKey',
                                  'rangeSize','roundness','centralness',     
                                  'distEdge','absPolarwardness',
                                  'relPolarwardness','elevation','biome',
                                  'bodyMass')])
}

#meanT
sps_meanT_SHAP_info <- list()
for(j in 1:length(sps_meanT_pr))
{
  sps_meanT_SHAP_info[[j]] <- cbind(sps_meanT_pr[[j]], 
                                    sps_measures[[j]][,c('key','datasetKey',
                                                    'rangeSize','roundness',
                                                    'centralness',     
                                                    'distEdge',
                                                    'absPolarwardness',
                                                    'relPolarwardness',
                                                    'elevation','biome',
                                                    'bodyMass')])
}

#maxT
sps_maxT_SHAP_info <- list()
for(j in 1:length(sps_maxT_pr))
{
  sps_maxT_SHAP_info[[j]] <- cbind(sps_maxT_pr[[j]], 
                                   sps_measures[[j]][,c('key','datasetKey',
                                                    'rangeSize','roundness',
                                                    'centralness',     
                                                    'distEdge',
                                                    'absPolarwardness',
                                                    'relPolarwardness',
                                                    'elevation','biome',
                                                    'bodyMass')])
}

#create vectors to store the results for each species
species <- character()
order <- character()
bodymass <- numeric()
n_biomes <- numeric()
min_elev <- numeric()
mean_elev <- numeric()
max_elev <- numeric()
n_records <- numeric()
hemisphere <- character()
minTR2_RP <- numeric()
minTslope_RP <- numeric()
minTSE_RP <- numeric()
minTp_value_RP <- numeric()
meanTR2_RP <- numeric()
meanTslope_RP <- numeric()
meanTSE_RP <- numeric()
meanTp_value_RP <- numeric()
maxTR2_RP <- numeric()
maxTslope_RP <- numeric()
maxTSE_RP <- numeric()
maxTp_value_RP <- numeric()
minTR2_C <- numeric()
minTslope_C <- numeric()
minTSE_C <- numeric()
minTp_value_C <- numeric()
meanTR2_C <- numeric()
meanTslope_C <- numeric()
meanTSE_C <- numeric()
meanTp_value_C <- numeric()
maxTR2_C <- numeric()
maxTslope_C <- numeric()
maxTSE_C <- numeric()
maxTp_value_C <- numeric()
rangeSize <- numeric()
rangeRoundness <- numeric()
maxLat <- numeric()
minLat <- numeric()

 
## Loop through species
for(j in 1:length(sps_list))
{
  #skip sps with only one point
  if(nrow(sps_minT_SHAP_info[[j]]) > 1){
    #list spp
    species[length(species) + 1] <- sps_list[j]
    
    #list order
    order[length(order) + 1] <-
      sps_orders$order[sps_orders$species == sps_list[j]]
    
    #list bodymass
    bodymass[length(bodymass) + 1] <- unique(sps_minT_SHAP_info[[j]]$bodyMass)
    
    #list number of biomes
    n_biomes[length(n_biomes) + 1] <-
      length(unique(sps_minT_SHAP_info[[j]]$biome))
    
    #list min_elev
    min_elev[length(min_elev) + 1] <- min(sps_minT_SHAP_info[[j]]$elevation)
    
    #list mean_elev
    mean_elev[length(mean_elev) + 1] <- mean(sps_minT_SHAP_info[[j]]$elevation)
    
    #list max_elev
    max_elev[length(max_elev) + 1] <- max(sps_minT_SHAP_info[[j]]$elevation)
    
    #check in which hemisphere the species is
    if(length(unique(sps_minT_SHAP_info[[j]]$decimalLatitude < 0)) == 1){
      if(unique(sps_minT_SHAP_info[[j]]$decimalLatitude < 0)){
        hemisphere[length(hemisphere) + 1] <- 'Southern'
      }
      if(unique(sps_minT_SHAP_info[[j]]$decimalLatitude > 0)){
        hemisphere[length(hemisphere) + 1] <- 'Northern'
    }
    }
    
    if(length(unique(sps_minT_SHAP_info[[j]]$decimalLatitude < 0)) == 2){
      hemisphere[length(hemisphere) + 1] <- 'Both'
    }
    
    ## Create a folder for the species results
    wd_sps <- paste0(wd_res, '/', gsub(' ', '_', sps_list[j]))
    dir.create(wd_sps)
    
    ## Get SHAP values of the target variable for each point
    
    #minT
    cont_minT <- sps_minT_SHAP_info[[j]]$Min_T_SHAP
    
    #meanT
    cont_meanT <- sps_meanT_SHAP_info[[j]]$Mean_T_SHAP
    
    #maxT
    cont_maxT <- sps_maxT_SHAP_info[[j]]$Max_T_SHAP
    
    ## Save plot temperature variables contribution per relative polewardness
    
    #only for species not crossing the equator
    
    if(hemisphere[length(hemisphere)] != 'Both'){
      setwd(wd_sps)
      pdf(file = paste0(gsub(' ', '_', 
                             sps_list[j]),'_Temperature_RelPolewardness.pdf'))
      
      #set parametres for plotting
      par(mar = c(5,5,5,5))
      
      #get y and x lims
      ylim <- range(c(cont_minT, cont_meanT, cont_maxT))
      xlim <- range(sps_minT_SHAP_info[[j]]$relPolarwardness) #no matter which
      
      #minT
      plot(sps_minT_SHAP_info[[j]]$relPolarwardness, cont_minT, 
           pch = 19, cex = 0.4, col = '#0000FF15',
           ylab = 'Temperature contribution',
           xlab = 'Relative polewardness',
           ylim = c(ylim[1], ylim[2]))
      
      #fit linear model
      lin_mod_minT <- lm(cont_minT ~ sps_minT_SHAP_info[[j]]$relPolarwardness)
      abline(lin_mod_minT, col = '#0000FF', lwd = 2)
      
      #get and plot R^2
      summa <- summary(lin_mod_minT)
      r2 <- summa$r.squared
      text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
           ylim[1], 
           expression(bold(paste('R'^2*' = '), sep='')),
           col = '#0000FF', pos = 4)
      text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
           ylim[1] - ((ylim[2] - ylim[1]) / 300), 
           round(r2,3), col = '#0000FF', font = 2, pos = 4)
      
      #prepare model results and populate vectors
      summa <- as.data.frame(summa$coefficients)
      
      minTR2_RP[length(minTR2_RP) + 1] <- r2
      minTslope_RP[length(minTslope_RP) + 1] <- summa$Estimate[2]
      minTSE_RP[length(minTSE_RP) + 1] <- summa$`Std. Error`[2]
      minTp_value_RP[length(minTp_value_RP) + 1] <- summa$`Pr(>|t|)`[2]
      
      #meanT
      points(sps_meanT_SHAP_info[[j]]$relPolarwardness, cont_meanT, 
             pch = 19, cex = 0.4, col = '#80008015')
      
      #fit linear model
      lin_mod_meanT <- lm(cont_meanT ~ sps_meanT_SHAP_info[[j]]$relPolarwardness)
      abline(lin_mod_meanT, col = '#800080', lwd = 2)
      
      #get and plot R^2
      summa <- summary(lin_mod_meanT)
      r2 <- summa$r.squared
      text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
           ylim[1] + ((ylim[2] - ylim[1]) / 22), 
           expression(bold(paste('R'^2*' = '), sep='')), col = '#800080', pos = 4)
      text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
           ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
           round(r2,3), col = '#800080', font = 2, pos = 4)
      
      #prepare model results and populate vectors
      summa <- as.data.frame(summa$coefficients)
      
      meanTR2_RP[length(meanTR2_RP) + 1] <- r2
      meanTslope_RP[length(meanTslope_RP) + 1] <- summa$Estimate[2]
      meanTSE_RP[length(meanTSE_RP) + 1] <- summa$`Std. Error`[2]
      meanTp_value_RP[length(meanTp_value_RP) + 1] <- summa$`Pr(>|t|)`[2]
      
      #maxT
      points(sps_maxT_SHAP_info[[j]]$relPolarwardness, cont_maxT, 
             pch = 19, cex = 0.4, col = '#FF000015')
      
      #fit linear model
      lin_mod_maxT <- lm(cont_maxT ~ sps_maxT_SHAP_info[[j]]$relPolarwardness)
      abline(lin_mod_maxT, col = '#FF0000', lwd = 2)
      
      #get and plot R^2
      summa <- summary(lin_mod_maxT)
      r2 <- summa$r.squared
      text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
           ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
           expression(bold(paste('R'^2*' = '), sep='')), col = '#FF0000', pos = 4)
      text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
           ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
           round(r2,3), col = '#FF0000', font = 2, pos = 4)
      
      #prepare model results and populate vectors
      summa <- as.data.frame(summa$coefficients)
      
      maxTR2_RP[length(maxTR2_RP) + 1] <- r2
      maxTslope_RP[length(maxTslope_RP) + 1]  <- summa$Estimate[2]
      maxTSE_RP[length(maxTSE_RP) + 1] <- summa$`Std. Error`[2]
      maxTp_value_RP[length(maxTp_value_RP) + 1] <- summa$`Pr(>|t|)`[2]
      
      
      #plot species range size and roundness
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5), 
           ylim[2] - ((ylim[2] - ylim[1]) / 25), 
           'Range', font = 2, cex = 0.8, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 2 * ((ylim[2] - ylim[1]) / 25), 
           paste0('Size  ', unique(round(sps_maxT_SHAP_info[[j]]$rangeSize))), 
           cex = 0.7, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 2.8 * ((ylim[2] - ylim[1]) / 25),
           paste0('Roundness  ',
                  unique(round(sps_maxT_SHAP_info[[j]]$roundness, 2))),
           cex = 0.7, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 3.6 * ((ylim[2] - ylim[1]) / 25),
           paste0('Max lat  ',
                  round(max(sps_maxT_SHAP_info[[j]]$decimalLatitude), 2)),
           cex = 0.7, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 4.4 * ((ylim[2] - ylim[1]) / 25),
           paste0('Min lat  ',
                  round(min(sps_minT_SHAP_info[[j]]$decimalLatitude), 2)),
           cex = 0.7, pos = 4)
      
      #plot species name
      text(xlim[1] + ((xlim[2] - xlim[1]) / 20),
           ylim[2] - ((ylim[2] - ylim[1]) / 25),
           sps_list[j], font = 4, pos = 4)
      
      dev.off()
    }else{
      
      minTR2_RP[length(minTR2_RP) + 1] <- NA
      minTslope_RP[length(minTslope_RP) + 1] <- NA
      minTSE_RP[length(minTSE_RP) + 1] <- NA
      minTp_value_RP[length(minTp_value_RP) + 1] <- NA
      
      meanTR2_RP[length(meanTR2_RP) + 1] <- NA
      meanTslope_RP[length(meanTslope_RP) + 1] <- NA
      meanTSE_RP[length(meanTSE_RP) + 1] <- NA
      meanTp_value_RP[length(meanTp_value_RP) + 1] <- NA
      
      maxTR2_RP[length(maxTR2_RP) + 1] <- NA
      maxTslope_RP[length(maxTslope_RP) + 1] <- NA
      maxTSE_RP[length(maxTSE_RP) + 1] <- NA
      maxTp_value_RP[length(maxTp_value_RP) + 1] <- NA
    }
    
    ## Save plot temperature variables contribution per centralness
    setwd(wd_sps)
    pdf(file = paste0(gsub(' ', '_',
                           sps_list[j]),'_Temperature_DistEdge.pdf'))
    
    #set parametres for plotting
    par(mar = c(5,5,5,5))
    
    #get y and x lims
    ylim <- range(c(cont_minT, cont_meanT, cont_maxT))
    xlim <- range(sps_minT_SHAP_info[[j]]$distEdge) #does not matter which I use
    
    #minT
    plot(sps_minT_SHAP_info[[j]]$distEdge, cont_minT, 
         pch = 19, cex = 0.4, col = '#0000FF15',
         ylab = 'Temperature contribution',
         xlab = 'Distance from edge (km)',
         ylim = c(ylim[1], ylim[2]))
    
    #fit linear model
    lin_mod_minT <- lm(cont_minT ~ sps_minT_SHAP_info[[j]]$distEdge)
    abline(lin_mod_minT, col = '#0000FF', lwd = 2)
    
    #get and plot R^2
    summa <- summary(lin_mod_minT)
    r2 <- summa$r.squared
    text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
         ylim[1], 
         expression(bold(paste('R'^2*' = '), sep='')), col = '#0000FF', pos = 4)
    text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
         ylim[1] - ((ylim[2] - ylim[1]) / 300), 
         round(r2,3), col = '#0000FF', font = 2, pos = 4)
    
    #prepare model results and populate vectors
    summa <- as.data.frame(summa$coefficients)
    
    minTR2_C[length(minTR2_C) + 1] <- r2
    minTslope_C[length(minTslope_C) + 1]  <- summa$Estimate[2]
    minTSE_C[length(minTSE_C) + 1] <- summa$`Std. Error`[2]
    minTp_value_C[length(minTp_value_C) + 1] <- summa$`Pr(>|t|)`[2]
  
    
    #meanT
    points(sps_meanT_SHAP_info[[j]]$distEdge, cont_meanT, 
           pch = 19, cex = 0.4, col = '#80008010')
    
    #fit linear model
    lin_mod_meanT <- lm(cont_meanT ~ sps_meanT_SHAP_info[[j]]$distEdge)
    abline(lin_mod_meanT, col = '#800080', lwd = 2)
    
    #get and plot R^2
    summa <- summary(lin_mod_meanT)
    r2 <- summa$r.squared
    text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
         ylim[1] + ((ylim[2] - ylim[1]) / 22), 
         expression(bold(paste('R'^2*' = '), sep='')), col = '#800080', pos = 4)
    text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
         ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
         round(r2,3), col = '#800080', font = 2, pos = 4)
    
    #prepare model results and populate vectors
    summa <- as.data.frame(summa$coefficients)
    
    meanTR2_C[length(meanTR2_C) + 1] <- r2
    meanTslope_C[length(meanTslope_C) + 1]  <- summa$Estimate[2]
    meanTSE_C[length(meanTSE_C) + 1] <- summa$`Std. Error`[2]
    meanTp_value_C[length(meanTp_value_C) + 1] <- summa$`Pr(>|t|)`[2]

    
    #maxT
    points(sps_maxT_SHAP_info[[j]]$distEdge, cont_maxT, 
           pch = 19, cex = 0.4, col = '#FF000010')
    
    #fit linear model
    lin_mod_maxT <- lm(cont_maxT ~ sps_maxT_SHAP_info[[j]]$distEdge)
    abline(lin_mod_maxT, col = '#FF0000', lwd = 2)
    
    #get and plot R^2
    summa <- summary(lin_mod_maxT)
    r2 <- summa$r.squared
    text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
         ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
         expression(bold(paste('R'^2*' = '), sep='')), col = '#FF0000', pos = 4)
    text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
         ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
         round(r2,3), col = '#FF0000', font = 2, pos = 4)
    
    #prepare model results and populate vectors
    summa <- as.data.frame(summa$coefficients)
    
    maxTR2_C[length(maxTR2_C) + 1] <- r2
    maxTslope_C[length(maxTslope_C) + 1]  <- summa$Estimate[2]
    maxTSE_C[length(maxTSE_C) + 1] <- summa$`Std. Error`[2]
    maxTp_value_C[length(maxTp_value_C) + 1] <- summa$`Pr(>|t|)`[2]
    
    #plot species range size and roundness
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5), 
         ylim[2] - ((ylim[2] - ylim[1]) / 25), 
         'Range', font = 2, cex = 0.8, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 2 * ((ylim[2] - ylim[1]) / 25), 
         paste0('Size  ', unique(round(sps_maxT_SHAP_info[[j]]$rangeSize))), 
         cex = 0.7, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 2.8 * ((ylim[2] - ylim[1]) / 25),
         paste0('Roundness  ',
                unique(round(sps_maxT_SHAP_info[[j]]$roundness, 2))),
         cex = 0.7, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 3.6 * ((ylim[2] - ylim[1]) / 25),
         paste0('Max lat  ',
                round(max(sps_maxT_SHAP_info[[j]]$decimalLatitude), 2)),
         cex = 0.7, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 4.4 * ((ylim[2] - ylim[1]) / 25),
         paste0('Min lat  ',
                round(min(sps_minT_SHAP_info[[j]]$decimalLatitude), 2)),
         cex = 0.7, pos = 4)
    
    #plot species name
    text(xlim[1] + ((xlim[2] - xlim[1]) / 20),
         ylim[2] - ((ylim[2] - ylim[1]) / 25),
         sps_list[j], font = 4, pos = 4)
    
    dev.off()
    
    rangeSize[length(rangeSize) + 1] <- unique(sps_maxT_SHAP_info[[j]]$rangeSize)
    rangeRoundness[length(rangeRoundness) + 1] <- 
      unique(sps_maxT_SHAP_info[[j]]$roundness) 
    maxLat[length(maxLat) + 1] <- max(abs(sps_maxT_SHAP_info[[j]]$decimalLatitude))
    minLat[length(minLat) + 1] <- min(abs(sps_maxT_SHAP_info[[j]]$decimalLatitude))
    n_records[length(n_records) + 1] <- length(cont_minT)
    
  }
  
  print(j)
}

#create a dataframe with the results for all species
temp_all_pts <- data.frame(species = species,
                           order = order,
                           bodymass = bodymass,
                           n_biomes = n_biomes,
                           min_elev = min_elev,
                           mean_elev = mean_elev,
                           max_elev = max_elev,
                           n_records = n_records,
                           hemisphere = hemisphere,
                           rangeSize = rangeSize,
                           rangeRoundness = rangeRoundness,
                           maxLat = maxLat,
                           minLat = minLat,
                           minTR2_RP = minTR2_RP,
                           minTslope_RP = minTslope_RP,
                           minTSE_RP = minTSE_RP,
                           minTp_value_RP = minTp_value_RP,
                           meanTR2_RP = meanTR2_RP, 
                           meanTslope_RP = meanTslope_RP,
                           meanTSE_RP = meanTSE_RP,
                           meanTp_value_RP = meanTp_value_RP,
                           maxTR2_RP = maxTR2_RP,
                           maxTslope_RP = maxTslope_RP,
                           maxTSE_RP = maxTSE_RP,
                           maxTp_value_RP = maxTp_value_RP,
                           minTR2_C = minTR2_C,
                           minTslope_C = minTslope_C, 
                           minTSE_C = minTSE_C,
                           minTp_value_C = minTp_value_C,
                           meanTR2_C = meanTR2_C,
                           meanTslope_C = meanTslope_C,
                           meanTSE_C = meanTSE_C,
                           meanTp_value_C = meanTp_value_C,
                           maxTR2_C = maxTR2_C,
                           maxTslope_C = maxTslope_C,
                           maxTSE_C = maxTSE_C,
                           maxTp_value_C = maxTp_value_C)


#save table
setwd(wd_tables)
write.csv(temp_all_pts, 'Temperature_Rel_Polar_all_points.csv', row.names = F)



##############################################################################

###############################################
############## EACH SPECIES PREC ##############
###############################################


## Load results for PPT variables 
setwd(wd_res_shap)

#minPPT
sps_minPPT <- lapply(sps_list, function(x){read.csv(paste0(x,'_minPPT_.csv'))})
names(sps_minPPT ) <- sps_list

#meanPPT
sps_meanPPT <- lapply(sps_list,
                      function(x){read.csv(paste0(x,'_meanPPT_.csv'))})
names(sps_meanPPT) <- sps_list

#maxPPT
sps_maxPPT <- lapply(sps_list, function(x){read.csv(paste0(x,'_maxPPT_.csv'))})
names(sps_maxPPT) <- sps_list


## Select only rows representing presence
sps_minPPT_pr <- lapply(sps_minPPT, function(x){x[x$Occurrence == 1,]}) #minPPT
sps_meanPPT_pr <- lapply(sps_meanPPT, function(x){x[x$Occurrence == 1,]}) #meanPPT
sps_maxPPT_pr <- lapply(sps_maxPPT, function(x){x[x$Occurrence == 1,]}) #maxPPT


## Concatenate point and range info to SHAP results

#minPPT
sps_minPPT_SHAP_info <- list()
for(j in 1:length(sps_minPPT_pr))
{
  sps_minPPT_SHAP_info[[j]] <- cbind(sps_minPPT_pr[[j]], 
                                     sps_measures[[j]][,c('key','datasetKey',
                                                     'rangeSize','roundness',
                                                     'centralness',     
                                                     'distEdge',
                                                     'absPolarwardness',
                                                     'relPolarwardness',
                                                     'elevation','biome',
                                                     'bodyMass')])
}

#meanPPT
sps_meanPPT_SHAP_info <- list()
for(j in 1:length(sps_meanPPT_pr))
{
  sps_meanPPT_SHAP_info[[j]] <- cbind(sps_meanPPT_pr[[j]], 
                                      sps_measures[[j]][,c('key','datasetKey',
                                                      'rangeSize','roundness',
                                                      'centralness',     
                                                      'distEdge',
                                                      'absPolarwardness',
                                                      'relPolarwardness',
                                                      'elevation','biome',
                                                      'bodyMass')])
}

#maxPPT
sps_maxPPT_SHAP_info <- list()
for(j in 1:length(sps_maxPPT_pr))
{
  sps_maxPPT_SHAP_info[[j]] <- cbind(sps_maxPPT_pr[[j]], 
                                     sps_measures[[j]][,c('key','datasetKey',
                                                     'rangeSize','roundness',
                                                     'centralness',     
                                                     'distEdge',
                                                     'absPolarwardness',
                                                     'relPolarwardness',
                                                     'elevation','biome',
                                                     'bodyMass')])
}


#create vectors to store the results for each species
species <- character()
order <- character()
bodymass <- numeric()
n_biomes <- numeric()
min_elev <- numeric()
mean_elev <- numeric()
max_elev <- numeric()
order <- character()
bodymass <- numeric()
n_biomes <- numeric()
min_elev <- numeric()
max_elev <- numeric()
n_records <- numeric()
hemisphere <- character()
minPPTR2_RP <- numeric()
minPPTslope_RP <- numeric()
minPPTSE_RP <- numeric()
minPPTp_value_RP <- numeric()
meanPPTR2_RP <- numeric()
meanPPTslope_RP <- numeric()
meanPPTSE_RP <- numeric()
meanPPTp_value_RP <- numeric()
maxPPTR2_RP <- numeric()
maxPPTslope_RP <- numeric()
maxPPTSE_RP <- numeric()
maxPPTp_value_RP <- numeric()
minPPTR2_C <- numeric()
minPPTslope_C <- numeric()
minPPTSE_C <- numeric()
minPPTp_value_C <- numeric()
meanPPTR2_C <- numeric()
meanPPTslope_C <- numeric()
meanPPTSE_C <- numeric()
meanPPTp_value_C <- numeric()
maxPPTR2_C <- numeric()
maxPPTslope_C <- numeric()
maxPPTSE_C <- numeric()
maxPPTp_value_C <- numeric()
rangeSize <- numeric()
rangeRoundness <- numeric()
maxLat <- numeric()
minLat <- numeric()

## Loop through species
for(j in 1:length(sps_list))
{
  #skip sps with only one point
  if(nrow(sps_minT_SHAP_info[[j]]) > 1){
    #list spp
    species[length(species) + 1] <- sps_list[j]
    
    #list order
    order[length(order) + 1] <-
      sps_orders$order[sps_orders$species == sps_list[j]]
    
    #list bodymass
    bodymass[length(bodymass) + 1] <- unique(sps_minPPT_SHAP_info[[j]]$bodyMass)
    
    #list number of biomes
    n_biomes[length(n_biomes) + 1] <-
      length(unique(sps_minPPT_SHAP_info[[j]]$biome))
    
    #list min_elev
    min_elev[length(min_elev) + 1] <- min(sps_minPPT_SHAP_info[[j]]$elevation)
    
    #list mean_elev
    mean_elev[length(mean_elev) + 1] <- mean(sps_minPPT_SHAP_info[[j]]$elevation)
    
    #list max_elev
    max_elev[length(max_elev) + 1] <- max(sps_minPPT_SHAP_info[[j]]$elevation)
    
    #check in which hemisphere the species is
    if(length(unique(sps_minPPT_SHAP_info[[j]]$decimalLatitude < 0)) == 1){
      if(unique(sps_minT_SHAP_info[[j]]$decimalLatitude < 0)){
        hemisphere[length(hemisphere) + 1] <- 'Southern'
      }
      if(unique(sps_minPPT_SHAP_info[[j]]$decimalLatitude > 0)){
        hemisphere[length(hemisphere) + 1] <- 'Northern'
      }
    }
    
    if(length(unique(sps_minPPT_SHAP_info[[j]]$decimalLatitude < 0)) == 2){
      hemisphere[length(hemisphere) + 1] <- 'Both'
    }
    
    ## Get SHAP values of the target variable for each point
    
    #minPPT
    cont_minPPT <- sps_minPPT_SHAP_info[[j]]$Min_PPT_SHAP
    
    #meanPPT
    cont_meanPPT <- sps_meanPPT_SHAP_info[[j]]$Mean_PPT_SHAP
    
    #maxPPT
    cont_maxPPT <- sps_maxPPT_SHAP_info[[j]]$Max_PPT_SHAP
    
    ## Inform folder for the species results
    wd_sps <- paste0(wd_res, '/', gsub(' ', '_', sps_list[j]))
    
    ## Save plot precipitation variables contribution per relative polewardness
    
    #only for species not crossing the equator
    
    if(hemisphere[length(hemisphere)] != 'Both'){
      setwd(wd_sps)
      pdf(file = paste0(gsub(' ', '_', sps_list[j]),
                        '_Precipitation_RelPolewardness.pdf'))
      
      #set parametres for plotting
      par(mar = c(5,5,5,5))
      
      #get y and x lims
      ylim <- range(c(cont_minPPT, cont_meanPPT, cont_maxPPT))
      xlim <- range(sps_minPPT_SHAP_info[[j]]$relPolarwardness) #does not matter which I use
      
      #minPPT
      plot(sps_minPPT_SHAP_info[[j]]$relPolarwardness, cont_minPPT, 
           pch = 19, cex = 0.4, col = '#fc8d5910',
           ylab = 'Precipitation contribution',
           xlab = 'Relative Polewardness',
           ylim = c(ylim[1], ylim[2]))
      
      #fit linear model
      lin_mod_minPPT <- lm(cont_minPPT ~ sps_minPPT_SHAP_info[[j]]$relPolarwardness)
      abline(lin_mod_minPPT, col = '#fc8d59', lwd = 2)
      
      #get and plot R^2
      summa <- summary(lin_mod_minPPT)
      r2 <- summa$r.squared
      text(xlim[1] + ((xlim[2] - xlim[1]) / 25),
           ylim[1],
           expression(bold(paste('R'^2*' = '), sep='')), col = '#fc8d59', pos = 4)
      text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
           ylim[1] - ((ylim[2] - ylim[1]) / 300),
           round(r2,3), col = '#fc8d59', font = 2, pos = 4)
      
      #prepare model results and populate vectors
      summa <- as.data.frame(summa$coefficients)
      
      minPPTR2_RP[length(minPPTR2_RP) + 1] <- r2
      minPPTslope_RP[length(minPPTslope_RP) + 1]  <- summa$Estimate[2]
      minPPTSE_RP[length(minPPTSE_RP) + 1] <- summa$`Std. Error`[2]
      minPPTp_value_RP[length(minPPTp_value_RP) + 1] <- summa$`Pr(>|t|)`[2]
      
      #meanPPT
      points(sps_meanPPT_SHAP_info[[j]]$relPolarwardness, cont_meanPPT, 
             pch = 19, cex = 0.4, col = '#8c510a10')
      
      #fit linear model
      lin_mod_meanPPT <- lm(cont_meanPPT ~ sps_meanPPT_SHAP_info[[j]]$relPolarwardness)
      abline(lin_mod_meanPPT, col = '#8c510a', lwd = 2)
      
      #get and plot R^2
      summa <- summary(lin_mod_meanPPT)
      r2 <- summa$r.squared
      text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
           ylim[1] + ((ylim[2] - ylim[1]) / 22), 
           expression(bold(paste('R'^2*' = '), sep='')), col = '#8c510a', pos = 4)
      text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
           ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
           round(r2,3), col = '#8c510a', font = 2, pos = 4)
      
      #prepare model results and populate vectors
      summa <- as.data.frame(summa$coefficients)
      
      meanPPTR2_RP[length(meanPPTR2_RP) + 1] <- r2
      meanPPTslope_RP[length(meanPPTslope_RP) + 1]  <- summa$Estimate[2]
      meanPPTSE_RP[length(meanPPTSE_RP) + 1] <- summa$`Std. Error`[2]
      meanPPTp_value_RP[length(meanPPTp_value_RP) + 1] <- summa$`Pr(>|t|)`[2]
      
      #maxPPT
      points(sps_maxPPT_SHAP_info[[j]]$relPolarwardness, cont_maxPPT, 
             pch = 19, cex = 0.4, col = '#1a985010')
      
      #fit linear model
      lin_mod_maxPPT <- lm(cont_maxPPT ~ sps_maxPPT_SHAP_info[[j]]$relPolarwardness)
      abline(lin_mod_maxPPT, col = '#1a9850', lwd = 2)
      
      #get and plot R^2
      summa <- summary(lin_mod_maxPPT)
      r2 <- summa$r.squared
      text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
           ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
           expression(bold(paste('R'^2*' = '), sep='')), col = '#1a9850', pos = 4)
      text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
           ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
           round(r2,3), col = '#1a9850', font = 2, pos = 4)
      
      #prepare model results and populate vectors
      summa <- as.data.frame(summa$coefficients)
      
      maxPPTR2_RP[length(maxPPTR2_RP) + 1] <- r2
      maxPPTslope_RP[length(maxPPTslope_RP) + 1]  <- summa$Estimate[2]
      maxPPTSE_RP[length(maxPPTSE_RP) + 1] <- summa$`Std. Error`[2]
      maxPPTp_value_RP[length(maxPPTp_value_RP) + 1] <- summa$`Pr(>|t|)`[2]
      
      #plot species range size and roundness
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5), 
           ylim[2] - ((ylim[2] - ylim[1]) / 25), 
           'Range', font = 2, cex = 0.8, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 2 * ((ylim[2] - ylim[1]) / 25), 
           paste0('Size  ', unique(round(sps_maxT_SHAP_info[[j]]$rangeSize))), 
           cex = 0.7, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 2.8 * ((ylim[2] - ylim[1]) / 25),
           paste0('Roundness  ',
                  unique(round(sps_maxT_SHAP_info[[j]]$roundness, 2))),
           cex = 0.7, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 3.6 * ((ylim[2] - ylim[1]) / 25),
           paste0('Max lat  ',
                  round(max(sps_maxT_SHAP_info[[j]]$decimalLatitude), 2)),
           cex = 0.7, pos = 4)
      text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
           ylim[2] - 4.4 * ((ylim[2] - ylim[1]) / 25),
           paste0('Min lat  ',
                  round(min(sps_minT_SHAP_info[[j]]$decimalLatitude), 2)),
           cex = 0.7, pos = 4)
      
      #plot species name
      text(xlim[1] + ((xlim[2] - xlim[1]) / 20),
           ylim[2] - ((ylim[2] - ylim[1]) / 25),
           sps_list[j], font = 4, pos = 4)
      
      dev.off()
    }else{
      
      minPPTR2_RP[length(minPPTR2_RP) + 1] <- NA
      minPPTslope_RP[length(minPPTslope_RP) + 1] <- NA
      minPPTSE_RP[length(minPPTSE_RP) + 1] <- NA
      minPPTp_value_RP[length(minPPTp_value_RP) + 1] <- NA
      
      meanPPTR2_RP[length(meanPPTR2_RP) + 1] <- NA
      meanPPTslope_RP[length(meanPPTslope_RP) + 1] <- NA
      meanPPTSE_RP[length(meanPPTSE_RP) + 1] <- NA
      meanPPTp_value_RP[length(meanPPTp_value_RP) + 1] <- NA
      
      maxPPTR2_RP[length(maxPPTR2_RP) + 1] <- NA
      maxPPTslope_RP[length(maxPPTslope_RP) + 1] <- NA
      maxPPTSE_RP[length(maxPPTSE_RP) + 1] <- NA
      maxPPTp_value_RP[length(maxPPTp_value_RP) + 1] <- NA
    }
    
    
    ## Save plot precipitation variables contribution per centralness
    
    #### TEST #### change wd after
    #setwd(wd_test)
    setwd(wd_sps)
    pdf(file = paste0(gsub(' ', '_',
                           sps_list[j]),'_Precipitation_DistEdge.pdf'))
    
    #set parametres for plotting
    par(mar = c(5,5,5,5))
    
    #get y and x lims
    ylim <- range(c(cont_minPPT, cont_meanPPT, cont_maxPPT))
    xlim <- range(sps_minPPT_SHAP_info[[j]]$distEdge) #does not matter which I use
    
    #minPPT
    plot(sps_minPPT_SHAP_info[[j]]$distEdge, cont_minPPT, 
         pch = 19, cex = 0.4, col = '#fc8d5910',
         ylab = 'Precipitation contribution',
         xlab = 'Distance from edge (km)',
         ylim = c(ylim[1], ylim[2]))
    
    #fit linear model
    lin_mod_minPPT <- lm(cont_minPPT ~ sps_minPPT_SHAP_info[[j]]$distEdge)
    abline(lin_mod_minPPT, col = '#fc8d59', lwd = 2)
    
    #get and plot R^2
    summa <- summary(lin_mod_minPPT)
    r2 <- summa$r.squared
    text(xlim[1] + ((xlim[2] - xlim[1]) / 25),
         ylim[1],
         expression(bold(paste('R'^2*' = '), sep='')), col = '#fc8d59', pos = 4)
    text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
         ylim[1] - ((ylim[2] - ylim[1]) / 300),
         round(r2,3), col = '#fc8d59', font = 2, pos = 4)
    
    #prepare model results and populate vectors
    summa <- as.data.frame(summa$coefficients)
    
    minPPTR2_C[length(minPPTR2_C) + 1] <- r2
    minPPTslope_C[length(minPPTslope_C) + 1]  <- summa$Estimate[2]
    minPPTSE_C[length(minPPTSE_C) + 1] <- summa$`Std. Error`[2]
    minPPTp_value_C[length(minPPTp_value_C) + 1] <- summa$`Pr(>|t|)`[2]
    
    #meanPPT
    points(sps_meanPPT_SHAP_info[[j]]$distEdge, cont_meanPPT, 
           pch = 19, cex = 0.4, col = '#8c510a10')
    
    #fit linear model
    lin_mod_meanPPT <- lm(cont_meanPPT ~ sps_meanPPT_SHAP_info[[j]]$distEdge)
    abline(lin_mod_meanPPT, col = '#8c510a', lwd = 2)
    
    #get and plot R^2
    summa <- summary(lin_mod_meanPPT)
    r2 <- summa$r.squared
    text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
         ylim[1] + ((ylim[2] - ylim[1]) / 22), 
         expression(bold(paste('R'^2*' = '), sep='')), col = '#8c510a', pos = 4)
    text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
         ylim[1] + ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
         round(r2,3), col = '#8c510a', font = 2, pos = 4)
    
    #prepare model results and populate vectors
    summa <- as.data.frame(summa$coefficients)
    
    meanPPTR2_C[length(meanPPTR2_C) + 1] <- r2
    meanPPTslope_C[length(meanPPTslope_C) + 1]  <- summa$Estimate[2]
    meanPPTSE_C[length(meanPPTSE_C) + 1] <- summa$`Std. Error`[2]
    meanPPTp_value_C[length(meanPPTp_value_C) + 1] <- summa$`Pr(>|t|)`[2]
    
    #maxPPT
    points(sps_maxPPT_SHAP_info[[j]]$distEdge, cont_maxPPT, 
           pch = 19, cex = 0.4, col = '#1a985010')
    
    #fit linear model
    lin_mod_maxPPT <- lm(cont_maxPPT ~ sps_maxPPT_SHAP_info[[j]]$distEdge)
    abline(lin_mod_maxPPT, col = '#1a9850', lwd = 2)
    
    #get and plot R^2
    summa <- summary(lin_mod_maxPPT)
    r2 <- summa$r.squared
    text(xlim[1] + ((xlim[2] - xlim[1]) / 25), 
         ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22), 
         expression(bold(paste('R'^2*' = '), sep='')), col = '#1a9850', pos = 4)
    text(xlim[1] + ((xlim[2] - xlim[1]) / 8),
         ylim[1] + 2 * ((ylim[2] - ylim[1]) / 22) - ((ylim[2] - ylim[1]) / 300), 
         round(r2,3), col = '#1a9850', font = 2, pos = 4)
    
    #prepare model results and populate vectors
    summa <- as.data.frame(summa$coefficients)
    
    maxPPTR2_C[length(maxPPTR2_C) + 1] <- r2
    maxPPTslope_C[length(maxPPTslope_C) + 1]  <- summa$Estimate[2]
    maxPPTSE_C[length(maxPPTSE_C) + 1] <- summa$`Std. Error`[2]
    maxPPTp_value_C[length(maxPPTp_value_C) + 1] <- summa$`Pr(>|t|)`[2]
    
    #plot species range size and roundness
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5), 
         ylim[2] - ((ylim[2] - ylim[1]) / 25), 
         'Range', font = 2, cex = 0.8, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 2 * ((ylim[2] - ylim[1]) / 25), 
         paste0('Size  ', unique(round(sps_maxT_SHAP_info[[j]]$rangeSize))), 
         cex = 0.7, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 2.8 * ((ylim[2] - ylim[1]) / 25),
         paste0('Roundness  ',
                unique(round(sps_maxT_SHAP_info[[j]]$roundness, 2))),
         cex = 0.7, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 3.6 * ((ylim[2] - ylim[1]) / 25),
         paste0('Max lat  ',
                round(max(sps_maxT_SHAP_info[[j]]$decimalLatitude), 2)),
         cex = 0.7, pos = 4)
    text(xlim[2] - ((xlim[2] - xlim[1]) / 4.5),
         ylim[2] - 4.4 * ((ylim[2] - ylim[1]) / 25),
         paste0('Min lat  ',
                round(min(sps_minT_SHAP_info[[j]]$decimalLatitude), 2)),
         cex = 0.7, pos = 4)
    
    #plot species name
    text(xlim[1] + ((xlim[2] - xlim[1]) / 20),
         ylim[2] - ((ylim[2] - ylim[1]) / 25),
         sps_list[j], font = 4, pos = 4)
    
    dev.off()
    
    rangeSize[length(rangeSize) + 1] <- unique(sps_maxPPT_SHAP_info[[j]]$rangeSize)
    rangeRoundness[length(rangeRoundness) + 1] <- 
      unique(sps_maxPPT_SHAP_info[[j]]$roundness) 
    maxLat[length(maxLat) + 1] <- max(abs(sps_maxPPT_SHAP_info[[j]]$decimalLatitude))
    minLat[length(minLat) + 1] <- min(abs(sps_maxPPT_SHAP_info[[j]]$decimalLatitude))
    n_records[length(n_records) + 1] <- length(cont_minPPT)
  }
}


#### COMPLETE THE LIST BELOW


#create a dataframe with the results for all species
prec_all_pts <- data.frame(species = species,
                           order = order,
                           bodymass = bodymass,
                           n_biomes = n_biomes,
                           min_elev = min_elev,
                           mean_elev = mean_elev,
                           max_elev = max_elev,
                           n_records = n_records,
                           hemisphere = hemisphere,
                           rangeSize = rangeSize,
                           rangeRoundness = rangeRoundness,
                           maxLat = maxLat,
                           minLat = minLat,
                           minPPTR2_RP = minPPTR2_RP,
                           minPPTslope_RP = minPPTslope_RP,
                           minPPTSE_RP = minPPTSE_RP,
                           minPPTp_value_RP = minPPTp_value_RP,
                           meanPPTR2_RP = meanPPTR2_RP,
                           meanPPTslope_RP = meanPPTslope_RP,
                           meanPPTSE_RP = meanPPTSE_RP,
                           meanPPTp_value_RP = meanPPTp_value_RP,
                           maxPPTR2_RP = maxPPTR2_RP,
                           maxPPTslope_RP = maxPPTslope_RP,
                           maxPPTSE_RP = maxPPTSE_RP,
                           maxPPTp_value_RP = maxPPTp_value_RP,
                           minPPTR2_C = minPPTR2_C,
                           minPPTslope_C = minPPTslope_C,
                           minPPTSE_C = minPPTSE_C,
                           minPPTp_value_C = minPPTp_value_C,
                           meanPPTR2_C = meanPPTR2_C, 
                           meanPPTslope_C = meanPPTslope_C,
                           meanPPTSE_C = meanPPTSE_C,
                           meanPPTp_value_C = meanPPTp_value_C,
                           maxPPTR2_C = maxPPTR2_C,
                           maxPPTslope_C = maxPPTslope_C,
                           maxPPTSE_C = maxPPTSE_C,
                           maxPPTp_value_C = maxPPTp_value_C)

#save table
setwd(wd_tables)
write.csv(prec_all_pts, 'Precipitation_Rel_Polar_all_points.csv', row.names = F)


