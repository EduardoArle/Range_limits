#load packages 
library(raster)

#list wds
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Thinned Occ'
wd_variables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_biases <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Biases'

#list orders
setwd(wd_thinned_occ)
orders <- gsub('_thin.csv', '', list.files())

#load one variable to get info for the bias raster
setwd(wd_variables)
var <- raster(list.files()[1])

#create coarser layer for the bias raster
var_bias <- aggregate(var, fact = 12, fun = mean)

for(i in 1:length(orders))
{
  #select order
  order <- orders[i]
  
  #load order thinned occurrences
  setwd(wd_thinned_occ)
  ord_occ <- read.csv(paste0(order, '_thin.csv'))
  
  #drop entries with missing coords
  ord_occ <- ord_occ[!is.na(ord_occ$decimalLongitude) &
                     !is.na(ord_occ$decimalLatitude), ]
  
  #select only points with species info
  ord_occ <- ord_occ[ord_occ$species != "",]
  
  #rasterise counts of other-species occurrences
  bias_r <- rasterize(x = ord_occ[,c("decimalLongitude", "decimalLatitude")],
                      y = var_bias,
                      fun = "count",
                      background = 0)
  
  #smooth bias raster
  bias_s <- raster::focal(bias_r, w = matrix(1, 3, 3), fun = sum,
                          na.rm = TRUE, pad = TRUE, padValue = 0)
  
  #normalise to 0–1
  bias_s <- bias_s - minValue(bias_s)
  bias_s <- bias_s / maxValue(bias_s)
  
  #write bias raster
  setwd(wd_biases)
  writeRaster(bias_s,
              filename = paste0(order, "_bias_0.5deg.tif"), overwrite = T)
}
