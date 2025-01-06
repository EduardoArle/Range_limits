#load libraries
library(mgcv); library(itsadug)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20241210_GAMs/Models'
wd_models_rest_dis <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20241210_GAMs/Restricted_distance/Models'

#read results table
setwd(wd_tables)
results <- read.csv('20241210_Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,12)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)

#select only entries up to 200 km from the edge
results_rest_100 <- results[results$distEdge <= 100,]
results_rest_200 <- results[results$distEdge <= 200,]

###########.  GAM   ################

########################
########## T ###########
########################


### minT


#load full model GS
setwd(wd_models)
minT_distEdge_GS <- readRDS('minT_distEdge_GS')

#run model GS with 100 km limit
minT_distEdge_GS_100 <- gam(avg_Min_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_rest_100,
                        method="REML")


summary(minT_distEdge_GS_100)

setwd(wd_models_rest_dis)
saveRDS(minT_distEdge_GS_100, 'minT_distEdge_GS_100')

#run model GS with 200 km limit
minT_distEdge_GS_200 <- gam(avg_Min_T_SHAP ~
                              s(distEdge, k=4, m=2) 
                            + s(distEdge, species, k=4, bs="fs", m=2),
                            data = results_rest_200,
                            method="REML")


summary(minT_distEdge_GS_200)

setwd(wd_models_rest_dis)
saveRDS(minT_distEdge_GS_200, 'minT_distEdge_GS_200')


#plot restricted 100 km model GS 
plot.gam(minT_distEdge_GS_100, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot restricted 200 km model GS 
plot.gam(minT_distEdge_GS_200, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot full model GS restricting to 100 km
plot.gam(minT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         xlim = c(0, 100),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot full model GS restricting to 200 km
plot.gam(minT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         xlim = c(0, 200),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanT


#load full model GS
setwd(wd_models)
meanT_distEdge_GS <- readRDS('meanT_distEdge_GS')

#run model GS with 100 km limit
meanT_distEdge_GS_100 <- gam(avg_Mean_T_SHAP ~
                              s(distEdge, k=4, m=2) 
                            + s(distEdge, species, k=4, bs="fs", m=2),
                            data = results_rest_100,
                            method="REML")


summary(meanT_distEdge_GS_100)

setwd(wd_models_rest_dis)
saveRDS(meanT_distEdge_GS_100, 'meanT_distEdge_GS_100')

#run model GS with 200 km limit
meanT_distEdge_GS_200 <- gam(avg_Mean_T_SHAP ~
                              s(distEdge, k=4, m=2) 
                            + s(distEdge, species, k=4, bs="fs", m=2),
                            data = results_rest_200,
                            method="REML")


summary(meanT_distEdge_GS_200)

setwd(wd_models_rest_dis)
saveRDS(meanT_distEdge_GS_200, 'meanT_distEdge_GS_200')

#plot restricted 100 km model GS
plot.gam(meanT_distEdge_GS_100, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot restricted 200 km model GS 
plot.gam(meanT_distEdge_GS_200, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot full model GS restricting to 100 km
plot.gam(meanT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         xlim = c(0, 100),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot full model GS restricting to 200 km
plot.gam(meanT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         xlim = c(0, 200),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxT


#load full model GS
setwd(wd_models)
maxT_distEdge_GS <- readRDS('maxT_distEdge_GS')

#run model GS with 100 km limit
maxT_distEdge_GS_100 <- gam(avg_Max_T_SHAP ~
                               s(distEdge, k=4, m=2) 
                             + s(distEdge, species, k=4, bs="fs", m=2),
                             data = results_rest_100,
                             method="REML")


summary(maxT_distEdge_GS_100)

setwd(wd_models_rest_dis)
saveRDS(maxT_distEdge_GS_100, 'maxT_distEdge_GS_100')

#run model GS with 200 km limit
maxT_distEdge_GS_200 <- gam(avg_Max_T_SHAP ~
                               s(distEdge, k=4, m=2) 
                             + s(distEdge, species, k=4, bs="fs", m=2),
                             data = results_rest_200,
                             method="REML")


summary(maxT_distEdge_GS_200)

setwd(wd_models_rest_dis)
saveRDS(maxT_distEdge_GS_200, 'maxT_distEdge_GS_200')

#plot restricted 100 km model GS
plot.gam(maxT_distEdge_GS_100, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot restricted 200 km model GS 
plot.gam(maxT_distEdge_GS_200, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot full model GS restricting to 100 km
plot.gam(maxT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         xlim = c(0, 100),
         cex.lab = 2, cex.axis = 1.5) #save 800

#plot full model GS restricting to 200 km
plot.gam(maxT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.07),
         xlim = c(0, 200),
         cex.lab = 2, cex.axis = 1.5) #save 800