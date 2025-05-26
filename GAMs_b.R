#load libraries
library(mgcv); library(itsadug)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20250305_GAMs/Models'

#read results table
setwd(wd_tables)
results <- read.csv('20241210_Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,12)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)



########################
########## T ###########
########################


#select only species that had lower correl between vars
results_minT <- results[abs(results$Cor_vars_minT) <= 0.7,]
results_minT <- results_minT[complete.cases(results_minT$Cor_vars_minT),]
n_sps_minT <- length(unique(results_minT$species))

results_meanT <- results[abs(results$Cor_vars_meanT) <= 0.7,]
results_meanT <- results_meanT[complete.cases(results_meanT$Cor_vars_meanT),]
n_sps_meanT <- length(unique(results_meanT$species))

results_maxT <- results[abs(results$Cor_vars_maxT) <= 0.7,]
results_maxT <- results_maxT[complete.cases(results_maxT$Cor_vars_maxT),]
n_sps_maxT <- length(unique(results_maxT$species))

#make a version of the results with absolute shap values
results_minT_abs <- results_minT
results_meanT_abs <- results_meanT
results_maxT_abs <- results_maxT

results_minT_abs$avg_Min_T_SHAP <- abs(results_minT_abs$avg_Min_T_SHAP)

results_meanT_abs$avg_Mean_T_SHAP <- abs(results_meanT_abs$avg_Mean_T_SHAP)

results_maxT_abs$avg_Max_T_SHAP <- abs(results_maxT_abs$avg_Max_T_SHAP)

#make a version of the results keeping only points up to 250km from edges









# SHAP values X distance from edge (GAM)

#set par for plotting
par(mar = c(6,6,6,6), pty='m')

### minT

#model G
minT_distEdge_G <- gam(avg_Min_T_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_minT, bs="re"),
                       data = results_minT_abs,
                       method="REML",
                       family="gaussian")


summary(minT_distEdge_G)

#AIC
AIC_minT_distEdge_G <- AIC(minT_distEdge_G)

setwd(wd_models)
saveRDS(minT_distEdge_G, 'minT_distEdge_G')


plot.gam(minT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_distEdge_GS <- gam(avg_Min_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_minT_abs,
                        method="REML")


summary(minT_distEdge_GS)

#AIC
AIC_minT_distEdge_GS <- AIC(minT_distEdge_GS)  

setwd(wd_models)
saveRDS(minT_distEdge_GS, 'minT_distEdge_GS')


plot.gam(minT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

summary(results_minT_abs$avg_Min_T_SHAP)


### meanT

#model G
meanT_distEdge_G <- gam(avg_Mean_T_SHAP ~
                          s(distEdge, k=4, bs="tp") 
                        + s(species, k = n_sps_meanT, bs="re"),
                        data = results_meanT,
                        method="REML",
                        family="gaussian")


summary(meanT_distEdge_G)

#AIC
AIC_meanT_distEdge_G <- AIC(meanT_distEdge_G)

setwd(wd_models)
saveRDS(meanT_distEdge_G, 'meanT_distEdge_G')


plot.gam(meanT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanT_distEdge_GS <- gam(avg_Mean_T_SHAP ~
                           s(distEdge, k=4, m=2) 
                         + s(distEdge, species, k=4, bs="fs", m=2),
                         data = results_meanT_abs,
                         method="REML")


summary(meanT_distEdge_GS)

#AIC
AIC_meanT_distEdge_GS <- AIC(meanT_distEdge_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_GS, 'meanT_distEdge_GS')


plot.gam(meanT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


summary(results_meanT_abs$avg_Mean_T_SHAP)

### maxT

#model G
maxT_distEdge_G <- gam(avg_Max_T_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_maxT, bs="re"),
                       data = results_maxT,
                       method="REML",
                       family="gaussian")


summary(maxT_distEdge_G)

#AIC
AIC_maxT_distEdge_G <- AIC(maxT_distEdge_G)

setwd(wd_models)
saveRDS(maxT_distEdge_G, 'maxT_distEdge_G')


plot.gam(maxT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxT_distEdge_GS <- gam(avg_Max_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_maxT_abs,
                        method="REML")


summary(maxT_distEdge_GS)

#AIC
AIC_maxT_distEdge_GS <- AIC(maxT_distEdge_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_GS, 'maxT_distEdge_GS')


plot.gam(maxT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

summary(results_maxT_abs$avg_Max_T_SHAP)



# SHAP values X relative polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G
minT_relPol_G <- gam(avg_Min_T_SHAP ~
                         s(relPolewardness, k=4, bs="tp") 
                       + s(species, k = n_sps_minT, bs="re"),
                       data = results_minT,
                       method="REML",
                       family="gaussian")


summary(minT_relPol_G)

#AIC
AIC_minT_relPol_G <- AIC(minT_relPol_G)

setwd(wd_models)
saveRDS(minT_relPol_G, 'minT_relPol_G')


plot.gam(minT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_relPol_GS <- gam(avg_Min_T_SHAP ~
                          s(relPolewardness, k=4, m=2) 
                        + s(relPolewardness, species, k=4, bs="fs", m=2),
                        data = results_minT,
                        method="REML")


summary(minT_relPol_GS)

#AIC
AIC_minT_relPol_GS <- AIC(minT_relPol_GS)

setwd(wd_models)
saveRDS(minT_relPol_GS, 'minT_relPol_GS')


plot.gam(minT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanT

#model G
meanT_relPol_G <- gam(avg_Mean_T_SHAP ~
                          s(relPolewardness, k=4, bs="tp") 
                        + s(species, k = n_sps_meanT, bs="re"),
                        data = results_meanT,
                        method="REML",
                        family="gaussian")


summary(meanT_relPol_G)

#AIC
AIC_meanT_relPol_G <- AIC(meanT_relPol_G)

setwd(wd_models)
saveRDS(meanT_relPol_G, 'meanT_relPol_G')


plot.gam(meanT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanT_relPol_GS <- gam(avg_Mean_T_SHAP ~
                           s(relPolewardness, k=4, m=2) 
                         + s(relPolewardness, species, k=4, bs="fs", m=2),
                         data = results_meanT,
                         method="REML")


summary(meanT_relPol_GS)

#AIC
AIC_meanT_relPol_GS <- AIC(meanT_relPol_GS)

setwd(wd_models)
saveRDS(meanT_relPol_GS, 'meanT_relPol_GS')


plot.gam(meanT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxT

#model G
maxT_relPol_G <- gam(avg_Max_T_SHAP ~
                         s(relPolewardness, k=4, bs="tp") 
                       + s(species, k = n_sps_maxT, bs="re"),
                       data = results_maxT,
                       method="REML",
                       family="gaussian")


summary(maxT_relPol_G)

#AIC
AIC_maxT_relPol_G <- AIC(maxT_relPol_G)

setwd(wd_models)
saveRDS(maxT_relPol_G, 'maxT_relPol_G')


plot.gam(maxT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxT_relPol_GS <- gam(avg_Max_T_SHAP ~
                          s(relPolewardness, k=4, m=2) 
                        + s(relPolewardness, species, k=4, bs="fs", m=2),
                        data = results_maxT,
                        method="REML")


summary(maxT_relPol_GS)

#AIC
AIC_maxT_relPol_GS <- AIC(maxT_relPol_GS)

setwd(wd_models)
saveRDS(maxT_relPol_GS, 'maxT_relPol_GS')


plot.gam(maxT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800




# SHAP values X absolute polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G
minT_absPol_G <- gam(avg_Min_T_SHAP ~
                       s(absPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_minT, bs="re"),
                     data = results_minT,
                     method="REML",
                     family="gaussian")


summary(minT_absPol_G)

#AIC
AIC_minT_absPol_G <- AIC(minT_absPol_G)

setwd(wd_models)
saveRDS(minT_absPol_G, 'minT_absPol_G')


plot.gam(minT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_absPol_GS <- gam(avg_Min_T_SHAP ~
                        s(absPolewardness, k=4, m=2) 
                      + s(absPolewardness, species, k=4, bs="fs", m=2),
                      data = results_minT,
                      method="REML")


summary(minT_absPol_GS)

#AIC
AIC_minT_absPol_GS <- AIC(minT_absPol_GS)

setwd(wd_models)
saveRDS(minT_absPol_GS, 'minT_absPol_GS')


plot.gam(minT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanT

#model G
meanT_absPol_G <- gam(avg_Mean_T_SHAP ~
                        s(absPolewardness, k=4, bs="tp") 
                      + s(species, k = n_sps_meanT, bs="re"),
                      data = results_meanT,
                      method="REML",
                      family="gaussian")


summary(meanT_absPol_G)

#AIC
AIC_meanT_absPol_G <- AIC(meanT_absPol_G)

setwd(wd_models)
saveRDS(meanT_absPol_G, 'meanT_absPol_G')


plot.gam(meanT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanT_absPol_GS <- gam(avg_Mean_T_SHAP ~
                         s(absPolewardness, k=4, m=2) 
                       + s(absPolewardness, species, k=4, bs="fs", m=2),
                       data = results_meanT,
                       method="REML")


summary(meanT_absPol_GS)

#AIC
AIC_meanT_absPol_GS <- AIC(meanT_absPol_GS)

setwd(wd_models)
saveRDS(meanT_absPol_GS, 'meanT_absPol_GS')


plot.gam(meanT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxT

#model G
maxT_absPol_G <- gam(avg_Max_T_SHAP ~
                       s(absPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_maxT, bs="re"),
                     data = results_maxT,
                     method="REML",
                     family="gaussian")


summary(maxT_absPol_G)

#AIC
AIC_maxT_absPol_G <- AIC(maxT_absPol_G)

setwd(wd_models)
saveRDS(maxT_absPol_G, 'maxT_absPol_G')


plot.gam(maxT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxT_absPol_GS <- gam(avg_Max_T_SHAP ~
                        s(absPolewardness, k=4, m=2) 
                      + s(absPolewardness, species, k=4, bs="fs", m=2),
                      data = results_maxT,
                      method="REML")


summary(maxT_absPol_GS)

#AIC
AIC_maxT_absPol_GS <- AIC(maxT_absPol_GS)

setwd(wd_models)
saveRDS(maxT_absPol_GS, 'maxT_absPol_GS')


plot.gam(maxT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800




# SHAP values X elevation (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G
minT_elev_G <- gam(avg_Min_T_SHAP ~
                       s(elevation, k=4, bs="tp") 
                     + s(species, k = n_sps_minT, bs="re"),
                     data = results_minT,
                     method="REML",
                     family="gaussian")


summary(minT_elev_G)

#AIC
AIC_minT_elev_G <- AIC(minT_elev_G)

setwd(wd_models)
saveRDS(minT_elev_G, 'minT_elev_G')


plot.gam(minT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.25),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minT_elev_GS <- gam(avg_Min_T_SHAP ~
                        s(elevation, k=4, m=2) 
                      + s(elevation, species, k=4, bs="fs", m=2),
                      data = results_minT,
                      method="REML")


summary(minT_elev_GS)

#AIC
AIC_minT_elev_GS <- AIC(minT_elev_GS)

setwd(wd_models)
saveRDS(minT_elev_GS, 'minT_elev_GS')


plot.gam(minT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.25),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanT

#model G
meanT_elev_G <- gam(avg_Mean_T_SHAP ~
                        s(elevation, k=4, bs="tp") 
                      + s(species, k = n_sps_meanT, bs="re"),
                      data = results_meanT,
                      method="REML",
                      family="gaussian")


summary(meanT_elev_G)

#AIC
AIC_meanT_elev_G <- AIC(meanT_elev_G)

setwd(wd_models)
saveRDS(meanT_elev_G, 'meanT_elev_G')


plot.gam(meanT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.05, 0.25),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanT_elev_GS <- gam(avg_Mean_T_SHAP ~
                         s(elevation, k=4, m=2) 
                       + s(elevation, species, k=4, bs="fs", m=2),
                       data = results_meanT,
                       method="REML")


summary(meanT_elev_GS)

#AIC
AIC_meanT_elev_GS <- AIC(meanT_elev_GS)

setwd(wd_models)
saveRDS(meanT_elev_GS, 'meanT_elev_GS')


plot.gam(meanT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.5, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxT

#model G
maxT_elev_G <- gam(avg_Max_T_SHAP ~
                       s(elevation, k=4, bs="tp") 
                     + s(species, k = n_sps_maxT, bs="re"),
                     data = results_maxT,
                     method="REML",
                     family="gaussian")


summary(maxT_elev_G)

#AIC
AIC_maxT_elev_G <- AIC(maxT_elev_G)

setwd(wd_models)
saveRDS(maxT_elev_G, 'maxT_elev_G')


plot.gam(maxT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxT_elev_GS <- gam(avg_Max_T_SHAP ~
                        s(elevation, k=4, m=2) 
                      + s(elevation, species, k=4, bs="fs", m=2),
                      data = results_maxT,
                      method="REML")


summary(maxT_elev_GS)

#AIC
AIC_maxT_elev_GS <- AIC(maxT_elev_GS)

setwd(wd_models)
saveRDS(maxT_elev_GS, 'maxT_elev_GS')


plot.gam(maxT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



########################
######### PPT ##########
########################


#select only species that had lower correl between vars
results_minPPT <- results[abs(results$Cor_vars_minPPT) <= 0.7,]
results_minPPT <- results_minPPT[complete.cases
                                 (results_minPPT$Cor_vars_minPPT),]
n_sps_minPPT <- length(unique(results_minPPT$species))

results_meanPPT <- results[abs(results$Cor_vars_meanPPT) <= 0.7,]
results_meanPPT <- results_meanPPT[complete.cases
                                   (results_meanPPT$Cor_vars_meanPPT),]
n_sps_meanPPT <- length(unique(results_meanPPT$species))

results_maxPPT <- results[abs(results$Cor_vars_maxPPT) <= 0.7,]
results_maxPPT <- results_maxPPT[complete.cases
                                 (results_maxPPT$Cor_vars_maxPPT),]
n_sps_maxPPT <- length(unique(results_maxPPT$species))

#make a version of the results with absolute shap values
results_minPPT_abs <- results_minPPT
results_meanPPT_abs <- results_meanPPT
results_maxPPT_abs <- results_maxPPT

results_minPPT_abs$avg_Min_PPT_SHAP <- abs(results_minPPT_abs$avg_Min_PPT_SHAP)

results_meanPPT_abs$avg_Mean_PPT_SHAP<-abs(results_meanPPT_abs$avg_Mean_PPT_SHAP)

results_maxPPT_abs$avg_Max_PPT_SHAP <- abs(results_maxPPT_abs$avg_Max_PPT_SHAP)

# SHAP values X distance from edge (GAM)

#set par for plotting
par(mar = c(6,6,6,6), pty='m')

### minPPT

#model G
minPPT_distEdge_G <- gam(avg_Min_PPT_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_minPPT, bs="re"),
                       data = results_minPPT,
                       method="REML",
                       family="gaussian")


summary(minPPT_distEdge_G)

#AIC
AIC_minPPT_distEdge_G <- AIC(minPPT_distEdge_G)

setwd(wd_models)
saveRDS(minPPT_distEdge_G, 'minPPT_distEdge_G')


plot.gam(minPPT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.05, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_distEdge_GS <- gam(avg_Min_PPT_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_minPPT_abs,
                        method="REML")


summary(minPPT_distEdge_GS)

#AIC
AIC_minPPT_distEdge_GS <- AIC(minPPT_distEdge_GS)

setwd(wd_models)
saveRDS(minPPT_distEdge_GS, 'minPPT_distEdge_GS')


plot.gam(minPPT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.05, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

summary(results_minPPT_abs$avg_Min_PPT_SHAP)


### meanPPT

#model G
meanPPT_distEdge_G <- gam(avg_Mean_PPT_SHAP ~
                          s(distEdge, k=4, bs="tp") 
                        + s(species, k = n_sps_meanPPT, bs="re"),
                        data = results_meanPPT,
                        method="REML",
                        family="gaussian")


summary(meanPPT_distEdge_G)

#AIC
AIC_meanPPT_distEdge_G <- AIC(meanPPT_distEdge_G)

setwd(wd_models)
saveRDS(meanPPT_distEdge_G, 'meanPPT_distEdge_G')


plot.gam(meanPPT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.05),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanPPT_distEdge_GS <- gam(avg_Mean_PPT_SHAP ~
                           s(distEdge, k=4, m=2) 
                         + s(distEdge, species, k=4, bs="fs", m=2),
                         data = results_meanPPT_abs,
                         method="REML")


summary(meanPPT_distEdge_GS)

#AIC
AIC_meanPPT_distEdge_GS <- AIC(meanPPT_distEdge_GS)

setwd(wd_models) 
saveRDS(meanPPT_distEdge_GS, 'meanPPT_distEdge_GS')


plot.gam(meanPPT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.05),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxPPT

#model G
maxPPT_distEdge_G <- gam(avg_Max_PPT_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_maxPPT, bs="re"),
                       data = results_maxPPT,
                       method="REML",
                       family="gaussian")


summary(maxPPT_distEdge_G)

#AIC
AIC_maxPPT_distEdge_G <- AIC(maxPPT_distEdge_G)

setwd(wd_models)
saveRDS(maxPPT_distEdge_G, 'maxPPT_distEdge_G')


plot.gam(maxPPT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxPPT_distEdge_GS <- gam(avg_Max_PPT_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_maxPPT,
                        method="REML")


summary(maxPPT_distEdge_GS)

#AIC
AIC_maxPPT_distEdge_GS <- AIC(maxPPT_distEdge_GS)

setwd(wd_models)
saveRDS(maxPPT_distEdge_GS, 'maxPPT_distEdge_GS')


plot.gam(maxPPT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800





# SHAP values X relative polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minPPT

#model G
minPPT_relPol_G <- gam(avg_Min_PPT_SHAP ~
                       s(relPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_minPPT, bs="re"),
                     data = results_minPPT,
                     method="REML",
                     family="gaussian")


summary(minPPT_relPol_G)

#AIC
AIC_minPPT_relPol_G <- AIC(minPPT_relPol_G)

setwd(wd_models)
saveRDS(minPPT_relPol_G, 'minPPT_relPol_G')


plot.gam(minPPT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_relPol_GS <- gam(avg_Min_PPT_SHAP ~
                        s(relPolewardness, k=4, m=2) 
                      + s(relPolewardness, species, k=4, bs="fs", m=2),
                      data = results_minPPT,
                      method="REML")


summary(minPPT_relPol_GS)

#AIC
AIC_minPPT_relPol_GS <- AIC(minPPT_relPol_GS)

setwd(wd_models)
saveRDS(minPPT_relPol_GS, 'minPPT_relPol_GS')


plot.gam(minPPT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanPPT

#model G
meanPPT_relPol_G <- gam(avg_Mean_PPT_SHAP ~
                        s(relPolewardness, k=4, bs="tp") 
                      + s(species, k = n_sps_meanPPT, bs="re"),
                      data = results_meanPPT,
                      method="REML",
                      family="gaussian")


summary(meanPPT_relPol_G)

#AIC
AIC_meanPPT_relPol_G <- AIC(meanPPT_relPol_G)

setwd(wd_models)
saveRDS(meanPPT_relPol_G, 'meanPPT_relPol_G')


plot.gam(meanPPT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanPPT_relPol_GS <- gam(avg_Mean_PPT_SHAP ~
                         s(relPolewardness, k=4, m=2) 
                       + s(relPolewardness, species, k=4, bs="fs", m=2),
                       data = results_meanPPT,
                       method="REML")


summary(meanPPT_relPol_GS)

#AIC
AIC_meanPPT_relPol_GS <- AIC(meanPPT_relPol_GS)

setwd(wd_models)
saveRDS(meanPPT_relPol_GS, 'meanPPT_relPol_GS')


plot.gam(meanPPT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxPPT

#model G
maxPPT_relPol_G <- gam(avg_Max_PPT_SHAP ~
                       s(relPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_maxPPT, bs="re"),
                     data = results_maxPPT,
                     method="REML",
                     family="gaussian")


summary(maxPPT_relPol_G)

#AIC
AIC_maxPPT_relPol_G <- AIC(maxPPT_relPol_G)

setwd(wd_models)
saveRDS(maxPPT_relPol_G, 'maxPPT_relPol_G')


plot.gam(maxPPT_relPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxPPT_relPol_GS <- gam(avg_Max_PPT_SHAP ~
                        s(relPolewardness, k=4, m=2) 
                      + s(relPolewardness, species, k=4, bs="fs", m=2),
                      data = results_maxPPT,
                      method="REML")


summary(maxPPT_relPol_GS)

#AIC
AIC_maxPPT_relPol_GS <- AIC(maxPPT_relPol_GS)

setwd(wd_models)
saveRDS(maxPPT_relPol_GS, 'maxPPT_relPol_GS')


plot.gam(maxPPT_relPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800




# SHAP values X absolute polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minPPT

#model G
minPPT_absPol_G <- gam(avg_Min_PPT_SHAP ~
                       s(absPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_minPPT, bs="re"),
                     data = results_minPPT,
                     method="REML",
                     family="gaussian")


summary(minPPT_absPol_G)

#AIC
AIC_minPPT_absPol_G <- AIC(minPPT_absPol_G)

setwd(wd_models)
saveRDS(minPPT_absPol_G, 'minPPT_absPol_G')


plot.gam(minPPT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.08, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_absPol_GS <- gam(avg_Min_PPT_SHAP ~
                        s(absPolewardness, k=4, m=2) 
                      + s(absPolewardness, species, k=4, bs="fs", m=2),
                      data = results_minPPT,
                      method="REML")


summary(minPPT_absPol_GS)

#AIC
AIC_minPPT_absPol_GS <- AIC(minPPT_absPol_GS)

setwd(wd_models)
saveRDS(minPPT_absPol_GS, 'minPPT_absPol_GS')


plot.gam(minPPT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.08, 0.08),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanPPT

#model G
meanPPT_absPol_G <- gam(avg_Mean_PPT_SHAP ~
                        s(absPolewardness, k=4, bs="tp") 
                      + s(species, k = n_sps_meanPPT, bs="re"),
                      data = results_meanPPT,
                      method="REML",
                      family="gaussian")


summary(meanPPT_absPol_G)

#AIC
AIC_meanPPT_absPol_G <- AIC(meanPPT_absPol_G)

setwd(wd_models)
saveRDS(meanPPT_absPol_G, 'meanPPT_absPol_G')


plot.gam(meanPPT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.25, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanPPT_absPol_GS <- gam(avg_Mean_PPT_SHAP ~
                         s(absPolewardness, k=4, m=2) 
                       + s(absPolewardness, species, k=4, bs="fs", m=2),
                       data = results_meanPPT,
                       method="REML")


summary(meanPPT_absPol_GS)

#AIC
AIC_meanPPT_absPol_GS <- AIC(meanPPT_absPol_GS)

setwd(wd_models)
saveRDS(meanPPT_absPol_GS, 'meanPPT_absPol_GS')


plot.gam(meanPPT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.25, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxPPT

#model G
maxPPT_absPol_G <- gam(avg_Max_PPT_SHAP ~
                       s(absPolewardness, k=4, bs="tp") 
                     + s(species, k = n_sps_maxPPT, bs="re"),
                     data = results_maxPPT,
                     method="REML",
                     family="gaussian")


summary(maxPPT_absPol_G)

#AIC
AIC_maxPPT_absPol_G <- AIC(maxPPT_absPol_G)

setwd(wd_models)
saveRDS(maxPPT_absPol_G, 'maxPPT_absPol_G')


plot.gam(maxPPT_absPol_G, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxPPT_absPol_GS <- gam(avg_Max_PPT_SHAP ~
                        s(absPolewardness, k=4, m=2) 
                      + s(absPolewardness, species, k=4, bs="fs", m=2),
                      data = results_maxPPT,
                      method="REML")


summary(maxPPT_absPol_GS)

#AIC
AIC_maxPPT_absPol_GS <- AIC(maxPPT_absPol_GS)

setwd(wd_models)
saveRDS(maxPPT_absPol_GS, 'maxPPT_absPol_GS')


plot.gam(maxPPT_absPol_GS, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



# SHAP values X elevation (GAM)


### minT

#model G
minPPT_elev_G <- gam(avg_Min_PPT_SHAP ~
                     s(elevation, k=4, bs="tp") 
                   + s(species, k = n_sps_minPPT, bs="re"),
                   data = results_minPPT,
                   method="REML",
                   family="gaussian")


summary(minPPT_elev_G)

#AIC
AIC_minPPT_elev_G <- AIC(minPPT_elev_G)

#gam.check
check <- gam.check(minPPT_elev_G)

setwd(wd_models)
saveRDS(minPPT_elev_G, 'minPPT_elev_G')


plot.gam(minPPT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.05, 0.05),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_elev_GS <- gam(avg_Min_PPT_SHAP ~
                      s(elevation, k=4, m=2) 
                    + s(elevation, species, k=4, bs="fs", m=2),
                    data = results_minPPT,
                    method="REML")


summary(minPPT_elev_GS)

#AIC
AIC_minPPT_elev_GS <- AIC(minPPT_elev_GS)

#gam.check
check <- gam.check(minPPT_elev_GS)

setwd(wd_models)
saveRDS(minPPT_elev_GS, 'minPPT_elev_GS')


plot.gam(minPPT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.08, 0.45),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanPPT

#model G
meanPPT_elev_G <- gam(avg_Mean_PPT_SHAP ~
                      s(elevation, k=4, bs="tp") 
                    + s(species, k = n_sps_meanPPT, bs="re"),
                    data = results_meanPPT,
                    method="REML",
                    family="gaussian")


summary(meanPPT_elev_G)

#AIC
AIC_meanPPT_elev_G <- AIC(meanPPT_elev_G)

#gam.check
check <- gam.check(meanPPT_elev_G)

setwd(wd_models)
saveRDS(meanPPT_elev_G, 'meanPPT_elev_G')


plot.gam(meanPPT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.05, 0.25),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanPPT_elev_GS <- gam(avg_Mean_PPT_SHAP ~
                       s(elevation, k=4, m=2) 
                     + s(elevation, species, k=4, bs="fs", m=2),
                     data = results_meanPPT,
                     method="REML")


summary(meanPPT_elev_GS)

#AIC
AIC_meanPPT_elev_GS <- AIC(meanPPT_elev_GS)

setwd(wd_models)
saveRDS(meanPPT_elev_GS, 'meanPPT_elev_GS')


plot.gam(meanT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.5, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxPPT

#model G
maxPPT_elev_G <- gam(avg_Max_PPT_SHAP ~
                     s(elevation, k=4, bs="tp") 
                   + s(species, k = n_sps_maxPPT, bs="re"),
                   data = results_maxPPT,
                   method="REML",
                   family="gaussian")


summary(maxPPT_elev_G)

#AIC
AIC_maxPPT_elev_G <- AIC(maxPPT_elev_G)

setwd(wd_models)
saveRDS(maxPPT_elev_G, 'maxPPT_elev_G')


plot.gam(maxPPT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxPPT_elev_GS <- gam(avg_Max_PPT_SHAP ~
                      s(elevation, k=4, m=2) 
                    + s(elevation, species, k=4, bs="fs", m=2),
                    data = results_maxPPT,
                    method="REML")



summary(maxPPT_elev_GS)

#AIC
AIC_maxPPT_elev_GS <- AIC(maxPPT_elev_GS)

setwd(wd_models)
saveRDS(maxPPT_elev_GS, 'maxPPT_elev_GS')


plot.gam(maxPPT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800








######################################################################
# SHAP values X relative polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minPREC

#model G
minPPT_relPolewarness_G <- gam(avg_Min_PPT_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")

summary(minPPT_relPolewarness_G)

setwd(wd_models)
saveRDS(minPPT_relPolewarness_G, 'minPPT_relPolewarness_G')

plot.gam(minPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_relPolewarness_GS <- gam(avg_Min_PPT_SHAP ~
                                s(relPolewardness, k=4, m=2) 
                              + s(relPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")

summary(minPPT_relPolewarness_GS)

setwd(wd_models)
saveRDS(minPPT_relPolewarness_GS, 'minPPT_relPolewarness_GS')

plot.gam(minPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


### meanPPT

#model G
meanPPT_relPolewarness_G <- gam(avg_Mean_PPT_SHAP ~
                               s(relPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")

summary(meanPPT_relPolewarness_G)

setwd(wd_models)
saveRDS(meanPPT_relPolewarness_G, 'meanPPT_relPolewarness_G')

plot.gam(meanPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
meanPPT_relPolewarness_GS <- gam(avg_Mean_PPT_SHAP ~
                                s(relPolewardness, k=4, m=2) 
                              + s(relPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")

summary(meanPPT_relPolewarness_GS)

setwd(wd_models)
saveRDS(meanPPT_relPolewarness_GS, 'meanPPT_relPolewarness_GS')

plot.gam(meanPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.04, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


### maxPPT

#model G
maxPPT_relPolewarness_G <- gam(avg_Max_PPT_SHAP ~
                                s(relPolewardness, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")

summary(maxPPT_relPolewarness_G)

setwd(wd_models)
saveRDS(maxPPT_relPolewarness_G, 'maxPPT_relPolewarness_G')


#########################################
######### Fix colour (to green) #########
#########################################



plot.gam(maxPPT_relPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
maxPPT_relPolewarness_GS <- gam(avg_Max_PPT_SHAP ~
                                 s(relPolewardness, k=4, m=2) 
                               + s(relPolewardness, species, k=4, bs="fs", m=2),
                               data = results,
                               method="REML")


summary(maxPPT_relPolewarness_GS)

setwd(wd_models)
saveRDS(maxPPT_relPolewarness_GS, 'maxPPT_relPolewarness_GS')

#########################################
######### Fix colour (to green) #########
#########################################


plot.gam(maxPPT_relPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.03, 0.03),
         cex.lab = 2, cex.axis = 1.5) #save 800



# SHAP values X absolute polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minPPT

#model G
minPPT_absPolewarness_G <- gam(avg_Min_PPT_SHAP ~
                               s(absPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(minPPT_absPolewarness_G)

setwd(wd_models)
saveRDS(minPPT_absPolewarness_G, 'minPPT_absPolewarness_G')

#########################################
##### Fix colour (to '#fc8d5930') #######
#########################################

plot.gam(minPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.11, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS
minPPT_absPolewarness_GS <- gam(avg_Min_PPT_SHAP ~
                                s(absPolewardness, k=4, m=2) 
                              + s(absPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(minPPT_absPolewarness_GS)
 
setwd(wd_models)
saveRDS(minPPT_absPolewarness_GS, 'minPPT_absPolewarness_GS')


#########################################
##### Fix colour (to '#fc8d5930') #######
#########################################


plot.gam(minPPT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.11, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800



### meanPPT

#model G
meanPPT_absPolewarness_G <- gam(avg_Mean_PPT_SHAP ~
                                s(absPolewardness, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")


summary(meanPPT_absPolewarness_G)

setwd(wd_models)
saveRDS(meanPPT_absPolewarness_G, 'meanPPT_absPolewarness_G')


plot.gam(meanPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.11, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
meanPPT_absPolewarness_GS <- gam(avg_Mean_PPT_SHAP ~
                                 s(absPolewardness, k=4, m=2) 
                               + s(absPolewardness, species, k=4, bs="fs", m=2),
                               data = results,
                               method="REML")

summary(meanPPT_absPolewarness_GS)

setwd(wd_models)
saveRDS(meanPPT_absPolewarness_GS, 'meanPPT_absPolewarness_GS')

plot.gam(meanPPT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.11, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800








### maxT

#model G
maxPPT_absPolewarness_G <- gam(Max_PPT_SHAP ~
                               s(absPolewardness, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(maxPPT_absPolewarness_G)

setwd(wd_models)
saveRDS(maxPPT_absPolewarness_G, 'maxT_absPPPolewarness_G')


plot.gam(maxPPT_absPolewarness_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.11, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_absPolewarness_G, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS
maxPPT_absPolewarness_GS <- gam(Max_PPT_SHAP ~
                                s(absPolewardness, k=4, m=2) 
                              + s(absPolewardness, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(maxT_absPolewarness_GS)

setwd(wd_models)
saveRDS(maxT_absPolewarness_GS, 'maxT_absPolewarness_GS')

draw(maxT_absPolewarness_GS, page = 1)


plot.gam(maxT_absPolewarness_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.11, 0.11),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_absPolewarness_GS, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model S
maxT_absPolewarness_S <- gam(Max_T_SHAP ~
                               s(absPolewardness, species, k=4, bs="fs", m=2),
                             data = results,
                             method="REML")

summary(maxT_absPolewarness_S)

setwd(wd_models)
saveRDS(maxT_absPolewarness_S, 'maxT_absPolewarness_S')

draw(maxT_absPolewarness_S, page = 1)

plot.gam(maxT_absPolewarness_S, select = 1, residuals = F, shade = T,
         cex.lab = 2, cex.axis = 1.5) #save 800




# SHAP values X elevation (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G
minT_elev_G <- gam(Min_T_SHAP ~
                         s(elevation, k=4, bs="tp") 
                       + s(species, k=503, bs="re"),
                       data = results,
                       method="REML",
                       family="gaussian")


summary(minT_elev_G)

setwd(wd_models)
saveRDS(minT_elev_G, 'minT_elev_G')

draw(minT_elev_G, page = 1)

plot.gam(minT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_elev_G, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS
minT_elev_GS <- gam(Min_T_SHAP ~
                          s(elevation, k=4, m=2) 
                        + s(elevation, species, k=4, bs="fs", m=2),
                        data = results,
                        method="REML")


summary(minT_elev_GS)

setwd(wd_models)
saveRDS(minT_elev_GS, 'minT_elev_GS')

draw(minT_elev_GS, page = 1)


plot.gam(minT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_elev_GS, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model S
minT_elev_S <- gam(Min_T_SHAP ~
                         s(elevation, species, k=4, bs="fs", m=2),
                       data = results,
                       method="REML")

summary(minT_elev_S)

setwd(wd_models)
saveRDS(minT_elev_S, 'minT_elev_S')

draw(minT_elev_S, page = 1)

plot.gam(minT_elev_S, select = 1, residuals = F, shade = T,
         cex.lab = 2, cex.axis = 1.5) #save 800

### meanT

#model G
meanT_elev_G <- gam(Mean_T_SHAP ~
                          s(elevation, k=4, bs="tp") 
                        + s(species, k=503, bs="re"),
                        data = results,
                        method="REML",
                        family="gaussian")


summary(meanT_elev_G)

setwd(wd_models)
saveRDS(meanT_elev_G, 'meanT_elev_G')

draw(meanT_elev_G, page = 1)

plot.gam(meanT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_elev_G, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS
meanT_elev_GS <- gam(Mean_T_SHAP ~
                           s(elevation, k=4, m=2) 
                         + s(elevation, species, k=4, bs="fs", m=2),
                         data = results,
                         method="REML")


summary(meanT_elev_GS)

setwd(wd_models)
saveRDS(meanT_elev_GS, 'meanT_elev_GS')

draw(meanT_elev_GS, page = 1)


plot.gam(meanT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_elev_GS, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model S
meanT_elev_S <- gam(Mean_T_SHAP ~
                          s(elevation, species, k=4, bs="fs", m=2),
                        data = results,
                        method="REML")

summary(meanT_elev_S)

setwd(wd_models)
saveRDS(meanT_elev_S, 'meanT_elev_S')

draw(meanT_elev_S, page = 1)

plot.gam(meanT_elev_S, select = 1, residuals = F, shade = T,
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxT

#model G
maxT_elev_G <- gam(Max_T_SHAP ~
                         s(elevation, k=4, bs="tp") 
                       + s(species, k=503, bs="re"),
                       data = results,
                       method="REML",
                       family="gaussian")


summary(maxT_elev_G)

setwd(wd_models)
saveRDS(maxT_elev_G, 'maxT_elev_G')

draw(maxT_elev_G, page = 1)

plot.gam(maxT_elev_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_elev_G, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS
maxT_elev_GS <- gam(Max_T_SHAP ~
                          s(elevation, k=4, m=2) 
                        + s(elevation, species, k=4, bs="fs", m=2),
                        data = results,
                        method="REML")


summary(maxT_elev_GS)

setwd(wd_models)
saveRDS(maxT_elev_GS, 'maxT_elev_GS')

draw(maxT_elev_GS, page = 1)


plot.gam(maxT_elev_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_elev_GS, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model S
maxT_elev_S <- gam(Max_T_SHAP ~
                         s(elevation, species, k=4, bs="fs", m=2),
                       data = results,
                       method="REML")

summary(maxT_elev_S)

setwd(wd_models)
saveRDS(maxT_elev_S, 'maxT_elev_S')

draw(maxT_elev_S, page = 1)

plot.gam(maxT_elev_S, select = 1, residuals = F, shade = T,
         cex.lab = 2, cex.axis = 1.5) #save 800


# ### meanT
# gam_meanT_relPolewarness <- gam(Mean_T_SHAP ~ s(relPolewardness, k = 4),
#                                 na.action = 'na.omit',
#                                 data = results)
# 
# summary(gam_meanT_relPolewarness)
# 
# 
# plot.gam(gam_meanT_relPolewarness, select = 1, residuals = F, shade = T,
#          shade.col = '#80008030', ylab = 'SHAP value',
#          ylim = c(-0.06, 0.06),
#          cex.lab = 2, cex.axis = 1.5)
# 
# #including random effects
# gam_meanT_relPolewarness_RE <- gam(Mean_T_SHAP ~ s(relPolewardness, k = 4) +
#                                     s(species, bs = "re"),
#                                   data = results)
# 
# summary(gam_meanT_relPolewarness_RE)
# 
# plot.gam(gam_meanT_relPolewarness_RE, select = 1, residuals = F, shade = T,
#          shade.col = '#80008030', ylab = 'SHAP value',
#          ylim = c(-0.06, 0.06),
#          cex.lab = 2, cex.axis = 1.5)
# 
# 
# ### maxT
# gam_maxT_relPolewarness <- gam(Max_T_SHAP ~ s(relPolewardness, k = 4),
#                                na.action = 'na.omit',
#                                data = results)
# 
# summary(gam_maxT_relPolewarness)
# 
# plot.gam(gam_maxT_relPolewarness, select = 1, residuals = F, shade = T,
#          shade.col = '#FF000030', ylab = 'SHAP value',
#          ylim = c(-0.06, 0.06),
#          cex.lab = 2, cex.axis = 1.5)
# 
# #including random effects
# gam_maxT_relPolewarness_RE <- gam(Max_T_SHAP ~ s(relPolewardness, k = 4) +
#                                      s(species, bs = "re"),
#                                    data = results)
# 
# summary(gam_maxT_relPolewarness_RE)
# 
# plot.gam(gam_maxT_relPolewarness_RE, select = 1, residuals = F, shade = T,
#          shade.col = '#FF000030', ylab = 'SHAP value',
#          ylim = c(-0.06, 0.06),
#          cex.lab = 2, cex.axis = 1.5)
# 
# 
# # abs latitude
# 
# 
# ### minT
# gam_minT_absPolewarness <- gam(Min_T_SHAP ~ s(abs(decimalLatitude), k = 4),
#                                na.action = 'na.omit',
#                                data = results)
# 
# plot.gam(gam_minT_absPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#0000FF30', ylab = 'SHAP value', xlab = 'Latitude',
#          ylim = c(-0.05, 0.05),
#          cex.lab = 2, cex.axis = 2)
# 
# 
# 
# ### meanT
# gam_meanT_absPolewarness <- gam(Mean_T_SHAP ~ s(abs(decimalLatitude), k = 4),
#                                 na.action = 'na.omit',
#                                 data = results)
# 
# plot.gam(gam_meanT_absPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#80008030', ylab = 'SHAP value', xlab = 'Latitude',
#          ylim = c(-0.05, 0.05),
#          cex.lab = 2, cex.axis = 2)
# 
# ### maxT
# gam_maxT_absPolewarness <- gam(Max_T_SHAP ~ s(abs(decimalLatitude), k = 4),
#                                na.action = 'na.omit',
#                                data = results)
# 
# plot.gam(gam_maxT_absPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF000030', ylab = 'SHAP value', xlab = 'Latitude',
#          ylim = c(-0.05, 0.05),
#          cex.lab = 2, cex.axis = 2)
# 
# 
# 
# # dist edge 
# 
# ### minT
# gam_minT_distEdge <- gam(Min_T_SHAP ~ s(distEdge, k = 4),
#                          na.action = 'na.omit',
#                          data = results)
# 
# plot.gam(gam_minT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#0000FF30', ylab = 'SHAP value')
# 
# 
# 
# ### meanT
# gam_meanT_distEdge <- gam(Mean_T_SHAP ~ s(distEdge, k = 4),
#                           na.action = 'na.omit',
#                           data = results)
# 
# plot.gam(gam_meanT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#80008030', ylab = 'SHAP value')
# 
# ### maxT
# gam_maxT_distEdge <- gam(Max_T_SHAP ~ s(distEdge, k = 4),
#                          na.action = 'na.omit',
#                          data = results)
# 
# plot.gam(gam_maxT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF000030', ylab = 'SHAP value')
# 
# 
# # dist edge (consider only points not too far from the edge)
# 
# results2 <- results[results$distEdge <= 100,]
# 
# ### minT
# gam_minT_distEdge <- gam(Min_T_SHAP ~ s(distEdge, k = 4),
#                          na.action = 'na.omit',
#                          data = results2)
# 
# plot.gam(gam_minT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#0000FF30', ylab = 'SHAP value')
# 
# 
# 
# ### meanT
# gam_meanT_distEdge <- gam(Mean_T_SHAP ~ s(distEdge, k = 4),
#                           na.action = 'na.omit',
#                           data = results2)
# 
# plot.gam(gam_meanT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#80008030', ylab = 'SHAP value')
# 
# ### maxT
# gam_maxT_distEdge <- gam(Max_T_SHAP ~ s(distEdge, k = 4),
#                          na.action = 'na.omit',
#                          data = results2)
# 
# plot.gam(gam_maxT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF000030', ylab = 'SHAP value')
# 
# 
# 
# # elevation
# 
# ### minT
# gam_minT_elev <- gam(Min_T_SHAP ~ s(elevation, k = 4),
#                      na.action = 'na.omit',
#                      data = results)
# 
# plot.gam(gam_minT_elev, pages = 1, residuals = F, shade = T,
#          shade.col = '#0000FF30', ylab = 'SHAP value')
# 
# 
# 
# ### meanT
# gam_meanT_elev <- gam(Mean_T_SHAP ~ s(elevation, k = 4),
#                       na.action = 'na.omit',
#                       data = results)
# 
# plot.gam(gam_meanT_elev, pages = 1, residuals = F, shade = T,
#          shade.col = '#80008030', ylab = 'SHAP value')
# 
# ### maxT
# gam_maxT_elev <- gam(Max_T_SHAP ~ s(elevation, k = 4),
#                      na.action = 'na.omit',
#                      data = results)
# 
# plot.gam(gam_maxT_elev, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF000030', ylab = 'SHAP value')
# 
# 
# 
# 
# #####  PPT
# 
# 
# # rel polewardness
# 
# 
# ### minPPT
# gam_minPPT_relPolewarness <- gam(Min_PPT_SHAP ~ s(relPolewardness, k = 4),
#                                  na.action = 'na.omit',
#                                  data = results)
# 
# plot.gam(gam_minPPT_relPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#65432130', ylab = 'SHAP value')
# 
# 
# 
# ### meanPPT
# gam_meanPPT_relPolewarness <- gam(Mean_PPT_SHAP ~ s(relPolewardness, k = 4),
#                                 na.action = 'na.omit',
#                                 data = results)
# 
# plot.gam(gam_meanPPT_relPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF752D30', ylab = 'SHAP value')
# 
# ### maxPPT
# gam_maxPPT_relPolewarness <- gam(Max_PPT_SHAP ~ s(relPolewardness, k = 4),
#                                na.action = 'na.omit',
#                                data = results)
# 
# plot.gam(gam_maxPPT_relPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#00FF0030', ylab = 'SHAP value')
# 
# 
# # abs polewardness
# 
# 
# ### minPPT
# gam_minPPT_absPolewarness <- gam(Min_PPT_SHAP ~ s(absPolewardness, k = 4),
#                                na.action = 'na.omit',
#                                data = results)
# 
# plot.gam(gam_minPPT_absPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#65432130', ylab = 'SHAP value')
# 
# 
# 
# ### meanPPT
# gam_meanPPT_absPolewarness <- gam(Mean_PPT_SHAP ~ s(absPolewardness, k = 4),
#                                 na.action = 'na.omit',
#                                 data = results)
# 
# plot.gam(gam_meanPPT_absPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF752D30', ylab = 'SHAP value')
# 
# ### maxPPT
# gam_maxPPT_absPolewarness <- gam(Max_PPT_SHAP ~ s(absPolewardness, k = 4),
#                                na.action = 'na.omit',
#                                data = results)
# 
# plot.gam(gam_maxPPT_absPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#00FF0030', ylab = 'SHAP value')
# 
# 
# 
# # dist edge
# 
# ### minPPT
# gam_minPPT_distEdge <- gam(Min_PPT_SHAP ~ s(log(distEdge + 1), k = 4),
#                          na.action = 'na.omit',
#                          data = results)
# 
# plot.gam(gam_minPPT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#65432130', ylab = 'SHAP value')
# 
# 
# ### meanPPT
# gam_meanPPT_distEdge <- gam(Mean_PPT_SHAP ~ s(log(distEdge + 1), k = 4),
#                           na.action = 'na.omit',
#                           data = results)
# 
# plot.gam(gam_meanPPT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF752D30', ylab = 'SHAP value')
# 
# ### maxPPT
# gam_maxPPT_distEdge <- gam(Max_PPT_SHAP ~ s(log(distEdge + 1), k = 4),
#                          na.action = 'na.omit',
#                          data = results)
# 
# plot.gam(gam_maxPPT_distEdge, pages = 1, residuals = F, shade = T,
#          shade.col = '#00FF0030', ylab = 'SHAP value')
# 
# 
# # elevation
# 
# ### minPPT
# gam_minPPT_elev <- gam(Min_PPT_SHAP ~ s(elevation, k = 4),
#                      na.action = 'na.omit',
#                      data = results)
# 
# plot.gam(gam_minPPT_elev, pages = 1, residuals = F, shade = T,
#          shade.col = '#65432130', ylab = 'SHAP value')
# 
# 
# 
# ### meanPPT
# gam_meanPPT_elev <- gam(Mean_PPT_SHAP ~ s(elevation, k = 4),
#                       na.action = 'na.omit',
#                       data = results)
# 
# plot.gam(gam_meanPPT_elev, pages = 1, residuals = F, shade = T,
#          shade.col = '#FF752D30', ylab = 'SHAP value')
# 
# ### maxPPT
# gam_maxPPT_elev <- gam(Max_PPT_SHAP ~ s(elevation, k = 4),
#                      na.action = 'na.omit',
#                      data = results)
# 
# plot.gam(gam_maxPPT_elev, pages = 1, residuals = F, shade = T,
#          shade.col = '#00FF0030', ylab = 'SHAP value')
# 
# 
# 
# #########
# 
# install.packages("lme4")
# library(lme4)
# data("sleepstudy")
# 
# head(sleepstudy)
# 
# #example without random effects
# gam_model <- gam(Reaction ~ s(Days, k = 5), data = sleepstudy)
# 
# summary(gam_model)
# 
# plot(gam_model, pages = 1, residuals = F, shade = T, shade.col = '#0000FF30')
# 
# #example with random effects
# gam_model_re <- gam(Reaction ~ s(Days, k = 5) +
#                       s(Subject, bs = "re"), data = sleepstudy)
# 
# summary(gam_model_re)
# 
# plot(gam_model_re, pages = 1, residuals = F, shade = T, shade.col = '#0000FF30')
# 
# 
# # gam_model <- gam(Reaction ~ s(Days, k = 10), data = sleepstudy)
# # summary(gam_model)
# # 
# # plot(gam_model, pages = 1)
# 
# 
# ###
# plot.gam(gam_minT_relPolewarness, pages = 1, residuals = F, shade = T,
#          shade.col = '#0000FF30', ylab = 'SHAP value',
#          ylim = c(-0.06, 0.06),
#          cex.lab = 2, cex.axis = 1.5)
# 
