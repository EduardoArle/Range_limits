#load libraries
library(mgcv); library(itsadug)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20241208_All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20241210_GAMs'

#read results table
setwd(wd_tables)
results <- read.csv('20241210_Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,12)] <- c("absPolewardness","relPolewardness" )

#transform species names in factor
results$species <- as.factor(results$species)

###########.  GAM   ################

########################
######### PPT ##########
########################


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


summary(minT_absPolewarness_G)

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



# SHAP values X distance from edge (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G
minT_distEdge_G <- gam(Min_T_SHAP ~
                               s(distEdge, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(minT_distEdge_G)

setwd(wd_models)
saveRDS(minT_distEdge_G, 'minT_distEdge_G')

draw(minT_distEdge_G, page = 1)

plot.gam(minT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_distEdge_G, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS
minT_distEdge_GS <- gam(Min_T_SHAP ~
                                s(distEdge, k=4, m=2) 
                              + s(distEdge, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(minT_distEdge_GS)

setwd(wd_models)
saveRDS(minT_distEdge_GS, 'minT_distEdge_GS')

draw(minT_distEdge_GS, page = 1)


plot.gam(minT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(minT_distEdge_GS, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model S
minT_distEdge_S <- gam(Min_T_SHAP ~
                               s(distEdge, species, k=4, bs="fs", m=2),
                             data = results,
                             method="REML")

summary(minT_distEdge_S)

setwd(wd_models)
saveRDS(minT_distEdge_S, 'minT_distEdge_S')

draw(minT_distEdge_S, page = 1)

plot.gam(minT_distEdge_S, select = 1, residuals = F, shade = T,
         cex.lab = 2, cex.axis = 1.5) #save 800

### meanT

#model G
meanT_distEdge_G <- gam(Mean_T_SHAP ~
                                s(distEdge, k=4, bs="tp") 
                              + s(species, k=503, bs="re"),
                              data = results,
                              method="REML",
                              family="gaussian")


summary(meanT_distEdge_G)

setwd(wd_models)
saveRDS(meanT_distEdge_G, 'meanT_distEdge_G')

draw(meanT_distEdge_G, page = 1)

plot.gam(meanT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_distEdge_G, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS
meanT_distEdge_GS <- gam(Mean_T_SHAP ~
                                 s(distEdge, k=4, m=2) 
                               + s(distEdge, species, k=4, bs="fs", m=2),
                               data = results,
                               method="REML")


summary(meanT_distEdge_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_GS, 'meanT_distEdge_GS')

draw(meanT_distEdge_GS, page = 1)


plot.gam(meanT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(meanT_distEdge_GS, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model S
meanT_distEdge_S <- gam(Mean_T_SHAP ~
                                s(distEdge, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")

summary(meanT_distEdge_S)

setwd(wd_models)
saveRDS(meanT_distEdge_S, 'meanT_distEdge_S')

draw(meanT_distEdge_S, page = 1)

plot.gam(meanT_distEdge_S, select = 1, residuals = F, shade = T,
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxT

#model G
maxT_distEdge_G <- gam(Max_T_SHAP ~
                               s(distEdge, k=4, bs="tp") 
                             + s(species, k=503, bs="re"),
                             data = results,
                             method="REML",
                             family="gaussian")


summary(maxT_distEdge_G)

setwd(wd_models)
saveRDS(maxT_distEdge_G, 'maxT_distEdge_G')

draw(maxT_distEdge_G, page = 1)

plot.gam(maxT_distEdge_G, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_distEdge_G, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800

#model GS
maxT_distEdge_GS <- gam(Max_T_SHAP ~
                                s(distEdge, k=4, m=2) 
                              + s(distEdge, species, k=4, bs="fs", m=2),
                              data = results,
                              method="REML")


summary(maxT_distEdge_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_GS, 'maxT_distEdge_GS')

draw(maxT_distEdge_GS, page = 1)


plot.gam(maxT_distEdge_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

plot.gam(maxT_distEdge_GS, select = 2, residuals = F, 
         cex.lab = 2, cex.axis = 1.5)  #save 800


#model S
maxT_distEdge_S <- gam(Max_T_SHAP ~
                               s(distEdge, species, k=4, bs="fs", m=2),
                             data = results,
                             method="REML")

summary(maxT_distEdge_S)

setwd(wd_models)
saveRDS(maxT_distEdge_S, 'maxT_distEdge_S')

draw(maxT_distEdge_S, page = 1)

plot.gam(maxT_distEdge_S, select = 1, residuals = F, shade = T,
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
