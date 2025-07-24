#load libraries
library(mgcv); library(itsadug); library(gratia)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/20250504_All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/20250504_GAMs/Models'
wd_sig <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Results/Significance_GAMs'

#read results table
setwd(wd_tables)
results <- read.csv('20250504_Results_all_sps.csv')

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

results_meanT <- results[abs(results$Cor_vars_meanT) <= 0.7,]
results_meanT <- results_meanT[complete.cases(results_meanT$Cor_vars_meanT),]

results_maxT <- results[abs(results$Cor_vars_maxT) <= 0.7,]
results_maxT <- results_maxT[complete.cases(results_maxT$Cor_vars_maxT),]

#make a version of the results with absolute shap values
results_minT_abs <- results_minT
results_meanT_abs <- results_meanT
results_maxT_abs <- results_maxT

results_minT_abs$avg_Min_T_SHAP <- abs(results_minT_abs$avg_Min_T_SHAP)
results_meanT_abs$avg_Mean_T_SHAP <- abs(results_meanT_abs$avg_Mean_T_SHAP)
results_maxT_abs$avg_Max_T_SHAP <- abs(results_maxT_abs$avg_Max_T_SHAP)

#make a version of the results keeping only points up to 250km from edges
results_minT_abs_250 <- results_minT_abs[results_minT_abs$distEdge <= 250,]
results_meanT_abs_250 <- results_meanT_abs[results_meanT_abs$distEdge <= 250,]
results_maxT_abs_250 <- results_maxT_abs[results_maxT_abs$distEdge <= 250,]

#count species in each table (to be used in the latitudinal analyses)
n_sps_minT <- length(unique(results_minT$species))
n_sps_meanT <- length(unique(results_meanT$species))
n_sps_maxT <- length(unique(results_maxT$species))

#count species in each modified table (to be used in the centre vs edge analyses)
n_sps_minT_mod <- length(unique(results_minT_abs_250$species))
n_sps_meanT_mod <- length(unique(results_meanT_abs_250$species))
n_sps_maxT_mod <- length(unique(results_maxT_abs_250$species))



# SHAP values X distance from edge (GAM)


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


#model G (250)
minT_distEdge_250_G <- gam(avg_Min_T_SHAP ~
                           s(distEdge, k=4, bs="tp") 
                         + s(species, k = n_sps_minT_mod, bs="re"),
                           data = results_minT_abs_250,
                           method="REML",
                           family="gaussian")


summary(minT_distEdge_250_G)


#AIC
AIC_minT_distEdge_250_G <- AIC(minT_distEdge_250_G)

setwd(wd_models)
saveRDS(minT_distEdge_250_G, 'minT_distEdge_250_G')


plot.gam(minT_distEdge_250_G, select = 1, residuals = F, shade = T,
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
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250)
minT_distEdge_250_GS <- gam(avg_Min_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_minT_abs_250,
                        method="REML")


summary(minT_distEdge_250_GS)

#AIC
AIC_minT_distEdge_250_GS <- AIC(minT_distEdge_250_GS)

setwd(wd_models)
saveRDS(minT_distEdge_250_GS, 'minT_distEdge_250_GS')


plot.gam(minT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#0000FF30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(minT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minT_distEdge_250_GS.csv', row.names = F)



### meanT

#model G
meanT_distEdge_G <- gam(avg_Mean_T_SHAP ~
                          s(distEdge, k=4, bs="tp") 
                        + s(species, k = n_sps_meanT, bs="re"),
                        data = results_meanT_abs,
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


#model G (250)
meanT_distEdge_250_G <- gam(avg_Mean_T_SHAP ~
                             s(distEdge, k=4, bs="tp") 
                           + s(species, k = n_sps_meanT_mod, bs="re"),
                           data = results_meanT_abs_250,
                           method="REML",
                           family="gaussian")


summary(meanT_distEdge_250_G)


#AIC
AIC_meanT_distEdge_250_G <- AIC(meanT_distEdge_250_G)

setwd(wd_models)
saveRDS(meanT_distEdge_250_G, 'meanT_distEdge_250_G')


plot.gam(meanT_distEdge_250_G, select = 1, residuals = F, shade = T,
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
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250)
meanT_distEdge_250_GS <- gam(avg_Mean_T_SHAP ~
                           s(distEdge, k=4, m=2) 
                         + s(distEdge, species, k=4, bs="fs", m=2),
                         data = results_meanT_abs_250,
                         method="REML")


summary(meanT_distEdge_250_GS)

#AIC
AIC_meanT_distEdge_250_GS <- AIC(meanT_distEdge_250_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_250_GS, 'meanT_distEdge_250_GS')


plot.gam(meanT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#80008030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(meanT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_distEdge_250_GS.csv', row.names = F)


### maxT

#model G
maxT_distEdge_G <- gam(avg_Max_T_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_maxT, bs="re"),
                       data = results_maxT_abs,
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


#model G (250)
maxT_distEdge_250_G <- gam(avg_Max_T_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_maxT_mod, bs="re"),
                       data = results_maxT_abs_250,
                       method="REML",
                       family="gaussian")


summary(maxT_distEdge_250_G)

#AIC
AIC_maxT_distEdge_250_G <- AIC(maxT_distEdge_250_G)

setwd(wd_models)
saveRDS(maxT_distEdge_250_G, 'maxT_distEdge_250_G')


plot.gam(maxT_distEdge_250_G, select = 1, residuals = F, shade = T,
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
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250)
maxT_distEdge_250_GS <- gam(avg_Max_T_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_maxT_abs_250,
                        method="REML")


summary(maxT_distEdge_250_GS)

#AIC
AIC_maxT_distEdge_250_GS <- AIC(maxT_distEdge_250_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_250_GS, 'maxT_distEdge_250_GS')


plot.gam(maxT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#FF000030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(maxT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxT_distEdge_250_GS.csv', row.names = F)



# SHAP values X relative polewardness (GAM)

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


#test stats for significant change points from First Derivative of GAM smooth

#calculate derivatives
deriv <- derivatives(minT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minT_relPol_GS.csv', row.names = F)

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

#calculate derivatives
deriv <- derivatives(meanT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_relPol_GS.csv', row.names = F)


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

#calculate derivatives
deriv <- derivatives(maxT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxT_relPol_GS.csv', row.names = F)


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


#calculate derivatives
deriv <- derivatives(minT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minT_absPol_GS.csv', row.names = F)




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
         ylim = c(-0.5, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

#calculate derivatives
deriv <- derivatives(meanT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_absPol_GS.csv', row.names = F)



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


#calculate derivatives
deriv <- derivatives(maxT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxT_absPol_GS.csv', row.names = F)




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
results_minPPT <- results_minPPT[complete.cases(results_minPPT$Cor_vars_minPPT),]

results_meanPPT <- results[abs(results$Cor_vars_meanPPT) <= 0.7,]
results_meanPPT <- results_meanPPT[complete.cases(results_meanPPT$Cor_vars_meanPPT),]

results_maxPPT <- results[abs(results$Cor_vars_maxPPT) <= 0.7,]
results_maxPPT <- results_maxPPT[complete.cases(results_maxPPT$Cor_vars_maxPPT),]

#make a version of the results with absolute shap values
results_minPPT_abs <- results_minPPT
results_meanPPT_abs <- results_meanPPT
results_maxPPT_abs <- results_maxPPT

results_minPPT_abs$avg_Min_PPT_SHAP <- abs(results_minPPT_abs$avg_Min_PPT_SHAP)
results_meanPPT_abs$avg_Mean_PPT_SHAP <- abs(results_meanPPT_abs$avg_Mean_PPT_SHAP)
results_maxPPT_abs$avg_Max_PPT_SHAP <- abs(results_maxPPT_abs$avg_Max_PPT_SHAP)

#make a version of the results keeping only points up to 250km from edges
results_minPPT_abs_250 <- results_minPPT_abs[results_minPPT_abs$distEdge <= 250,]
results_meanPPT_abs_250 <- results_meanPPT_abs[results_meanPPT_abs$distEdge <= 250,]
results_maxPPT_abs_250 <- results_maxPPT_abs[results_maxPPT_abs$distEdge <= 250,]

#count species in each table (to be used in the latitudinal analyses)
n_sps_minPPT <- length(unique(results_minPPT$species))
n_sps_meanPPT <- length(unique(results_meanPPT$species))
n_sps_maxPPT <- length(unique(results_maxPPT$species))

#count species in each modified table (to be used in the centre vs edge analyses)
n_sps_minPPT_mod <- length(unique(results_minPPT_abs_250$species))
n_sps_meanPPT_mod <- length(unique(results_meanPPT_abs_250$species))
n_sps_maxPPT_mod <- length(unique(results_maxPPT_abs_250$species))



# SHAP values X distance from edge (GAM)

#set par for plotting
par(mar = c(6,6,6,6), pty='m')

### minPPT

#model G
minPPT_distEdge_G <- gam(avg_Min_PPT_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_minPPT, bs="re"),
                       data = results_minPPT_abs,
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


#model G (250)
minPPT_distEdge_250_G <- gam(avg_Min_PPT_SHAP ~
                             s(distEdge, k=4, bs="tp") 
                           + s(species, k = n_sps_minPPT_mod, bs="re"),
                           data = results_minPPT_abs_250,
                           method="REML",
                           family="gaussian")


summary(minPPT_distEdge_250_G)


#AIC
AIC_minPPT_distEdge_250_G <- AIC(minPPT_distEdge_250_G)

setwd(wd_models)
saveRDS(minPPT_distEdge_250_G, 'minPPT_distEdge_250_G')


plot.gam(minPPT_distEdge_250_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
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


#model GS (250)
minPPT_distEdge_250_GS <- gam(avg_Min_PPT_SHAP ~
                              s(distEdge, k=4, m=2) 
                            + s(distEdge, species, k=4, bs="fs", m=2),
                            data = results_minPPT_abs_250,
                            method="REML")


summary(minPPT_distEdge_250_GS)

#AIC
AIC_minPPT_distEdge_250_GS <- AIC(minPPT_distEdge_250_GS)

setwd(wd_models)
saveRDS(minPPT_distEdge_250_GS, 'minPPT_distEdge_250_GS')


plot.gam(minPPT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



### meanPPT

#model G
meanPPT_distEdge_G <- gam(avg_Mean_PPT_SHAP ~
                          s(distEdge, k=4, bs="tp") 
                        + s(species, k = n_sps_meanPPT, bs="re"),
                        data = results_meanPPT_abs,
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


#model G (250)
meanPPT_distEdge_250_G <- gam(avg_Mean_PPT_SHAP ~
                               s(distEdge, k=4, bs="tp") 
                             + s(species, k = n_sps_meanPPT_mod, bs="re"),
                             data = results_meanPPT_abs_250,
                             method="REML",
                             family="gaussian")


summary(meanPPT_distEdge_250_G)


#AIC
AIC_meanPPT_distEdge_250_G <- AIC(meanPPT_distEdge_250_G)

setwd(wd_models)
saveRDS(meanPPT_distEdge_250_G, 'meanPPT_distEdge_250_G')


plot.gam(meanPPT_distEdge_250_G, select = 1, residuals = F, shade = T,
         shade.col = '#fc8d5930', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
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


#model GS (250)
meanPPT_distEdge_250_GS <- gam(avg_Mean_PPT_SHAP ~
                             s(distEdge, k=4, m=2) 
                           + s(distEdge, species, k=4, bs="fs", m=2),
                           data = results_meanPPT_abs_250,
                           method="REML")


summary(meanPPT_distEdge_250_GS)

#AIC
AIC_meanPPT_distEdge_250_GS <- AIC(meanPPT_distEdge_250_GS)

setwd(wd_models)
saveRDS(meanPPT_distEdge_250_GS, 'meanPPT_distEdge_250_GS')


plot.gam(meanPPT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#8c510a30', ylab = 'SHAP value',
         ylim = c(-0.15, 0.05),
         cex.lab = 2, cex.axis = 1.5) #save 800



### maxPPT

#model G
maxPPT_distEdge_G <- gam(avg_Max_PPT_SHAP ~
                         s(distEdge, k=4, bs="tp") 
                       + s(species, k = n_sps_maxPPT, bs="re"),
                       data = results_maxPPT_abs,
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


#model G (250)
maxPPT_distEdge_250_G <- gam(avg_Max_PPT_SHAP ~
                                s(distEdge, k=4, bs="tp") 
                              + s(species, k = n_sps_maxPPT_mod, bs="re"),
                              data = results_maxPPT_abs_250,
                              method="REML",
                              family="gaussian")


summary(maxPPT_distEdge_250_G)


#AIC
AIC_maxPPT_distEdge_250_G <- AIC(maxPPT_distEdge_250_G)

setwd(wd_models)
saveRDS(maxPPT_distEdge_250_G, 'maxPPT_distEdge_250_G')


plot.gam(maxPPT_distEdge_250_G, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS
maxPPT_distEdge_GS <- gam(avg_Max_PPT_SHAP ~
                          s(distEdge, k=4, m=2) 
                        + s(distEdge, species, k=4, bs="fs", m=2),
                        data = results_maxPPT_abs,
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


#model GS (250)
maxPPT_distEdge_250_GS <- gam(avg_Max_PPT_SHAP ~
                                 s(distEdge, k=4, m=2) 
                               + s(distEdge, species, k=4, bs="fs", m=2),
                               data = results_maxPPT_abs_250,
                               method="REML")


summary(maxPPT_distEdge_250_GS)

#AIC
AIC_maxPPT_distEdge_250_GS <- AIC(maxPPT_distEdge_250_GS)

setwd(wd_models)
saveRDS(maxPPT_distEdge_250_GS, 'maxPPT_distEdge_250_GS')


plot.gam(maxPPT_distEdge_250_GS, select = 1, residuals = F, shade = T,
         shade.col = '#00FF0030', ylab = 'SHAP value',
         ylim = c(-0.15, 0.05),
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


plot.gam(meanPPT_elev_GS, select = 1, residuals = F, shade = T,
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







##############. meanT_relPol


meanT_relPol_GS <- readRDS('meanT_relPol_GS')

#estimate first derivatives

smooths(meanT_relPol_GS)


deriv <- derivatives(meanT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 200)  # number of points to evaluate

draw(deriv)

names(deriv)
head(deriv)


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

# Add significance to the derivatives data frame
deriv$sig <- sig_points

df <- deriv

str(deriv)

# Plot the derivative with its confidence interval
plot(deriv$relPolewardness, deriv$.derivative, type = "l",
     ylim = range(c(deriv$.lower_ci, deriv$.upper_ci), na.rm = TRUE),
     xlab = "Relative Polewardness", ylab = "First derivative",
     col = "blue", lwd = 2)

# Add confidence interval ribbon
polygon(c(deriv$relPolewardness, rev(deriv$relPolewardness)),
        c(deriv$.upper_ci, rev(deriv$.lower_ci)),
        col = rgb(0, 0, 1, 0.2), border = NA)

# Re-draw the derivative line (so it’s on top)
lines(deriv$relPolewardness, deriv$.derivative, col = "blue", lwd = 2)

# Add significant points in red
points(deriv$relPolewardness[deriv$sig],
       deriv$.derivative[deriv$sig],
       col = "red", pch = 19)




##############. maxT_relPol


maxT_relPol_GS <- readRDS('maxT_relPol_GS')

#estimate first derivatives

smooths(maxT_relPol_GS)


deriv <- derivatives(maxT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 200)  # number of points to evaluate

draw(deriv)

names(deriv)
head(deriv)


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

# Add significance to the derivatives data frame
deriv$sig <- sig_points

df <- deriv

str(deriv)

# Plot the derivative with its confidence interval
plot(deriv$relPolewardness, deriv$.derivative, type = "l",
     ylim = range(c(deriv$.lower_ci, deriv$.upper_ci), na.rm = TRUE),
     xlab = "Relative Polewardness", ylab = "First derivative",
     col = "blue", lwd = 2)

# Add confidence interval ribbon
polygon(c(deriv$relPolewardness, rev(deriv$relPolewardness)),
        c(deriv$.upper_ci, rev(deriv$.lower_ci)),
        col = rgb(0, 0, 1, 0.2), border = NA)

# Re-draw the derivative line (so it’s on top)
lines(deriv$relPolewardness, deriv$.derivative, col = "blue", lwd = 2)

# Add significant points in red
points(deriv$relPolewardness[deriv$sig],
       deriv$.derivative[deriv$sig],
       col = "red", pch = 19)




