#load libraries
library(mgcv); library(itsadug); library(gratia)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Models'
wd_sig <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Significance_GAMs'

#read results table
setwd(wd_tables)
results <- read.csv('20260126_Results_all_sps.csv')

#transform species names in factor
results$species <- as.factor(results$sps)

#select only presences (for now)
results <- results[results$Occurrence == 1,]

#thin results per max 1000 points per species
X <- 1000
set.seed(1)
results_thin <- do.call(rbind,
                        lapply(split(results, results$species),function(d){
                          if(nrow(d) > X) d[sort(sample.int(nrow(d), X)), ]
                          else d}))
row.names(results_thin) <- NULL


########################
########## T ###########
########################

#select only species that had lower correl between vars

res_minT <- results_thin[abs(results_thin$Cor_vars_minT) <= 0.7,]
res_minT <- res_minT[complete.cases(res_minT$Cor_vars_minT),]

res_meanT <- results[abs(results_thin$Cor_vars_meanT) <= 0.7,]
res_meanT <- res_meanT[complete.cases(res_meanT$Cor_vars_meanT),]

res_maxT <- results[abs(results_thin$Cor_vars_maxT) <= 0.7,]
res_maxT <- res_maxT[complete.cases(res_maxT$Cor_vars_maxT),]


#make a version of the results with absolute shap values
res_minT_abs <- res_minT
res_meanT_abs <- res_meanT
res_maxT_abs <- res_maxT

res_minT_abs$avg_Min_T_SHAP <- abs(res_minT_abs$avg_Min_T_SHAP)
res_meanT_abs$avg_Mean_T_SHAP <- abs(res_meanT_abs$avg_Mean_T_SHAP)
res_maxT_abs$avg_Max_T_SHAP <- abs(res_maxT_abs$avg_Max_T_SHAP)

#make a version of the results keeping only points up to 250km from edges
res_minT_abs_250 <- res_minT_abs[res_minT_abs$distEdge <= 250,]
res_meanT_abs_250 <- res_meanT_abs[res_meanT_abs$distEdge <= 250,]
res_maxT_abs_250 <- res_maxT_abs[res_maxT_abs$distEdge <= 250,]

#count species in each table (to be used in the latitudinal analyses)
n_sps_minT <- length(unique(res_minT$species))
n_sps_meanT <- length(unique(res_meanT$species))
n_sps_maxT <- length(unique(res_maxT$species))

#count species in each modified table (to be used in the centre vs edge analyses)
n_sps_minT_mod <- length(unique(res_minT_abs_250$species))
n_sps_meanT_mod <- length(unique(res_meanT_abs_250$species))
n_sps_maxT_mod <- length(unique(res_maxT_abs_250$species))


# SHAP values X distance from edge (GAM)

### minT

#model G (presence)
minT_distEdge_G <- bam(avg_Min_T_SHAP ~
                         s(distEdge, k = 4, bs = "tp")
                       + s(species, k = n_sps_minT, bs = "re"),
                       data = res_minT_abs,
                       method = "fREML",
                       discrete = TRUE,
                       family = gaussian())


summary(minT_distEdge_G)
gam.check(minT_distEdge_G)


#AIC
AIC_minT_distEdge_G <- AIC(minT_distEdge_G)

setwd(wd_models)
saveRDS(minT_distEdge_G, 'minT_distEdge_G')


plot.gam(minT_distEdge_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
minT_distEdge_250_G <- bam(avg_Min_T_SHAP ~
                             s(distEdge, k = 4, bs = "tp")
                           + s(species, k = n_sps_minT_mod, bs = "re"),
                           data = res_minT_abs_250,
                           method = "fREML",
                           discrete = TRUE,
                           family = gaussian())


summary(minT_distEdge_250_G)
gam.check(minT_distEdge_250_G)


#AIC
AIC_minT_distEdge_250_G <- AIC(minT_distEdge_250_G)

setwd(wd_models)
saveRDS(minT_distEdge_250_G, 'minT_distEdge_250_G')


plot.gam(minT_distEdge_250_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minT_distEdge_GS <- bam(avg_Min_T_SHAP ~
                          s(distEdge, k = 4, m = 2)
                        + s(distEdge, species, k = 4, bs = "fs", m = 2),
                        data = res_minT_abs,
                        method = "fREML", discrete = TRUE)

summary(minT_distEdge_GS)
gam.check(minT_distEdge_GS)

#AIC
AIC_minT_distEdge_GS <- AIC(minT_distEdge_GS)

setwd(wd_models)
saveRDS(minT_distEdge_GS, 'minT_distEdge_GS')

plot.gam(minT_distEdge_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250 presence)
minT_distEdge_250_GS <- bam(avg_Min_T_SHAP ~
                              s(distEdge, k = 4, m = 2)
                            + s(distEdge, species, k = 4, bs = "fs", m = 2),
                            data = res_minT_abs_250,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(minT_distEdge_250_GS)
gam.check(minT_distEdge_250_GS)

#AIC
AIC_minT_distEdge_250_GS <- AIC(minT_distEdge_250_GS)

setwd(wd_models)
saveRDS(minT_distEdge_250_GS, 'minT_distEdge_250_GS')


plot.gam(minT_distEdge_250_GS, select = 1, residuals = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.15),
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
 
#model G (presence)
meanT_distEdge_G <- bam(avg_Mean_T_SHAP ~
                         s(distEdge, k = 4, bs = "tp")
                       + s(species, k = n_sps_meanT, bs = "re"),
                       data = res_meanT_abs,
                       method = "fREML",
                       discrete = TRUE,
                       family = gaussian())


summary(meanT_distEdge_G)

#AIC
AIC_meanT_distEdge_G <- AIC(meanT_distEdge_G)

setwd(wd_models)
saveRDS(meanT_distEdge_G, 'meanT_distEdge_G')


plot.gam(meanT_distEdge_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.20, 0.30),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
meanT_distEdge_250_G <- bam(avg_Mean_T_SHAP ~
                             s(distEdge, k = 4, bs = "tp")
                           + s(species, k = n_sps_meanT_mod, bs = "re"),
                           data = res_meanT_abs_250,
                           method = "fREML",
                           discrete = TRUE,
                           family = gaussian())


summary(meanT_distEdge_250_G)
gam.check(meanT_distEdge_G)

#AIC
AIC_meanT_distEdge_250_G <- AIC(meanT_distEdge_250_G)

setwd(wd_models)
saveRDS(meanT_distEdge_250_G, 'meanT_distEdge_250_G')


plot.gam(meanT_distEdge_250_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence)
meanT_distEdge_GS <- bam(avg_Mean_T_SHAP ~
                          s(distEdge, k = 4, m = 2)
                        + s(distEdge, species, k = 4, bs = "fs", m = 2),
                        data = res_meanT_abs,
                        method = "fREML", discrete = TRUE)


summary(meanT_distEdge_GS)
gam.check(meanT_distEdge_GS)

#AIC
AIC_meanT_distEdge_GS <- AIC(meanT_distEdge_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_GS, 'meanT_distEdge_GS')


plot.gam(meanT_distEdge_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
meanT_distEdge_250_GS <- bam(avg_Mean_T_SHAP ~
                              s(distEdge, k = 4, m = 2)
                            + s(distEdge, species, k = 4, bs = "fs", m = 2),
                            data = res_meanT_abs_250,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(meanT_distEdge_250_GS)
gam.check(meanT_distEdge_250_GS)

#AIC
AIC_meanT_distEdge_250_GS <- AIC(meanT_distEdge_250_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_250_GS, 'meanT_distEdge_250_GS')


plot.gam(meanT_distEdge_250_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
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

#model G (presence) - fast bam
maxT_distEdge_G <- bam(avg_Max_T_SHAP ~
                          s(distEdge, k = 4, bs = "tp")
                        + s(species, k = n_sps_maxT, bs = "re"),
                        data = res_maxT_abs,
                        method = "fREML",
                        discrete = TRUE,
                        family = gaussian())


summary(maxT_distEdge_G)
gam.check(maxT_distEdge_G)

#AIC
AIC_maxT_distEdge_G <- AIC(maxT_distEdge_G)

setwd(wd_models)
saveRDS(maxT_distEdge_G, 'maxT_distEdge_G')


plot.gam(maxT_distEdge_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (presence 250 km) 
maxT_distEdge_250_G <- bam(avg_Max_T_SHAP ~
                              s(distEdge, k = 4, bs = "tp")
                            + s(species, k = n_sps_maxT_mod, bs = "re"),
                            data = res_maxT_abs_250,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(maxT_distEdge_250_G)
gam.check(maxT_distEdge_250_G)

#AIC
AIC_maxT_distEdge_250_G <- AIC(maxT_distEdge_250_G)

setwd(wd_models)
saveRDS(maxT_distEdge_250_G, 'maxT_distEdge_250_G')


plot.gam(maxT_distEdge_250_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS

#model GS (presence) 
maxT_distEdge_GS <- bam(avg_Max_T_SHAP ~
                           s(distEdge, k = 4, m = 2)
                         + s(distEdge, species, k = 4, bs = "fs", m = 2),
                         data = res_maxT_abs,
                         method = "fREML", discrete = TRUE)

summary(maxT_distEdge_GS)
gam.check(maxT_distEdge_GS)

#AIC
AIC_maxT_distEdge_GS <- AIC(maxT_distEdge_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_GS, 'maxT_distEdge_GS')


plot.gam(maxT_distEdge_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.5, 0.5),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
maxT_distEdge_250_GS <- bam(avg_Max_T_SHAP ~
                               s(distEdge, k = 4, m = 2)
                             + s(distEdge, species, k = 4, bs = "fs", m = 2),
                             data = res_maxT_abs_250,
                             method = "fREML",
                             discrete = TRUE,
                             family = gaussian())


summary(maxT_distEdge_250_GS)
gam.check(maxT_distEdge_250_GS)

#AIC
AIC_maxT_distEdge_250_GS <- AIC(maxT_distEdge_250_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_250_GS, 'maxT_distEdge_250_GS')


plot.gam(maxT_distEdge_250_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
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




######### Select only sps not too limited by the ocean  (max 20%) #########

res_minT_abs_20 <- res_minT_abs[
  as.numeric(gsub("%", "", res_minT_abs$perc_ocean)) <= 20, ]

res_minT_abs_250_20 <- res_minT_abs_250[
  as.numeric(gsub("%", "", res_minT_abs_250$perc_ocean)) <= 20, ]

res_meanT_abs_20 <- res_meanT_abs[
  as.numeric(gsub("%", "", res_meanT_abs$perc_ocean)) <= 20, ]

res_meanT_abs_250_20 <- res_meanT_abs_250[
  as.numeric(gsub("%", "", res_meanT_abs_250$perc_ocean)) <= 20, ]

res_maxT_abs_20 <- res_maxT_abs[
  as.numeric(gsub("%", "", res_maxT_abs$perc_ocean)) <= 20, ]

res_maxT_abs_250_20 <- res_maxT_abs_250[
  as.numeric(gsub("%", "", res_maxT_abs_250$perc_ocean)) <= 20, ]

# SHAP values X distance from edge (GAM)

### minT

#model G (presence)
minT_distEdge_20_G <- bam(avg_Min_T_SHAP ~
                         s(distEdge, k = 4, bs = "tp")
                       + s(species, k = n_sps_minT, bs = "re"),
                       data = res_minT_abs_20,
                       method = "fREML",
                       discrete = TRUE,
                       family = gaussian())


summary(minT_distEdge_20_G)
gam.check(minT_distEdge_20_G)


#AIC
AIC_minT_distEdge_20_G <- AIC(minT_distEdge_20_G)

setwd(wd_models)
saveRDS(minT_distEdge_20_G, 'minT_distEdge_20_G')


plot.gam(minT_distEdge_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
minT_distEdge_250_20_G <- bam(avg_Min_T_SHAP ~
                             s(distEdge, k = 4, bs = "tp")
                           + s(species, k = n_sps_minT_mod, bs = "re"),
                           data = res_minT_abs_250_20,
                           method = "fREML",
                           discrete = TRUE,
                           family = gaussian())


summary(minT_distEdge_250_20_G)
gam.check(minT_distEdge_250_20_G)


#AIC
AIC_minT_distEdge_250_20_G <- AIC(minT_distEdge_250_20_G)

setwd(wd_models)
saveRDS(minT_distEdge_250_20_G, 'minT_distEdge_250_20_G')


plot.gam(minT_distEdge_250_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minT_distEdge_20_GS <- bam(avg_Min_T_SHAP ~
                          s(distEdge, k = 4, m = 2)
                        + s(distEdge, species, k = 4, bs = "fs", m = 2),
                        data = res_minT_abs_20,
                        method = "fREML", discrete = TRUE)

summary(minT_distEdge_20_GS)
gam.check(minT_distEdge_20_GS)

#AIC
AIC_minT_distEdge_20_GS <- AIC(minT_distEdge_20_GS)

setwd(wd_models)
saveRDS(minT_distEdge_20_GS, 'minT_distEdge_20_GS')

plot.gam(minT_distEdge_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.4, 0.4),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250 presence)
minT_distEdge_250_20_GS <- bam(avg_Min_T_SHAP ~
                              s(distEdge, k = 4, m = 2)
                            + s(distEdge, species, k = 4, bs = "fs", m = 2),
                            data = res_minT_abs_250_20,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(minT_distEdge_250_20_GS)
gam.check(minT_distEdge_250_20_GS)

#AIC
AIC_minT_distEdge_250_20_GS <- AIC(minT_distEdge_250_20_GS)

setwd(wd_models)
saveRDS(minT_distEdge_250_20_GS, 'minT_distEdge_250_20_GS')


plot.gam(minT_distEdge_250_20_GS, select = 1, residuals = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.15),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(minT_distEdge_250_20_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minT_distEdge_250_20_GS.csv', row.names = F)



### meanT

#model G (presence)
meanT_distEdge_20_G <- bam(avg_Mean_T_SHAP ~
                          s(distEdge, k = 4, bs = "tp")
                        + s(species, k = n_sps_meanT, bs = "re"),
                        data = res_meanT_abs_20,
                        method = "fREML",
                        discrete = TRUE,
                        family = gaussian())


summary(meanT_distEdge_20_G)

#AIC
AIC_meanT_distEdge_20_G <- AIC(meanT_distEdge_20_G)

setwd(wd_models)
saveRDS(meanT_distEdge_20_G, 'meanT_distEdge_20_G')


plot.gam(meanT_distEdge_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.20, 0.30),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
meanT_distEdge_250_20_G <- bam(avg_Mean_T_SHAP ~
                              s(distEdge, k = 4, bs = "tp")
                            + s(species, k = n_sps_meanT_mod, bs = "re"),
                            data = res_meanT_abs_250_20,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(meanT_distEdge_250_20_G)
gam.check(meanT_distEdge_20_G)

#AIC
AIC_meanT_distEdge_250_20_G <- AIC(meanT_distEdge_250_20_G)

setwd(wd_models)
saveRDS(meanT_distEdge_250_20_G, 'meanT_distEdge_250_20_G')


plot.gam(meanT_distEdge_250_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence)
meanT_distEdge_20_GS <- bam(avg_Mean_T_SHAP ~
                           s(distEdge, k = 4, m = 2)
                         + s(distEdge, species, k = 4, bs = "fs", m = 2),
                         data = res_meanT_abs_20,
                         method = "fREML", discrete = TRUE)


summary(meanT_distEdge_20_GS)
gam.check(meanT_distEdge_20_GS)

#AIC
AIC_meanT_distEdge_20_GS <- AIC(meanT_distEdge_20_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_20_GS, 'meanT_distEdge_20_GS')


plot.gam(meanT_distEdge_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
meanT_distEdge_250_20_GS <- bam(avg_Mean_T_SHAP ~
                               s(distEdge, k = 4, m = 2)
                             + s(distEdge, species, k = 4, bs = "fs", m = 2),
                             data = res_meanT_abs_250_20,
                             method = "fREML",
                             discrete = TRUE,
                             family = gaussian())


summary(meanT_distEdge_250_20_GS)
gam.check(meanT_distEdge_250_20_GS)

#AIC
AIC_meanT_distEdge_250_20_GS <- AIC(meanT_distEdge_250_20_GS)

setwd(wd_models)
saveRDS(meanT_distEdge_250_20_GS, 'meanT_distEdge_250_20_GS')


plot.gam(meanT_distEdge_250_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.2, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(meanT_distEdge_250_20_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_distEdge_250_20_GS.csv', row.names = F)


### maxT

#model G (presence) - fast bam
maxT_distEdge_20_G <- bam(avg_Max_T_SHAP ~
                         s(distEdge, k = 4, bs = "tp")
                       + s(species, k = n_sps_maxT, bs = "re"),
                       data = res_maxT_abs_20,
                       method = "fREML",
                       discrete = TRUE,
                       family = gaussian())


summary(maxT_distEdge_20_G)
gam.check(maxT_distEdge_20_G)

#AIC
AIC_maxT_distEdge_20_G <- AIC(maxT_distEdge_20_G)

setwd(wd_models)
saveRDS(maxT_distEdge_20_G, 'maxT_distEdge_20_G')


plot.gam(maxT_distEdge_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (presence 250 km) 
maxT_distEdge_250_20_G <- bam(avg_Max_T_SHAP ~
                             s(distEdge, k = 4, bs = "tp")
                           + s(species, k = n_sps_maxT_mod, bs = "re"),
                           data = res_maxT_abs_250_20,
                           method = "fREML",
                           discrete = TRUE,
                           family = gaussian())


summary(maxT_distEdge_250_20_G)
gam.check(maxT_distEdge_250_20_G)

#AIC
AIC_maxT_distEdge_250_20_G <- AIC(maxT_distEdge_250_20_G)

setwd(wd_models)
saveRDS(maxT_distEdge_250_20_G, 'maxT_distEdge_250_20_G')


plot.gam(maxT_distEdge_250_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS

#model GS (presence) 
maxT_distEdge_20_GS <- bam(avg_Max_T_SHAP ~
                          s(distEdge, k = 4, m = 2)
                        + s(distEdge, species, k = 4, bs = "fs", m = 2),
                        data = res_maxT_abs_20,
                        method = "fREML", discrete = TRUE)

summary(maxT_distEdge_20_GS)
gam.check(maxT_distEdge_20_GS)

#AIC
AIC_maxT_distEdge_20_GS <- AIC(maxT_distEdge_20_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_20_GS, 'maxT_distEdge_20_GS')


plot.gam(maxT_distEdge_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.5, 0.5),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
maxT_distEdge_250_20_GS <- bam(avg_Max_T_SHAP ~
                              s(distEdge, k = 4, m = 2)
                            + s(distEdge, species, k = 4, bs = "fs", m = 2),
                            data = res_maxT_abs_250_20,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(maxT_distEdge_250_20_GS)
gam.check(maxT_distEdge_250_20_GS)

#AIC
AIC_maxT_distEdge_250_20_GS <- AIC(maxT_distEdge_250_20_GS)

setwd(wd_models)
saveRDS(maxT_distEdge_250_20_GS, 'maxT_distEdge_250_20_GS')


plot.gam(maxT_distEdge_250_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.15),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(maxT_distEdge_250_20_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxT_distEdge_250_20_GS.csv', row.names = F)




# SHAP values X relative polewardness (GAM)

### minT

#model G (presence) 
minT_relPol_G <- bam(avg_Min_T_SHAP ~
                       s(relPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_minT, bs = "re"),
                     data = res_minT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(minT_relPol_G)
gam.check(minT_relPol_G)

#AIC
AIC_minT_relPol_G <- AIC(minT_relPol_G)

setwd(wd_models)
saveRDS(minT_relPol_G, 'minT_relPol_G')


plot.gam(minT_relPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minT_relPol_GS <- bam(avg_Min_T_SHAP ~
                        s(relPolewardness, k = 4, m = 2)
                      + s(relPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_minT,
                      method = "fREML",
                      discrete = TRUE)


summary(minT_relPol_GS)
gam.check(minT_relPol_GS)

#AIC
AIC_minT_relPol_GS <- AIC(minT_relPol_GS)

setwd(wd_models)
saveRDS(minT_relPol_GS, 'minT_relPol_GS')


plot.gam(minT_relPol_GS, select = 1, residuals = F, shade = T,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.2),
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

#model G (presence) 
meanT_relPol_G <- bam(avg_Mean_T_SHAP ~
                       s(relPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_meanT, bs = "re"),
                     data = res_meanT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(meanT_relPol_G)
gam.check(meanT_relPol_G)

#AIC
AIC_meanT_relPol_G <- AIC(meanT_relPol_G)

setwd(wd_models)
saveRDS(meanT_relPol_G, 'meanT_relPol_G')


plot.gam(meanT_relPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
meanT_relPol_GS <- bam(avg_Mean_T_SHAP ~
                        s(relPolewardness, k = 4, m = 2)
                      + s(relPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_meanT,
                      method = "fREML",
                      discrete = TRUE)


summary(meanT_relPol_GS)
gam.check(meanT_relPol_GS)

#AIC
AIC_meanT_relPol_GS <- AIC(meanT_relPol_GS)

setwd(wd_models)
saveRDS(meanT_relPol_GS, 'meanT_relPol_GS')


plot.gam(meanT_relPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.25, 0.25),
         cex.lab = 2, cex.axis = 1.5) #save 800



#calculate derivatives
deriv <- derivatives(meanT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate

#identify points where derivative is significantly different from zero
sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanT_relPol_GS.csv', row.names = F)


### maxT

#model G (presence) 
maxT_relPol_G <- bam(avg_Max_T_SHAP ~
                        s(relPolewardness, k = 4, bs = "tp")
                      + s(species, k = n_sps_maxT, bs = "re"),
                      data = res_maxT,
                      method = "fREML",
                      discrete = TRUE,
                      family = gaussian())


summary(maxT_relPol_G)
gam.check(maxT_relPol_G)

#AIC
AIC_maxT_relPol_G <- AIC(maxT_relPol_G)

setwd(wd_models)
saveRDS(maxT_relPol_G, 'maxT_relPol_G')


plot.gam(maxT_relPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.06, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence)
maxT_relPol_GS <- bam(avg_Max_T_SHAP ~
                         s(relPolewardness, k = 4, m = 2)
                       + s(relPolewardness, species, k = 4, bs = "fs", m = 2),
                       data = res_maxT,
                       method = "fREML",
                       discrete = TRUE)


summary(maxT_relPol_GS)
gam.check(maxT_relPol_GS)

#AIC
AIC_maxT_relPol_GS <- AIC(maxT_relPol_GS)

setwd(wd_models)
saveRDS(maxT_relPol_GS, 'maxT_relPol_GS')


plot.gam(maxT_relPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
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

#model G (presence)
minT_absPol_G <- bam(avg_Min_T_SHAP ~
                       s(absPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_minT, bs = "re"),
                     data = res_minT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(minT_absPol_G)
gam.check(minT_absPol_G)

#AIC
AIC_minT_absPol_G <- AIC(minT_absPol_G)

setwd(wd_models)
saveRDS(minT_absPol_G, 'minT_absPol_G')


plot.gam(minT_absPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.5, 0.5),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minT_absPol_GS <- bam(avg_Min_T_SHAP ~
                        s(absPolewardness, k = 4, m = 2)
                      + s(absPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_minT,
                      method = "fREML",
                      discrete = TRUE)


summary(minT_absPol_GS)
gam.check(minT_absPol_GS)

#AIC
AIC_minT_absPol_GS <- AIC(minT_absPol_GS)

setwd(wd_models)
saveRDS(minT_absPol_GS, 'minT_absPol_GS')


plot.gam(minT_absPol_GS, select = 1, residuals = F,
         ylab = 'SHAP value',
         ylim = c(-0.8, 0.8),
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

#model G (presence)
meanT_absPol_G <- bam(avg_Mean_T_SHAP ~
                       s(absPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_meanT, bs = "re"),
                     data = res_meanT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(meanT_absPol_G)
gam.check(meanT_absPol_G)

#AIC
AIC_meanT_absPol_G <- AIC(meanT_absPol_G)

setwd(wd_models)
saveRDS(meanT_absPol_G, 'meanT_absPol_G')


plot.gam(meanT_absPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.8, 0.8),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
meanT_absPol_GS <- bam(avg_Mean_T_SHAP ~
                        s(absPolewardness, k = 4, m = 2)
                      + s(absPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_meanT,
                      method = "fREML",
                      discrete = TRUE)


summary(meanT_absPol_GS)
gam.check(meanT_absPol_GS)

#AIC
AIC_meanT_absPol_GS <- AIC(meanT_absPol_GS)

setwd(wd_models)
saveRDS(meanT_absPol_GS, 'meanT_absPol_GS')


plot.gam(meanT_absPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-4, 1.5),
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

#model G (presence)
maxT_absPol_G <- bam(avg_Max_T_SHAP ~
                        s(absPolewardness, k = 4, bs = "tp")
                      + s(species, k = n_sps_maxT, bs = "re"),
                      data = res_maxT,
                      method = "fREML",
                      discrete = TRUE,
                      family = gaussian())


summary(maxT_absPol_G)
gam.check(maxT_absPol_G)

#AIC
AIC_maxT_absPol_G <- AIC(maxT_absPol_G)

setwd(wd_models)
saveRDS(maxT_absPol_G, 'maxT_absPol_G')


plot.gam(maxT_absPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence)
maxT_absPol_GS <- bam(avg_Max_T_SHAP ~
                         s(absPolewardness, k = 4, m = 2)
                       + s(absPolewardness, species, k = 4, bs = "fs", m = 2),
                       data = res_maxT,
                       method = "fREML",
                       discrete = TRUE)


summary(maxT_absPol_GS)
gam.check(maxT_absPol_GS)

#AIC
AIC_maxT_absPol_GS <- AIC(maxT_absPol_GS)

setwd(wd_models)
saveRDS(maxT_absPol_GS, 'maxT_absPol_GS')


plot.gam(maxT_absPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-2.8, 1.3),
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






















########################
######### PPT ##########
########################

#select only species that had lower correl between vars

res_minPPT <- results_thin[abs(results_thin$Cor_vars_minPPT) <= 0.7,]
res_minPPT <- res_minPPT[complete.cases(res_minPPT$Cor_vars_minPPT),]

res_meanPPT <- results[abs(results_thin$Cor_vars_meanPPT) <= 0.7,]
res_meanPPT <- res_meanPPT[complete.cases(res_meanPPT$Cor_vars_meanPPT),]

res_maxPPT <- results[abs(results_thin$Cor_vars_maxPPT) <= 0.7,]
res_maxPPT <- res_maxPPT[complete.cases(res_maxPPT$Cor_vars_maxPPT),]


#make a version of the results with absolute shap values
res_minPPT_abs <- res_minPPT
res_meanPPT_abs <- res_meanPPT
res_maxPPT_abs <- res_maxPPT

res_minPPT_abs$avg_Min_PPT_SHAP <- abs(res_minPPT_abs$avg_Min_PPT_SHAP)
res_meanPPT_abs$avg_Mean_PPT_SHAP <- abs(res_meanPPT_abs$avg_Mean_PPT_SHAP)
res_maxPPT_abs$avg_Max_PPT_SHAP <- abs(res_maxPPT_abs$avg_Max_PPT_SHAP)

#make a version of the results keeping only points up to 250km from edges
res_minPPT_abs_250 <- res_minPPT_abs[res_minPPT_abs$distEdge <= 250,]
res_meanPPT_abs_250 <- res_meanPPT_abs[res_meanPPT_abs$distEdge <= 250,]
res_maxPPT_abs_250 <- res_maxPPT_abs[res_maxPPT_abs$distEdge <= 250,]

#count species in each table (to be used in the latitudinal analyses)
n_sps_minPPT <- length(unique(res_minPPT$species))
n_sps_meanPPT <- length(unique(res_meanPPT$species))
n_sps_maxPPT <- length(unique(res_maxPPT$species))

#count species in each modified table (to be used in the centre vs edge analyses)
n_sps_minPPT_mod <- length(unique(res_minPPT_abs_250$species))
n_sps_meanPPT_mod <- length(unique(res_meanPPT_abs_250$species))
n_sps_maxPPT_mod <- length(unique(res_maxPPT_abs_250$species))


# SHAP values X distance from edge (GAM)

### minPPT

#model G (presence)
minPPT_distEdge_G <- bam(avg_Min_PPT_SHAP ~
                         s(distEdge, k = 4, bs = "tp")
                       + s(species, k = n_sps_minPPT, bs = "re"),
                       data = res_minPPT_abs,
                       method = "fREML",
                       discrete = TRUE,
                       family = gaussian())


summary(minPPT_distEdge_G)
gam.check(minPPT_distEdge_G)


#AIC
AIC_minPPT_distEdge_G <- AIC(minPPT_distEdge_G)

setwd(wd_models)
saveRDS(minPPT_distEdge_G, 'minPPT_distEdge_G')


plot.gam(minPPT_distEdge_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
minPPT_distEdge_250_G <- bam(avg_Min_PPT_SHAP ~
                             s(distEdge, k = 4, bs = "tp")
                           + s(species, k = n_sps_minPPT_mod, bs = "re"),
                           data = res_minPPT_abs_250,
                           method = "fREML",
                           discrete = TRUE,
                           family = gaussian())


summary(minPPT_distEdge_250_G)
gam.check(minPPT_distEdge_250_G)


#AIC
AIC_minPPT_distEdge_250_G <- AIC(minPPT_distEdge_250_G)

setwd(wd_models)
saveRDS(minPPT_distEdge_250_G, 'minPPT_distEdge_250_G')


plot.gam(minPPT_distEdge_250_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minPPT_distEdge_GS <- bam(avg_Min_PPT_SHAP ~
                          s(distEdge, k = 4, m = 2)
                        + s(distEdge, species, k = 4, bs = "fs", m = 2),
                        data = res_minPPT_abs,
                        method = "fREML", discrete = TRUE,
                        family = gaussian())

summary(minPPT_distEdge_GS)
gam.check(minPPT_distEdge_GS)

#AIC
AIC_minPPT_distEdge_GS <- AIC(minPPT_distEdge_GS)

setwd(wd_models)
saveRDS(minPPT_distEdge_GS, 'minPPT_distEdge_GS')

plot.gam(minPPT_distEdge_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.2, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250 presence)
minPPT_distEdge_250_GS <- bam(avg_Min_PPT_SHAP ~
                              s(distEdge, k = 4, m = 2)
                            + s(distEdge, species, k = 4, bs = "fs", m = 2),
                            data = res_minPPT_abs_250,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(minPPT_distEdge_250_GS)
gam.check(minPPT_distEdge_250_GS)

#AIC
AIC_minPPT_distEdge_250_GS <- AIC(minPPT_distEdge_250_GS)

setwd(wd_models)
saveRDS(minPPT_distEdge_250_GS, 'minPPT_distEdge_250_GS')


plot.gam(minPPT_distEdge_250_GS, select = 1, residuals = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.15),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(minPPT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minPPT_distEdge_250_GS.csv', row.names = F)



### meanPPT

#model G (presence)
meanPPT_distEdge_G <- bam(avg_Mean_PPT_SHAP ~
                          s(distEdge, k = 4, bs = "tp")
                        + s(species, k = n_sps_meanPPT, bs = "re"),
                        data = res_meanPPT_abs,
                        method = "fREML",
                        discrete = TRUE,
                        family = gaussian())


summary(meanPPT_distEdge_G)

#AIC
AIC_meanPPT_distEdge_G <- AIC(meanPPT_distEdge_G)

setwd(wd_models)
saveRDS(meanPPT_distEdge_G, 'meanPPT_distEdge_G')


plot.gam(meanPPT_distEdge_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.20, 0.30),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
meanPPT_distEdge_250_G <- bam(avg_Mean_PPT_SHAP ~
                              s(distEdge, k = 4, bs = "tp")
                            + s(species, k = n_sps_meanPPT_mod, bs = "re"),
                            data = res_meanPPT_abs_250,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(meanPPT_distEdge_250_G)
gam.check(meanPPT_distEdge_G)

#AIC
AIC_meanPPT_distEdge_250_G <- AIC(meanPPT_distEdge_250_G)

setwd(wd_models)
saveRDS(meanPPT_distEdge_250_G, 'meanPPT_distEdge_250_G')


plot.gam(meanPPT_distEdge_250_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence 250 km)
meanPPT_distEdge_GS <- bam(avg_Mean_PPT_SHAP ~
                           s(distEdge, k = 4, m = 2)
                         + s(distEdge, species, k = 4, bs = "fs", m = 2),
                         data = res_meanPPT_abs,
                         method = "fREML", discrete = TRUE,
                         family = gaussian())


summary(meanPPT_distEdge_GS)
gam.check(meanPPT_distEdge_GS)

#AIC
AIC_meanPPT_distEdge_GS <- AIC(meanPPT_distEdge_GS)

setwd(wd_models)
saveRDS(meanPPT_distEdge_GS, 'meanPPT_distEdge_GS')


plot.gam(meanPPT_distEdge_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
meanPPT_distEdge_250_GS <- bam(avg_Mean_PPT_SHAP ~
                               s(distEdge, k = 4, m = 2)
                             + s(distEdge, species, k = 4, bs = "fs", m = 2),
                             data = res_meanPPT_abs_250,
                             method = "fREML",
                             discrete = TRUE,
                             family = gaussian())


summary(meanPPT_distEdge_250_GS)
gam.check(meanPPT_distEdge_250_GS)

#AIC
AIC_meanPPT_distEdge_250_GS <- AIC(meanPPT_distEdge_250_GS)

setwd(wd_models)
saveRDS(meanPPT_distEdge_250_GS, 'meanPPT_distEdge_250_GS')


plot.gam(meanPPT_distEdge_250_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.2, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(meanPPT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanPPT_distEdge_250_GS.csv', row.names = F)


### maxT

#model G (presence) - fast bam
maxPPT_distEdge_G <- bam(avg_Max_PPT_SHAP ~
                         s(distEdge, k = 4, bs = "tp")
                       + s(species, k = n_sps_maxPPT, bs = "re"),
                       data = res_maxPPT_abs,
                       method = "fREML",
                       discrete = TRUE,
                       family = gaussian())


summary(maxPPT_distEdge_G)
gam.check(maxPPT_distEdge_G)

#AIC
AIC_maxPPT_distEdge_G <- AIC(maxPPT_distEdge_G)

setwd(wd_models)
saveRDS(maxPPT_distEdge_G, 'maxPPT_distEdge_G')


plot.gam(maxPPT_distEdge_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (presence 250 km) 
maxPPT_distEdge_250_G <- bam(avg_Max_PPT_SHAP ~
                             s(distEdge, k = 4, bs = "tp")
                           + s(species, k = n_sps_maxPPT_mod, bs = "re"),
                           data = res_maxPPT_abs_250,
                           method = "fREML",
                           discrete = TRUE,
                           family = gaussian())


summary(maxPPT_distEdge_250_G)
gam.check(maxPPT_distEdge_250_G)

#AIC
AIC_maxPPT_distEdge_250_G <- AIC(maxPPT_distEdge_250_G)

setwd(wd_models)
saveRDS(maxPPT_distEdge_250_G, 'maxPPT_distEdge_250_G')


plot.gam(maxPPT_distEdge_250_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS

#model GS (presence) 
maxPPT_distEdge_GS <- bam(avg_Max_PPT_SHAP ~
                          s(distEdge, k = 4, m = 2)
                        + s(distEdge, species, k = 4, bs = "fs", m = 2),
                        data = res_maxPPT_abs,
                        method = "fREML", discrete = TRUE)

summary(maxPPT_distEdge_GS)
gam.check(maxPPT_distEdge_GS)

#AIC
AIC_maxPPT_distEdge_GS <- AIC(maxPPT_distEdge_GS)

setwd(wd_models)
saveRDS(maxPPT_distEdge_GS, 'maxPPT_distEdge_GS')


plot.gam(maxPPT_distEdge_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.5, 0.5),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
maxPPT_distEdge_250_GS <- bam(avg_Max_PPT_SHAP ~
                              s(distEdge, k = 4, m = 2)
                            + s(distEdge, species, k = 4, bs = "fs", m = 2),
                            data = res_maxPPT_abs_250,
                            method = "fREML",
                            discrete = TRUE,
                            family = gaussian())


summary(maxPPT_distEdge_250_GS)
gam.check(maxPPT_distEdge_250_GS)

#AIC
AIC_maxPPT_distEdge_250_GS <- AIC(maxPPT_distEdge_250_GS)

setwd(wd_models)
saveRDS(maxPPT_distEdge_250_GS, 'maxPPT_distEdge_250_GS')


plot.gam(maxPPT_distEdge_250_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(maxPPT_distEdge_250_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxPPT_distEdge_250_GS.csv', row.names = F)




######### Select only sps not too limited by the ocean  (max 20%) #########

res_minPPT_abs_20 <- res_minPPT_abs[
  as.numeric(gsub("%", "", res_minPPT_abs$perc_ocean)) <= 20, ]

res_minPPT_abs_250_20 <- res_minPPT_abs_250[
  as.numeric(gsub("%", "", res_minPPT_abs_250$perc_ocean)) <= 20, ]

res_meanPPT_abs_20 <- res_meanPPT_abs[
  as.numeric(gsub("%", "", res_meanPPT_abs$perc_ocean)) <= 20, ]

res_meanPPT_abs_250_20 <- res_meanPPT_abs_250[
  as.numeric(gsub("%", "", res_meanPPT_abs_250$perc_ocean)) <= 20, ]

res_maxPPT_abs_20 <- res_maxPPT_abs[
  as.numeric(gsub("%", "", res_maxPPT_abs$perc_ocean)) <= 20, ]

res_maxPPT_abs_250_20 <- res_maxPPT_abs_250[
  as.numeric(gsub("%", "", res_maxPPT_abs_250$perc_ocean)) <= 20, ]

# SHAP values X distance from edge (GAM)

### minT

#model G (presence)
minPPT_distEdge_20_G <- bam(avg_Min_PPT_SHAP ~
                            s(distEdge, k = 4, bs = "tp")
                          + s(species, k = n_sps_minPPT, bs = "re"),
                          data = res_minPPT_abs_20,
                          method = "fREML",
                          discrete = TRUE,
                          family = gaussian())


summary(minPPT_distEdge_20_G)
gam.check(minPPT_distEdge_20_G)


#AIC
AIC_minT_distEdge_20_G <- AIC(minT_distEdge_20_G)

setwd(wd_models)
saveRDS(minPPT_distEdge_20_G, 'minPPT_distEdge_20_G')


plot.gam(minPPT_distEdge_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
minPPT_distEdge_250_20_G <- bam(avg_Min_PPT_SHAP ~
                                s(distEdge, k = 4, bs = "tp")
                              + s(species, k = n_sps_minPPT_mod, bs = "re"),
                              data = res_minPPT_abs_250_20,
                              method = "fREML",
                              discrete = TRUE,
                              family = gaussian())


summary(minPPT_distEdge_250_20_G)
gam.check(minPPT_distEdge_250_20_G)


#AIC
AIC_minPPT_distEdge_250_20_G <- AIC(minPPT_distEdge_250_20_G)

setwd(wd_models)
saveRDS(minPPT_distEdge_250_20_G, 'minPPT_distEdge_250_20_G')


plot.gam(minPPT_distEdge_250_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minPPT_distEdge_20_GS <- bam(avg_Min_PPT_SHAP ~
                             s(distEdge, k = 4, m = 2)
                           + s(distEdge, species, k = 4, bs = "fs", m = 2),
                           data = res_minPPT_abs_20,
                           method = "fREML", discrete = TRUE)

summary(minPPT_distEdge_20_GS)
gam.check(minPPT_distEdge_20_GS)

#AIC
AIC_minPPT_distEdge_20_GS <- AIC(minPPT_distEdge_20_GS)

setwd(wd_models)
saveRDS(minPPT_distEdge_20_GS, 'minPPT_distEdge_20_GS')

plot.gam(minPPT_distEdge_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.4, 0.4),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (250 presence)
minPPT_distEdge_250_20_GS <- bam(avg_Min_PPT_SHAP ~
                                 s(distEdge, k = 4, m = 2)
                               + s(distEdge, species, k = 4, bs = "fs", m = 2),
                               data = res_minPPT_abs_250_20,
                               method = "fREML",
                               discrete = TRUE,
                               family = gaussian())


summary(minPPT_distEdge_250_20_GS)
gam.check(minPPT_distEdge_250_20_GS)

#AIC
AIC_minPPT_distEdge_250_20_GS <- AIC(minPPT_distEdge_250_20_GS)

setwd(wd_models)
saveRDS(minPPT_distEdge_250_20_GS, 'minPPT_distEdge_250_20_GS')


plot.gam(minPPT_distEdge_250_20_GS, select = 1, residuals = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.15),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(minPPT_distEdge_250_20_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minPPT_distEdge_250_20_GS.csv', row.names = F)



### meanT

#model G (presence)
meanPPT_distEdge_20_G <- bam(avg_Mean_PPT_SHAP ~
                             s(distEdge, k = 4, bs = "tp")
                           + s(species, k = n_sps_meanPPT, bs = "re"),
                           data = res_meanPPT_abs_20,
                           method = "fREML",
                           discrete = TRUE,
                           family = gaussian())


summary(meanPPT_distEdge_20_G)

#AIC
AIC_meanPPT_distEdge_20_G <- AIC(meanPPT_distEdge_20_G)

setwd(wd_models)
saveRDS(meanPPT_distEdge_20_G, 'meanPPT_distEdge_20_G')


plot.gam(meanPPT_distEdge_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.20, 0.30),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (250 presence)
meanPPT_distEdge_250_20_G <- bam(avg_Mean_PPT_SHAP ~
                                 s(distEdge, k = 4, bs = "tp")
                               + s(species, k = n_sps_meanPPT_mod, bs = "re"),
                               data = res_meanPPT_abs_250_20,
                               method = "fREML",
                               discrete = TRUE,
                               family = gaussian())


summary(meanPPT_distEdge_250_20_G)
gam.check(meanPPT_distEdge_20_G)

#AIC
AIC_meanPPT_distEdge_250_20_G <- AIC(meanPPT_distEdge_250_20_G)

setwd(wd_models)
saveRDS(meanPPT_distEdge_250_20_G, 'meanPPT_distEdge_250_20_G')


plot.gam(meanPPT_distEdge_250_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.08, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence)
meanPPT_distEdge_20_GS <- bam(avg_Mean_PPT_SHAP ~
                              s(distEdge, k = 4, m = 2)
                            + s(distEdge, species, k = 4, bs = "fs", m = 2),
                            data = res_meanPPT_abs_20,
                            method = "fREML", discrete = TRUE)


summary(meanPPT_distEdge_20_GS)
gam.check(meanPPT_distEdge_20_GS)

#AIC
AIC_meanPPT_distEdge_20_GS <- AIC(meanPPT_distEdge_20_GS)

setwd(wd_models)
saveRDS(meanPPT_distEdge_20_GS, 'meanPPT_distEdge_20_GS')


plot.gam(meanPPT_distEdge_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
meanPPT_distEdge_250_20_GS <- bam(avg_Mean_PPT_SHAP ~
                                  s(distEdge, k = 4, m = 2)
                                + s(distEdge, species, k = 4, bs = "fs", m = 2),
                                data = res_meanPPT_abs_250_20,
                                method = "fREML",
                                discrete = TRUE,
                                family = gaussian())


summary(meanPPT_distEdge_250_20_GS)
gam.check(meanPPT_distEdge_250_20_GS)

#AIC
AIC_meanPPT_distEdge_250_20_GS <- AIC(meanPPT_distEdge_250_20_GS)

setwd(wd_models)
saveRDS(meanPPT_distEdge_250_20_GS, 'meanPPT_distEdge_250_20_GS')


plot.gam(meanPPT_distEdge_250_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.2, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(meanPPT_distEdge_250_20_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanPPT_distEdge_250_20_GS.csv', row.names = F)


### maxT

#model G (presence) - fast bam
maxPPT_distEdge_20_G <- bam(avg_Max_PPT_SHAP ~
                            s(distEdge, k = 4, bs = "tp")
                          + s(species, k = n_sps_maxPPT, bs = "re"),
                          data = res_maxPPT_abs_20,
                          method = "fREML",
                          discrete = TRUE,
                          family = gaussian())


summary(maxPPT_distEdge_20_G)
gam.check(maxPPT_distEdge_20_G)

#AIC
AIC_maxPPT_distEdge_20_G <- AIC(maxPPT_distEdge_20_G)

setwd(wd_models)
saveRDS(maxPPT_distEdge_20_G, 'maxPPT_distEdge_20_G')


plot.gam(maxPPT_distEdge_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.5, 0.5),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model G (presence 250 km) 
maxPPT_distEdge_250_20_G <- bam(avg_Max_PPT_SHAP ~
                                s(distEdge, k = 4, bs = "tp")
                              + s(species, k = n_sps_maxPPT_mod, bs = "re"),
                              data = res_maxPPT_abs_250_20,
                              method = "fREML",
                              discrete = TRUE,
                              family = gaussian())


summary(maxPPT_distEdge_250_20_G)
gam.check(maxPPT_distEdge_250_20_G)

#AIC
AIC_maxPPT_distEdge_250_20_G <- AIC(maxPPT_distEdge_250_20_G)

setwd(wd_models)
saveRDS(maxPPT_distEdge_250_20_G, 'maxPPT_distEdge_250_20_G')


plot.gam(maxPPT_distEdge_250_20_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS

#model GS (presence) 
maxPPT_distEdge_20_GS <- bam(avg_Max_PPT_SHAP ~
                             s(distEdge, k = 4, m = 2)
                           + s(distEdge, species, k = 4, bs = "fs", m = 2),
                           data = res_maxPPT_abs_20,
                           method = "fREML", discrete = TRUE)

summary(maxPPT_distEdge_20_GS)
gam.check(maxPPT_distEdge_20_GS)

#AIC
AIC_maxPPT_distEdge_20_GS <- AIC(maxPPT_distEdge_20_GS)

setwd(wd_models)
saveRDS(maxPPT_distEdge_20_GS, 'maxPPT_distEdge_20_GS')


plot.gam(maxPPT_distEdge_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.5, 0.5),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence 250 km)
maxPPT_distEdge_250_20_GS <- bam(avg_Max_PPT_SHAP ~
                                 s(distEdge, k = 4, m = 2)
                               + s(distEdge, species, k = 4, bs = "fs", m = 2),
                               data = res_maxPPT_abs_250_20,
                               method = "fREML",
                               discrete = TRUE,
                               family = gaussian())


summary(maxPPT_distEdge_250_20_GS)
gam.check(maxPPT_distEdge_250_20_GS)

#AIC
AIC_maxPPT_distEdge_250_20_GS <- AIC(maxPPT_distEdge_250_20_GS)

setwd(wd_models)
saveRDS(maxPPT_distEdge_250_20_GS, 'maxPPT_distEdge_250_20_GS')


plot.gam(maxPPT_distEdge_250_20_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.15),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(maxPPT_distEdge_250_20_GS,
                     select = "s(distEdge)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxPPT_distEdge_250_20_GS.csv', row.names = F)




# SHAP values X relative polewardness (GAM)

### minPPT

#model G (presence) 
minPPT_relPol_G <- bam(avg_Min_PPT_SHAP ~
                       s(relPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_minPPT, bs = "re"),
                     data = res_minPPT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(minPPT_relPol_G)
gam.check(minPPT_relPol_G)

#AIC
AIC_minPPT_relPol_G <- AIC(minPPT_relPol_G)

setwd(wd_models)
saveRDS(minPPT_relPol_G, 'minPPT_relPol_G')


plot.gam(minPPT_relPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minPPT_relPol_GS <- bam(avg_Min_PPT_SHAP ~
                        s(relPolewardness, k = 4, m = 2)
                      + s(relPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_minPPT,
                      method = "fREML",
                      discrete = TRUE)


summary(minPPT_relPol_GS)
gam.check(minPPT_relPol_GS)

#AIC
AIC_minPPT_relPol_GS <- AIC(minPPT_relPol_GS)

setwd(wd_models)
saveRDS(minPPT_relPol_GS, 'minPPT_relPol_GS')


plot.gam(minPPT_relPol_GS, select = 1, residuals = F, shade = T,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.2),
         cex.lab = 2, cex.axis = 1.5) #save 800


#test stats for significant change points from First Derivative of GAM smooth

#calculate derivatives
deriv <- derivatives(minPPT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minPPT_relPol_GS.csv', row.names = F)

### meanPPT

#model G (presence) 
meanPPT_relPol_G <- bam(avg_Mean_PPT_SHAP ~
                        s(relPolewardness, k = 4, bs = "tp")
                      + s(species, k = n_sps_meanPPT, bs = "re"),
                      data = res_meanPPT,
                      method = "fREML",
                      discrete = TRUE,
                      family = gaussian())


summary(meanPPT_relPol_G)
gam.check(meanPPT_relPol_G)

#AIC
AIC_meanPPT_relPol_G <- AIC(meanPPT_relPol_G)

setwd(wd_models)
saveRDS(meanPPT_relPol_G, 'meanPPT_relPol_G')


plot.gam(meanPPT_relPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.15, 0.06),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
meanPPT_relPol_GS <- bam(avg_Mean_PPT_SHAP ~
                         s(relPolewardness, k = 4, m = 2)
                       + s(relPolewardness, species, k = 4, bs = "fs", m = 2),
                       data = res_meanPPT,
                       method = "fREML",
                       discrete = TRUE)


summary(meanPPT_relPol_GS)
gam.check(meanPPT_relPol_GS)

#AIC
AIC_meanPPT_relPol_GS <- AIC(meanPPT_relPol_GS)

setwd(wd_models)
saveRDS(meanPPT_relPol_GS, 'meanPPT_relPol_GS')


plot.gam(meanPPT_relPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.25, 0.25),
         cex.lab = 2, cex.axis = 1.5) #save 800




library(gratia)

#calculate derivatives
deriv <- derivatives(meanPPT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate

#identify points where derivative is significantly different from zero
sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanPPT_relPol_GS.csv', row.names = F)


### maxT

#model G (presence) 
maxPPT_relPol_G <- bam(avg_Max_PPT_SHAP ~
                       s(relPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_maxPPT, bs = "re"),
                     data = res_maxPPT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(maxPPT_relPol_G)
gam.check(maxPPT_relPol_G)

#AIC
AIC_maxPPT_relPol_G <- AIC(maxPPT_relPol_G)

setwd(wd_models)
saveRDS(maxPPT_relPol_G, 'maxPPT_relPol_G')


plot.gam(maxPPT_relPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.2, 0.15),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence)
maxPPT_relPol_GS <- bam(avg_Max_PPT_SHAP ~
                        s(relPolewardness, k = 4, m = 2)
                      + s(relPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_maxPPT,
                      method = "fREML",
                      discrete = TRUE)


summary(maxPPT_relPol_GS)
gam.check(maxPPT_relPol_GS)

#AIC
AIC_maxPPT_relPol_GS <- AIC(maxPPT_relPol_GS)

setwd(wd_models)
saveRDS(maxPPT_relPol_GS, 'maxPPT_relPol_GS')


plot.gam(maxPPT_relPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.1, 0.1),
         cex.lab = 2, cex.axis = 1.5) #save 800

#calculate derivatives
deriv <- derivatives(maxPPT_relPol_GS,
                     select = "s(relPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxPPT_relPol_GS.csv', row.names = F)


# SHAP values X absolute polewardness (GAM)

#set par for plotting
par(mar = c(6,6,6,6))

### minT

#model G (presence)
minPPT_absPol_G <- bam(avg_Min_PPT_SHAP ~
                       s(absPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_minPPT, bs = "re"),
                     data = res_minPPT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(minPPT_absPol_G)
gam.check(minPPT_absPol_G)

#AIC
AIC_minPPT_absPol_G <- AIC(minPPT_absPol_G)

setwd(wd_models)
saveRDS(minPPT_absPol_G, 'minPPT_absPol_G')


plot.gam(minPPT_absPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.5, 0.5),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
minPPT_absPol_GS <- bam(avg_Min_PPT_SHAP ~
                        s(absPolewardness, k = 4, m = 2)
                      + s(absPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_minPPT,
                      method = "fREML",
                      discrete = TRUE)


summary(minPPT_absPol_GS)
gam.check(minPPT_absPol_GS)

#AIC
AIC_minPPT_absPol_GS <- AIC(minPPT_absPol_GS)

setwd(wd_models)
saveRDS(minPPT_absPol_GS, 'minPPT_absPol_GS')


plot.gam(minPPT_absPol_GS, select = 1, residuals = F,
         ylab = 'SHAP value',
         ylim = c(-0.8, 0.8),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(minPPT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'minPPT_absPol_GS.csv', row.names = F)




### meanT

#model G (presence)
meanPPT_absPol_G <- bam(avg_Mean_PPT_SHAP ~
                        s(absPolewardness, k = 4, bs = "tp")
                      + s(species, k = n_sps_meanPPT, bs = "re"),
                      data = res_meanPPT,
                      method = "fREML",
                      discrete = TRUE,
                      family = gaussian())


summary(meanPPT_absPol_G)
gam.check(meanPPT_absPol_G)

#AIC
AIC_meanPPT_absPol_G <- AIC(meanPPT_absPol_G)

setwd(wd_models)
saveRDS(meanPPT_absPol_G, 'meanPPT_absPol_G')


plot.gam(meanPPT_absPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.8, 0.8),
         cex.lab = 2, cex.axis = 1.5) #save 800


#model GS (presence)
meanPPT_absPol_GS <- bam(avg_Mean_PPT_SHAP ~
                         s(absPolewardness, k = 4, m = 2)
                       + s(absPolewardness, species, k = 4, bs = "fs", m = 2),
                       data = res_meanPPT,
                       method = "fREML",
                       discrete = TRUE)


summary(meanPPT_absPol_GS)
gam.check(meanPPT_absPol_GS)

#AIC
AIC_meanPPT_absPol_GS <- AIC(meanPPT_absPol_GS)

setwd(wd_models)
saveRDS(meanPPT_absPol_GS, 'meanPPT_absPol_GS')


plot.gam(meanPPT_absPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-4, 1.5),
         cex.lab = 2, cex.axis = 1.5) #save 800

#calculate derivatives
deriv <- derivatives(meanPPT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'meanPPT_absPol_GS.csv', row.names = F)



### maxT

#model G (presence)
maxPPT_absPol_G <- bam(avg_Max_PPT_SHAP ~
                       s(absPolewardness, k = 4, bs = "tp")
                     + s(species, k = n_sps_maxPPT, bs = "re"),
                     data = res_maxPPT,
                     method = "fREML",
                     discrete = TRUE,
                     family = gaussian())


summary(maxPPT_absPol_G)
gam.check(maxPPT_absPol_G)

#AIC
AIC_maxPPT_absPol_G <- AIC(maxPPT_absPol_G)

setwd(wd_models)
saveRDS(maxPPT_absPol_G, 'maxPPT_absPol_G')


plot.gam(maxPPT_absPol_G, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-0.3, 0.3),
         cex.lab = 2, cex.axis = 1.5) #save 800



#model GS (presence)
maxPPT_absPol_GS <- bam(avg_Max_PPT_SHAP ~
                        s(absPolewardness, k = 4, m = 2)
                      + s(absPolewardness, species, k = 4, bs = "fs", m = 2),
                      data = res_maxPPT,
                      method = "fREML",
                      discrete = TRUE)


summary(maxPPT_absPol_GS)
gam.check(maxPPT_absPol_GS)

#AIC
AIC_maxPPT_absPol_GS <- AIC(maxPPT_absPol_GS)

setwd(wd_models)
saveRDS(maxPPT_absPol_GS, 'maxPPT_absPol_GS')


plot.gam(maxPPT_absPol_GS, select = 1, residuals = F, shade = F,
         ylab = 'SHAP value',
         ylim = c(-2.8, 1.3),
         cex.lab = 2, cex.axis = 1.5) #save 800


#calculate derivatives
deriv <- derivatives(maxPPT_absPol_GS,
                     select = "s(absPolewardness)",
                     type = "central",
                     n = 500)  # number of points to evaluate


sig_points <- with(deriv, !(.lower_ci <= 0 & .upper_ci >= 0))

#add significance to the derivatives data frame
deriv$sig <- sig_points
df <- deriv

#save significance table
setwd(wd_sig)
write.csv(df, 'maxPPT_absPol_GS.csv', row.names = F)
