#############################
### extract ENM performance ###
#############################

#load libraries
library(data.table)

#list WDs
wd_shap <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SHAP_results'
wd_out <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SI/Tables'

#list species folders
setwd(wd_shap)
sps_dirs <- list.dirs(wd_shap, recursive = FALSE, full.names = TRUE)

#check
length(sps_dirs)
head(basename(sps_dirs))

#list replicate result files inside species folders
rep_files <- list.files(sps_dirs, pattern = '\\.csv$', recursive = TRUE,
                        full.names = TRUE)

#check
length(rep_files)
head(basename(rep_files))


#############################
### empty results table ###
#############################



#empty table for ENM performance
perf_all <- data.frame()



#############################
### extract ENM performance ###
#############################

for(i in 1:length(rep_files))
{
  #progress
  print(i)
  
  #get file name
  file_i <- basename(rep_files[i])
  
  #remove .csv
  file_i2 <- sub('.csv$', '', file_i)
  
  #get species name
  sps_i <- basename(dirname(rep_files[i]))
  
  #get variable
  var_i <- sub(paste0(sps_i, '_'), '', file_i2)
  var_i <- sub('_run_.*$', '', var_i)
  
  #get run
  run_i <- sub('.*_run_', '', file_i2)
  run_i <- sub('_rep_.*$', '', run_i)
  
  #get replicate
  rep_i <- sub('.*_rep_', '', file_i2)
  
  #read file
  dat_i <- fread(rep_files[i])
  
  #extract metrics
  auc_i <- dat_i$AUC_test[1]
  tss_i <- dat_i$TSS_test[1]
  
  #store results
  tmp <- data.frame(
    species = sps_i,
    variable = var_i,
    run = as.numeric(run_i),
    replicate = as.numeric(rep_i),
    AUC = auc_i,
    TSS = tss_i
  )
  
  perf_all <- rbind(perf_all, tmp)
}

head(perf_all)
summary(perf_all$AUC)
summary(perf_all$TSS)
nrow(perf_all)



#############################
### save ENM performance ###
#############################


#set wd
setwd(wd_out)

#save R object
saveRDS(perf_all, 'ENM_performance_all_models.rds')

#save csv
write.csv(perf_all, 'ENM_performance_all_models.csv', row.names = FALSE)



#############################
### retained models ###
#############################



#subset retained models
perf_retained <- perf_all[perf_all$AUC > 0.7 &
                            perf_all$TSS > 0.4,]

#check retained models
nrow(perf_retained)

#proportion retained
nrow(perf_retained) / nrow(perf_all)

#summarise retained AUC
summary(perf_retained$AUC)

#summarise retained TSS
summary(perf_retained$TSS)



#############################
### histogram AUC and TSS ###
#############################



#set layout
par(mfrow = c(2,2),
    mar = c(4.5,4.5,2,2),
    mgp = c(2.5,0.8,0),
    cex.axis = 1.2,
    cex.lab = 1.3)

#AUC, all models
hist(perf_all$AUC, breaks = seq(0,1,0.025), xlim = c(0,1),
     main = '', xlab = 'AUC', ylab = 'Frequency',
     col = 'grey80', border = 'grey30')
mtext('All models', side = 3, line = 0.3, cex = 1.2)

#AUC, retained models
hist(perf_retained$AUC, breaks = seq(0.7,1,0.01), xlim = c(0.7,1),
     main = '', xlab = 'AUC', ylab = 'Frequency',
     col = 'grey80', border = 'grey30')
mtext('Retained models', side = 3, line = 0.3, cex = 1.2)

#TSS, all models
hist(perf_all$TSS, breaks = seq(-1,1,0.05), xlim = c(-1,1),
     main = '', xlab = 'TSS', ylab = 'Frequency',
     col = 'grey80', border = 'grey30')
mtext('All models', side = 3, line = 0.3, cex = 1.2)

#TSS, retained models
hist(perf_retained$TSS, breaks = seq(0.4,1,0.02), xlim = c(0.4,1),
     main = '', xlab = 'TSS', ylab = 'Frequency',
     col = 'grey80', border = 'grey30')
mtext('Retained models', side = 3, line = 0.3, cex = 1.2)



#############################
### retained replicates ###
#############################



#count total replicates per species
n_total <- aggregate(AUC ~ species, data = perf_all, FUN = length)
names(n_total)[2] <- 'n_total'

#count retained replicates per species
n_retained <- aggregate(AUC ~ species, data = perf_retained, FUN = length)
names(n_retained)[2] <- 'n_retained'

#merge tables
retention <- merge(n_total, n_retained, by = 'species', all.x = TRUE)

#replace NA with zero
retention$n_retained[is.na(retention$n_retained)] <- 0

#calculate retained proportion
retention$prop_retained <- retention$n_retained / retention$n_total



#############################
### summarise retention ###
#############################



#number retained
summary(retention$n_retained)

#proportion retained
summary(retention$prop_retained)

#species with no retained models
sum(retention$n_retained == 0)




