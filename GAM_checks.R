#load libraries
library(mgcv)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models'
wd_models_k <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models_k_tests'
wd_models_restr <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models_restricted'

wd_res_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Models/Console_results'
wd_plots_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Models/Plots'

wd_res_models_k <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Models_k_test/Console_results'
wd_plots_models_k <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Models_k_test/Plots'

wd_res_models_restr <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Models_restricted/Console_results'
wd_plots_models_restr <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Models_restricted/Plots'

#load models
setwd(wd_models)
models <- lapply(list.files(pattern = '_G'), readRDS)

#rename models
names(models) <- list.files(pattern = '_G')

#run checks
for(i in 1:length(models))
{
  #create text file to save the output printed in the console
  sink(file = paste0(wd_res_models, '/', names(models)[i], '.txt'))
  
  #set the WD for the plots
  setwd(wd_plots_models)
  #create pdf for plots
  pdf(file = paste0(names(models)[i], ".pdf"))
  #set the mfrow for a two by two grid
  par(mfrow = c(2,2))
  #run the gam.check that will produce plots and results in the console
  gam.check(models[[i]], k.rep = 1000)
  
  #save the result printed in the console
  sink(file = NULL)
  
  #close the connection to the pdf with the plots
  dev.off()
  
  print(i)
}


#load models_k
setwd(wd_models_k)
models_k <- lapply(list.files(pattern = '_G'), readRDS)

#rename models
names(models_k) <- list.files(pattern = '_G')

#run checks
for(i in 1:length(models_k))
{
  #create text file to save the output printed in the console
  sink(file = paste0(wd_res_models_k, '/', names(models_k)[i], '.txt'))
  
  #set the WD for the plots
  setwd(wd_plots_models_k)
  #create pdf for plots
  pdf(file = paste0(names(models_k)[i], ".pdf"))
  #set the mfrow for a two by two grid
  par(mfrow = c(2,2))
  #run the gam.check that will produce plots and results in the console
  gam.check(models_k[[i]], k.rep = 1000)
  
  #save the result printed in the console
  sink(file = NULL)
  
  #close the connection to the pdf with the plots
  dev.off()
  
  print(i)
}


#load models_restricted
setwd(wd_models_restr)
models_restr <- lapply(list.files(pattern = '_G'), readRDS)

#rename models
names(models_restr) <- list.files(pattern = '_G')

#run checks
for(i in 1:length(models_restr))
{
  #create text file to save the output printed in the console
  sink(file = paste0(wd_res_models_restr, '/', names(models_restr)[i], '.txt'))
  
  #set the WD for the plots
  setwd(wd_plots_models_restr)
  #create pdf for plots
  pdf(file = paste0(names(models_restr)[i], ".pdf"))
  #set the mfrow for a two by two grid
  par(mfrow = c(2,2))
  #run the gam.check that will produce plots and results in the console
  gam.check(models_restr[[i]], k.rep = 1000)
  
  #save the result printed in the console
  sink(file = NULL)
  
  #close the connection to the pdf with the plots
  dev.off()
  
  print(i)
}



### TESTS ###

#select one model for tests
setwd(wd_models)
meanT_relPolewarness_GS <- readRDS('meanT_relPolewarness_GS')

wd_test_k.rep_10000 <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Tests/K.rep_10000'
  
setwd(wd_test_k.rep_10000)

#run checks
for(i in 1:10)
{
  #create text file to save the output printed in the console
  sink(file = paste0(wd_test_k.rep_10000, '/Console_results/k.rep_10000_', i, '.txt'))
  
  #set the WD for the plots
  setwd(paste0(wd_test_k.rep_10000, '/Plots'))
  #create pdf for plots
  pdf(file = paste0('k.rep_10000_', i, '.pdf'))
  #set the mfrow for a two by two grid
  par(mfrow = c(2,2))
  #run the gam.check that will produce plots and results in the console
  gam.check(meanT_relPolewarness_GS, k.rep = 10000)
  
  #save the result printed in the console
  sink(file = NULL)
  
  #close the connection to the pdf with the plots
  dev.off()
  
  print(i)
}


wd_test_k.sample_10000 <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/gam.check_results/Tests/K.rep_10000'

setwd(wd_test_k.rep_10000)

#run checks
for(i in 1:10)
{
  #create text file to save the output printed in the console
  sink(file = paste0(wd_test_k.rep_10000, '/Console_results/k.rep_10000_', i, '.txt'))
  
  #set the WD for the plots
  setwd(paste0(wd_test_k.rep_10000, '/Plots'))
  #create pdf for plots
  pdf(file = paste0('k.rep_10000_', i, '.pdf'))
  #set the mfrow for a two by two grid
  par(mfrow = c(2,2))
  #run the gam.check that will produce plots and results in the console
  gam.check(meanT_relPolewarness_GS, k.rep = 10000)
  
  #save the result printed in the console
  sink(file = NULL)
  
  #close the connection to the pdf with the plots
  dev.off()
  
  print(i)
}


