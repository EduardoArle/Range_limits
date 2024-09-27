#load libraries
library(mgcv)

#list wds
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/GAMs/Models'

#load models
setwd(wd_models)
models <- lapply(list.files(pattern = '_G'), readRDS)

#rename models
names(models) <- list.files(pattern = '_G')

k.check(models[[2]], n.rep = 4000)

rsd <- residuals(models[[2]])
qq.gam(models[[2]], rep=100)
plot(fitted(models[[2]]), rsd)
plot(dat$x0,rsd); plot(dat$x1,rsd)




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