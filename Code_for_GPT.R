#load libraries
library(mgcv); library(itsadug)

#read results table with shap values results and covariates for all occurrences across species
setwd(wd_tables)
results <- read.csv('20241210_Results_all_sps.csv')

#change names that are wrong
names(results)[c(10,12)] <- c("absPolewardness","relPolewardness" )

#get names of all columns in tables
col_names <- names(results)

#model G
minPPT_relPolewarness_G <- gam(avg_Min_PPT_SHAP ~
                                 s(relPolewardness, k=4, bs="tp") 
                               + s(species, k=503, bs="re"),
                               data = results,
                               method="REML",
                               family="gaussian")

#model GS
minPPT_relPolewarness_GS <- gam(avg_Min_PPT_SHAP ~
                                  s(relPolewardness, k=4, m=2) 
                                + s(relPolewardness, species, k=4, bs="fs", m=2),
                                data = results,
                                method="REML")
