#load packages
library(brms)

#set WDs
wd_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All_species_analysis'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Bayesian_models'

#load data
setwd(wd_data)
table <- read.csv('Results_all_sps.csv')

#add absolute latitude to the table
table$abs_lat <- abs(table$decimalLatitude)

#add square root of number of occurrences to the table
table$sqrt_n_occ <- sqrt(table$n_occ)

#select only columns of interest
cols_int <- c('species', 'key', 'abs_lat', 'rangeSize', 'distEdge',
              'relPolarwardness', 'elevation', 'biome', 'bodyMass',
              'sqrt_n_occ', 'n_occ', 'order', 'Min_T_SHAP', 'Mean_T_SHAP',
              'Max_T_SHAP', 'Min_PPT_SHAP', 'Mean_PPT_SHAP', 'Max_PPT_SHAP',
              'NOTE')

table2 <- table[, which(names(table) %in% cols_int)]

#select only species whose range does not cross the equator
table3 <- table2[which(is.na(table2$NOTE)),]

#select only species which at least 20 occurrences
table4 <- table3[which(table3$n_occ >= 20),]

#transform biome into factor
table4$biome <- as.factor(table4$biome)

#Polewardness

##max_temp - prior 1
max_temp_model <- brm(data = table4,
                      formula = Max_T_SHAP ~ relPolarwardness + 
                                (1 + relPolarwardness | biome) + 
                                (1 + relPolarwardness | species),
                      chains = 3,
                      iter = 4000,
                      warmup = 2000,
                      cores = 3,
                      prior = set_prior('normal(0,5)',
                      class="b",
                      coef="relPolarwardness"),
                      family=gaussian()) #22 min

#save model
setwd(wd_models)
saveRDS(max_temp_model, file = "Max_temp_model.RDS") 

bayes_R2(max_temp_model)
summary(max_temp_model)
max_temp_model_gg=ggeffects::ggpredict(max_temp_model,terms="relPolarwardness[all]")

#Diagnosis
brms::pp_check(max_temp_model,ndraws=99)# posterior predictive checks
plot(max_temp_model)
