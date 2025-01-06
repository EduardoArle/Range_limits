#list WDs
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Tables'

#load mistressfile
setwd(wd_tables)
mistress_file <-  read.csv('Mistressfile.csv')

#count number of species with ok correl for each model
minT_meanPPT <- t(table(abs(mistress_file$Cor_vars_minT) <= 0.7))
meanT_meanPPT <- t(table(abs(mistress_file$Cor_vars_meanT) <= 0.7))
maxT_meanPPT <- t(table(abs(mistress_file$Cor_vars_maxT) <= 0.7))

minPPT_meanT <- t(table(abs(mistress_file$Cor_vars_minPPT) <= 0.7))
meanPPT_meanT <- t(table(abs(mistress_file$Cor_vars_meanPPT) <= 0.7))
maxPPT_meanT <- t(table(abs(mistress_file$Cor_vars_maxPPT) <= 0.7))



