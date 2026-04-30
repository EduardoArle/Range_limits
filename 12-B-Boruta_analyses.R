#load packages
library(Boruta)

#list wds
wd_slopes <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Slopes'

wd_minT_POL_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/minT_POL_relPol'
wd_minT_EQ_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/minT_EQ_relPol'
wd_meanT_POL_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/meanT_POL_relPol'
wd_meanT_EQ_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/meanT_EQ_relPol'
wd_maxT_POL_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/maxT_POL_relPol'
wd_maxT_EQ_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/maxT_EQ_relPol'

wd_minT_distEdge <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/minT_distEdge'
wd_meanT_distEdge <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/meanT_distEdge'
wd_maxT_distEdge <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/maxT_distEdge'

wd_minPPT_POL_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/minPPT_POL_relPol'
wd_minPPT_EQ_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/minPPT_EQ_relPol'
wd_meanPPT_POL_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/meanPPT_POL_relPol'
wd_meanPPT_EQ_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/meanPPT_EQ_relPol'
wd_maxPPT_POL_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/maxPPT_POL_relPol'
wd_maxPPT_EQ_relPol <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/maxPPT_EQ_relPol'

wd_minPPT_distEdge <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/minPPT_distEdge'
wd_meanPPT_distEdge <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/meanPPT_distEdge'
wd_maxPPT_distEdge <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta/maxPPT_distEdge'


#read slopes table
setwd(wd_slopes)
slopes <- read.csv('20260428_Slopes_split.csv')


###### MinT vs RELATIVE POLEWARDNESS ######


#select only species that had lower correl between vars
s_minT_relPol <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_relPol <- s_minT_relPol[complete.cases(s_minT_relPol$Cor_vars_minT),]

#separate responses
s_minT_relPol_EQ <- s_minT_relPol[
  complete.cases(s_minT_relPol$slope_EQ_minT_relPol), ]
s_minT_relPol_POL <- s_minT_relPol[
  complete.cases(s_minT_relPol$slope_POL_minT_relPol), ]

#transform categorical variable in factor
s_minT_relPol_EQ$order <- as.factor(s_minT_relPol_EQ$order)
s_minT_relPol_POL$order <- as.factor(s_minT_relPol_POL$order)

#create new columns with log values (boruta does not accept it...)
s_minT_relPol_EQ$rangeSize_log10 <- log(s_minT_relPol_EQ$rangeSize, 10)
s_minT_relPol_POL$rangeSize_log10 <- log(s_minT_relPol_POL$rangeSize, 10)

s_minT_relPol_EQ$nOcc_EQ_log <- log(s_minT_relPol_EQ$nOcc_EQ)
s_minT_relPol_POL$nOcc_POL_log <- log(s_minT_relPol_POL$nOcc_Pol)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.minT_EQ_relPol <- Boruta(slope_EQ_minT_relPol ~
                                    rangeSize_log10 +
                                    rangeLoc +
                                    latAmplitude +
                                    roundness +
                                    nOcc_EQ_log +
                                    order +
                                    elevMedian +
                                    bodyMass +
                                    elevAmplitude,
                                  data = s_minT_relPol_EQ, doTrace = 2, ntree = 500,
                                  maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.minT_EQ_relPol)
  setwd(wd_minT_EQ_relPol)
  write.csv(stats, paste0('minT_EQ_relPol_', i, '.csv'))
}


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.minT_POL_relPol <- Boruta(slope_POL_minT_relPol ~
                                     rangeSize_log10 +
                                     rangeLoc +
                                     latAmplitude +
                                     roundness +
                                     nOcc_POL_log +
                                     order +
                                     elevMedian +
                                     bodyMass +
                                     elevAmplitude,
                                data = s_minT_relPol_POL, doTrace = 2, ntree = 500,
                                maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.minT_POL_relPol)
  setwd(wd_minT_POL_relPol)
  write.csv(stats, paste0('minT_POL_relPol_', i, '.csv'))
}




###### meanT vs RELATIVE POLEWARDNESS ######


#select only species that had lower correl between vars
s_meanT_relPol <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,]
s_meanT_relPol <- s_meanT_relPol[complete.cases(s_meanT_relPol$Cor_vars_meanT),]

#separate responses
s_meanT_relPol_EQ <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_EQ_meanT_relPol), ]
s_meanT_relPol_POL <- s_meanT_relPol[
  complete.cases(s_meanT_relPol$slope_POL_meanT_relPol), ]

#transform categorical variable in factor
s_meanT_relPol_EQ$order <- as.factor(s_meanT_relPol_EQ$order)
s_meanT_relPol_POL$order <- as.factor(s_meanT_relPol_POL$order)

#create new columns with log values (boruta does not accept it...)
s_meanT_relPol_EQ$rangeSize_log10 <- log(s_meanT_relPol_EQ$rangeSize, 10)
s_meanT_relPol_POL$rangeSize_log10 <- log(s_meanT_relPol_POL$rangeSize, 10)

s_meanT_relPol_EQ$nOcc_EQ_log <- log(s_meanT_relPol_EQ$nOcc_EQ)
s_meanT_relPol_POL$nOcc_POL_log <- log(s_meanT_relPol_POL$nOcc_Pol)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.meanT_EQ_relPol <- Boruta(slope_EQ_meanT_relPol ~
                                     rangeSize_log10 +
                                     rangeLoc +
                                     latAmplitude +
                                     roundness +
                                     nOcc_EQ_log +
                                     order +
                                     elevMedian +
                                     bodyMass +
                                     elevAmplitude,
                                data = s_meanT_relPol_EQ, doTrace = 2, ntree = 500,
                                maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.meanT_EQ_relPol)
  setwd(wd_meanT_EQ_relPol)
  write.csv(stats, paste0('meanT_EQ_relPol_', i, '.csv'))
}


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.meanT_POL_relPol <- Boruta(slope_POL_meanT_relPol ~
                                      rangeSize_log10 +
                                      rangeLoc +
                                      latAmplitude +
                                      roundness +
                                      nOcc_POL_log +
                                      order +
                                      elevMedian +
                                      bodyMass +
                                      elevAmplitude,
                            data = s_meanT_relPol_POL, doTrace = 2, ntree = 500,
                            maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.meanT_POL_relPol)
  setwd(wd_meanT_POL_relPol)
  write.csv(stats, paste0('meanT_POL_relPol_', i, '.csv'))
}




###### maxT vs RELATIVE POLEWARDNESS ######

###  MAX  ###

#select only species that had lower correl between vars
s_maxT_relPol <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,]
s_maxT_relPol <- s_maxT_relPol[complete.cases(s_maxT_relPol$Cor_vars_maxT),]

#separate responses
s_maxT_relPol_EQ <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_EQ_maxT_relPol), ]
s_maxT_relPol_POL <- s_maxT_relPol[
  complete.cases(s_maxT_relPol$slope_POL_maxT_relPol), ]

#transform categorical variable in factor
s_maxT_relPol_EQ$order <- as.factor(s_maxT_relPol_EQ$order)
s_maxT_relPol_POL$order <- as.factor(s_maxT_relPol_POL$order)

#create new columns with log values (boruta does not accept it...)
s_maxT_relPol_EQ$rangeSize_log10 <- log(s_maxT_relPol_EQ$rangeSize, 10)
s_maxT_relPol_POL$rangeSize_log10 <- log(s_maxT_relPol_POL$rangeSize, 10)

s_maxT_relPol_EQ$nOcc_EQ_log <- log(s_maxT_relPol_EQ$nOcc_EQ)
s_maxT_relPol_POL$nOcc_POL_log <- log(s_maxT_relPol_POL$nOcc_Pol)



#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.maxT_EQ_relPol <- Boruta(slope_EQ_maxT_relPol ~
                                    rangeSize_log10 +
                                    rangeLoc +
                                    latAmplitude +
                                    roundness +
                                    nOcc_EQ_log +
                                    order +
                                    elevMedian +
                                    bodyMass +
                                    elevAmplitude,
                                  data = s_maxT_relPol_EQ, doTrace = 2, ntree = 500,
                                  maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.maxT_EQ_relPol)
  setwd(wd_maxT_EQ_relPol)
  write.csv(stats, paste0('maxT_EQ_relPol_', i, '.csv'))
}



#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.maxT_POL_relPol <- Boruta(slope_POL_maxT_relPol ~
                                     rangeSize_log10 +
                                     rangeLoc +
                                     latAmplitude +
                                     roundness +
                                     nOcc_POL_log +
                                     order +
                                     elevMedian +
                                     bodyMass +
                                     elevAmplitude,
                                data = s_maxT_relPol_POL, doTrace = 2, ntree = 500,
                                maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.maxT_POL_relPol)
  setwd(wd_maxT_POL_relPol)
  write.csv(stats, paste0('maxT_POL_relPol_', i, '.csv'))
}





###### TEMPERATURE vs CENTRALNESS ######

###  MIN  ###

#select only species that had lower correl between vars
s_minT_distEdge <- slopes[abs(slopes$Cor_vars_minT) <= 0.7,]
s_minT_distEdge <- s_minT_distEdge[complete.cases(s_minT_distEdge$Cor_vars_minT),]
s_minT_distEdge <- s_minT_distEdge[
  complete.cases(s_minT_distEdge$slope_minT_distEdge),]


#transform categorical variable in factor
s_minT_distEdge$order <- as.factor(s_minT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_minT_distEdge$rangeSize_log10 <- log(s_minT_distEdge$rangeSize, 10)
s_minT_distEdge$nOcc_log <- log(s_minT_distEdge$nOcc)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.minT_distEdge <- Boruta(slope_minT_distEdge ~
                                   rangeSize_log10 +
                                   rangeLoc +
                                   latAmplitude +
                                   roundness +
                                   nOcc_log +
                                   order +
                                   elevMedian +
                                   bodyMass +
                                   elevAmplitude,
                                 data = s_minT_distEdge, doTrace = 2, ntree = 500,
                                 maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.minT_distEdge)
  setwd(wd_minT_distEdge)
  write.csv(stats, paste0('minT_distEdge_', i, '.csv'))
}


###  MEAN  ###

#select only species that had lower correl between vars
s_meanT_distEdge <- slopes[abs(slopes$Cor_vars_meanT) <= 0.7,]
s_meanT_distEdge <- s_meanT_distEdge[
  complete.cases(s_meanT_distEdge$Cor_vars_meanT),]
s_meanT_distEdge <- s_meanT_distEdge[
  complete.cases(s_meanT_distEdge$slope_meanT_distEdge),]


#transform categorical variable in factor
s_meanT_distEdge$order <- as.factor(s_meanT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_meanT_distEdge$rangeSize_log10 <- log(s_meanT_distEdge$rangeSize, 10)
s_meanT_distEdge$nOcc_log <- log(s_meanT_distEdge$nOcc)

#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.meanT_distEdge <- Boruta(slope_meanT_distEdge ~
                                    rangeSize_log10 +
                                    rangeLoc +
                                    latAmplitude +
                                    roundness +
                                    nOcc_log +
                                    order +
                                    elevMedian +
                                    bodyMass +
                                    elevAmplitude,
                                  data = s_meanT_distEdge, doTrace = 2, ntree = 500,
                                  maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.meanT_distEdge)
  setwd(wd_meanT_distEdge)
  write.csv(stats, paste0('meanT_distEdge_', i, '.csv'))
}



###  MAX  ###

#select only species that had lower correl between vars
s_maxT_distEdge <- slopes[abs(slopes$Cor_vars_maxT) <= 0.7,]
s_maxT_distEdge <- s_maxT_distEdge[complete.cases(s_maxT_distEdge$Cor_vars_maxT),]
s_maxT_distEdge <- s_maxT_distEdge[
  complete.cases(s_maxT_distEdge$slope_maxT_distEdge),]

#transform categorical variable in factor
s_maxT_distEdge$order <- as.factor(s_maxT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_maxT_distEdge$rangeSize_log10 <- log(s_maxT_distEdge$rangeSize, 10)
s_maxT_distEdge$nOcc_log <- log(s_maxT_distEdge$nOcc)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.maxT_distEdge <- Boruta(slope_maxT_distEdge ~
                                   rangeSize_log10 +
                                   rangeLoc +
                                   latAmplitude +
                                   roundness +
                                   nOcc_log +
                                   order +
                                   elevMedian +
                                   bodyMass +
                                   elevAmplitude,
                                 data = s_maxT_distEdge, doTrace = 2, ntree = 500,
                                 maxRuns = 500)
  
  
  #save results
  stats <- attStats(boruta.maxT_distEdge)
  setwd(wd_maxT_distEdge)
  write.csv(stats, paste0('maxT_distEdge_', i, '.csv'))
}






###### PRECIPITATION vs RELATIVE POLEWARDNESS ######


###### MinPPT vs RELATIVE POLEWARDNESS ######


#select only species that had lower correl between vars
s_minPPT_relPol <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7,]
s_minPPT_relPol <- s_minPPT_relPol[complete.cases(s_minPPT_relPol$Cor_vars_minPPT),]

#separate responses
s_minPPT_relPol_EQ <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_EQ_minPPT_relPol), ]
s_minPPT_relPol_POL <- s_minPPT_relPol[
  complete.cases(s_minPPT_relPol$slope_POL_minPPT_relPol), ]

#transform categorical variable in factor
s_minPPT_relPol_EQ$order <- as.factor(s_minPPT_relPol_EQ$order)
s_minPPT_relPol_POL$order <- as.factor(s_minPPT_relPol_POL$order)

#create new columns with log values (boruta does not accept it...)
s_minPPT_relPol_EQ$rangeSize_log10 <- log(s_minPPT_relPol_EQ$rangeSize, 10)
s_minPPT_relPol_POL$rangeSize_log10 <- log(s_minPPT_relPol_POL$rangeSize, 10)

s_minPPT_relPol_EQ$nOcc_EQ_log <- log(s_minPPT_relPol_EQ$nOcc_EQ)
s_minPPT_relPol_POL$nOcc_POL_log <- log(s_minPPT_relPol_POL$nOcc_Pol)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.minPPT_EQ_relPol <- Boruta(slope_EQ_minPPT_relPol ~
                                    rangeSize_log10 +
                                    rangeLoc +
                                    latAmplitude +
                                    roundness +
                                    nOcc_EQ_log +
                                    order +
                                    elevMedian +
                                    bodyMass +
                                    elevAmplitude,
                                data = s_minPPT_relPol_EQ, doTrace = 2, ntree = 500,
                                maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.minPPT_EQ_relPol)
  setwd(wd_minPPT_EQ_relPol)
  write.csv(stats, paste0('minPPT_EQ_relPol_', i, '.csv'))
}


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.minPPT_POL_relPol <- Boruta(slope_POL_minPPT_relPol ~
                                     rangeSize_log10 +
                                     rangeLoc +
                                     latAmplitude +
                                     roundness +
                                     nOcc_POL_log +
                                     order +
                                     elevMedian +
                                     bodyMass +
                                     elevAmplitude,
                              data = s_minPPT_relPol_POL, doTrace = 2, ntree = 500,
                              maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.minPPT_POL_relPol)
  setwd(wd_minPPT_POL_relPol)
  write.csv(stats, paste0('minPPT_POL_relPol_', i, '.csv'))
}




###### meanPPT vs RELATIVE POLEWARDNESS ######


#select only species that had lower correl between vars
s_meanPPT_relPol <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,]
s_meanPPT_relPol <-
  s_meanPPT_relPol[complete.cases(s_meanPPT_relPol$Cor_vars_meanPPT),]

#separate responses
s_meanPPT_relPol_EQ <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_EQ_meanPPT_relPol), ]
s_meanPPT_relPol_POL <- s_meanPPT_relPol[
  complete.cases(s_meanPPT_relPol$slope_POL_meanPPT_relPol), ]

#transform categorical variable in factor
s_meanPPT_relPol_EQ$order <- as.factor(s_meanPPT_relPol_EQ$order)
s_meanPPT_relPol_POL$order <- as.factor(s_meanPPT_relPol_POL$order)

#create new columns with log values (boruta does not accept it...)
s_meanPPT_relPol_EQ$rangeSize_log10 <- log(s_meanPPT_relPol_EQ$rangeSize, 10)
s_meanPPT_relPol_POL$rangeSize_log10 <- log(s_meanPPT_relPol_POL$rangeSize, 10)

s_meanPPT_relPol_EQ$nOcc_EQ_log <- log(s_meanPPT_relPol_EQ$nOcc_EQ)
s_meanPPT_relPol_POL$nOcc_POL_log <- log(s_meanPPT_relPol_POL$nOcc_Pol)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.meanPPT_EQ_relPol <- Boruta(slope_EQ_meanPPT_relPol ~
                                     rangeSize_log10 +
                                     rangeLoc +
                                     latAmplitude +
                                     roundness +
                                     nOcc_EQ_log +
                                     order +
                                     elevMedian +
                                     bodyMass +
                                     elevAmplitude,
                              data = s_meanPPT_relPol_EQ, doTrace = 2, ntree = 500,
                              maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.meanPPT_EQ_relPol)
  setwd(wd_meanPPT_EQ_relPol)
  write.csv(stats, paste0('meanPPT_EQ_relPol_', i, '.csv'))
}


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.meanPPT_POL_relPol <- Boruta(slope_POL_meanPPT_relPol ~
                                      rangeSize_log10 +
                                      rangeLoc +
                                      latAmplitude +
                                      roundness +
                                      nOcc_POL_log +
                                      order +
                                      elevMedian +
                                      bodyMass +
                                      elevAmplitude,
                              data = s_meanPPT_relPol_POL, doTrace = 2, ntree = 500,
                              maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.meanPPT_POL_relPol)
  setwd(wd_meanPPT_POL_relPol)
  write.csv(stats, paste0('meanPPT_POL_relPol_', i, '.csv'))
}




###### maxPPT vs RELATIVE POLEWARDNESS ######

###  MAX  ###

#select only species that had lower correl between vars
s_maxPPT_relPol <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,]
s_maxPPT_relPol <- s_maxPPT_relPol[complete.cases(s_maxPPT_relPol$Cor_vars_maxPPT),]

#separate responses
s_maxPPT_relPol_EQ <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_EQ_maxPPT_relPol), ]
s_maxPPT_relPol_POL <- s_maxPPT_relPol[
  complete.cases(s_maxPPT_relPol$slope_POL_maxPPT_relPol), ]

#transform categorical variable in factor
s_maxPPT_relPol_EQ$order <- as.factor(s_maxPPT_relPol_EQ$order)
s_maxPPT_relPol_POL$order <- as.factor(s_maxPPT_relPol_POL$order)

#create new columns with log values (boruta does not accept it...)
s_maxPPT_relPol_EQ$rangeSize_log10 <- log(s_maxPPT_relPol_EQ$rangeSize, 10)
s_maxPPT_relPol_POL$rangeSize_log10 <- log(s_maxPPT_relPol_POL$rangeSize, 10)

s_maxPPT_relPol_EQ$nOcc_EQ_log <- log(s_maxPPT_relPol_EQ$nOcc_EQ)
s_maxPPT_relPol_POL$nOcc_POL_log <- log(s_maxPPT_relPol_POL$nOcc_Pol)



#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.maxPPT_EQ_relPol <- Boruta(slope_EQ_maxPPT_relPol ~
                                    rangeSize_log10 +
                                    rangeLoc +
                                    latAmplitude +
                                    roundness +
                                    nOcc_EQ_log +
                                    order +
                                    elevMedian +
                                    bodyMass +
                                    elevAmplitude,
                            data = s_maxPPT_relPol_EQ, doTrace = 2, ntree = 500,
                            maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.maxPPT_EQ_relPol)
  setwd(wd_maxPPT_EQ_relPol)
  write.csv(stats, paste0('maxPPT_EQ_relPol_', i, '.csv'))
}



#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.maxPPT_POL_relPol <- Boruta(slope_POL_maxPPT_relPol ~
                                     rangeSize_log10 +
                                     rangeLoc +
                                     latAmplitude +
                                     roundness +
                                     nOcc_POL_log +
                                     order +
                                     elevMedian +
                                     bodyMass +
                                     elevAmplitude,
                            data = s_maxPPT_relPol_POL, doTrace = 2, ntree = 500,
                            maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.maxPPT_POL_relPol)
  setwd(wd_maxPPT_POL_relPol)
  write.csv(stats, paste0('maxPPT_POL_relPol_', i, '.csv'))
}





###### PRECIPITATION vs CENTRALNESS ######

###  MIN  ###

#select only species that had lower correl between vars
s_minPPT_distEdge <- slopes[abs(slopes$Cor_vars_minPPT) <= 0.7,]
s_minPPT_distEdge <- s_minPPT_distEdge[
  complete.cases(s_minPPT_distEdge$Cor_vars_minPPT),]
s_minPPT_distEdge <- s_minPPT_distEdge[
  complete.cases(s_minPPT_distEdge$slope_minPPT_distEdge),]


#transform categorical variable in factor
s_minPPT_distEdge$order <- as.factor(s_minPPT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_minPPT_distEdge$rangeSize_log10 <- log(s_minPPT_distEdge$rangeSize, 10)
s_minPPT_distEdge$nOcc_log <- log(s_minPPT_distEdge$nOcc)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.minPPT_distEdge <- Boruta(slope_minPPT_distEdge ~
                                   rangeSize_log10 +
                                   rangeLoc +
                                   latAmplitude +
                                   roundness +
                                   nOcc_log +
                                   order +
                                   elevMedian +
                                   bodyMass +
                                   elevAmplitude,
                                 data = s_minPPT_distEdge, doTrace = 2, ntree = 500,
                                 maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.minPPT_distEdge)
  setwd(wd_minPPT_distEdge)
  write.csv(stats, paste0('minPPT_distEdge_', i, '.csv'))
}


###  MEAN  ###

#select only species that had lower correl between vars
s_meanPPT_distEdge <- slopes[abs(slopes$Cor_vars_meanPPT) <= 0.7,]
s_meanPPT_distEdge <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$Cor_vars_meanPPT),]
s_meanPPT_distEdge <- s_meanPPT_distEdge[
  complete.cases(s_meanPPT_distEdge$slope_meanPPT_distEdge),]


#transform categorical variable in factor
s_meanPPT_distEdge$order <- as.factor(s_meanPPT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_meanPPT_distEdge$rangeSize_log10 <- log(s_meanPPT_distEdge$rangeSize, 10)
s_meanPPT_distEdge$nOcc_log <- log(s_meanPPT_distEdge$nOcc)

#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.meanPPT_distEdge <- Boruta(slope_meanPPT_distEdge ~
                                    rangeSize_log10 +
                                    rangeLoc +
                                    latAmplitude +
                                    roundness +
                                    nOcc_log +
                                    order +
                                    elevMedian +
                                    bodyMass +
                                    elevAmplitude,
                                  data = s_meanPPT_distEdge, doTrace = 2, ntree = 500,
                                  maxRuns = 500)
  
  #save results
  stats <- attStats(boruta.meanPPT_distEdge)
  setwd(wd_meanPPT_distEdge)
  write.csv(stats, paste0('meanPPT_distEdge_', i, '.csv'))
}



###  MAX  ###

#select only species that had lower correl between vars
s_maxPPT_distEdge <- slopes[abs(slopes$Cor_vars_maxPPT) <= 0.7,]
s_maxPPT_distEdge <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$Cor_vars_maxPPT),]
s_maxPPT_distEdge <- s_maxPPT_distEdge[
  complete.cases(s_maxPPT_distEdge$slope_maxPPT_distEdge),]

#transform categorical variable in factor
s_maxPPT_distEdge$order <- as.factor(s_maxPPT_distEdge$order)

#create new columns with log values (boruta does not accept it...)
s_maxPPT_distEdge$rangeSize_log10 <- log(s_maxPPT_distEdge$rangeSize, 10)
s_maxPPT_distEdge$nOcc_log <- log(s_maxPPT_distEdge$nOcc)


#run RF through boruta 5 times
for(i in 1:5)
{
  boruta.maxPPT_distEdge <- Boruta(slope_maxPPT_distEdge ~
                                   rangeSize_log10 +
                                   rangeLoc +
                                   latAmplitude +
                                   roundness +
                                   nOcc_log +
                                   order +
                                   elevMedian +
                                   bodyMass +
                                   elevAmplitude,
                                 data = s_maxPPT_distEdge, doTrace = 2, ntree = 500,
                                 maxRuns = 500)
  
  
  #save results
  stats <- attStats(boruta.maxPPT_distEdge)
  setwd(wd_maxPPT_distEdge)
  write.csv(stats, paste0('maxPPT_distEdge_', i, '.csv'))
}




##### Prepare tables for figure

#base WD
WD_base <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Boruta'

#list folders
folders <- list.dirs(WD_base, recursive = FALSE)

#get names
basename(folders)

#define row names for temperature table
temp_rows <- c('minT relpol EQ',
               'minT relpol POL',
               'meanT relpol EQ',
               'meanT relpol POL',
               'maxT relpol EQ',
               'maxT relpol POL',
               '',
               '',
               'minT distEdge',
               'meanT distEdge',
               'maxT distEdge')

#define row names for precipitation table
ppt_rows <- c('minPPT relpol EQ',
              'minPPT relpol POL',
              'meanPPT relpol EQ',
              'meanPPT relpol POL',
              'maxPPT relpol EQ',
              'maxPPT relpol POL',
              '',
              '',
              'minPPT distEdge',
              'meanPPT distEdge',
              'maxPPT distEdge')

#create empty temperature table
boruta_temp_table <- matrix(NA,
                            nrow=length(temp_rows),
                            ncol=length(pred_internal))

#add row and column names to temperature table
rownames(boruta_temp_table) <- temp_rows
colnames(boruta_temp_table) <- pred_labels

#create empty precipitation table
boruta_ppt_table <- matrix(NA,
                           nrow=length(ppt_rows),
                           ncol=length(pred_internal))

#add row and column names to precipitation table
rownames(boruta_ppt_table) <- ppt_rows
colnames(boruta_ppt_table) <- pred_labels

#loop through temperature Boruta folders
for(i in 1:length(folders))
{
  #select current folder
  current_folder <- folders[i]
  
  #get folder name only
  folder_name <- basename(current_folder)
  
  #skip non-temperature folders
  if(!grepl('^(minT|meanT|maxT)_', folder_name)) next
  
  #list csv files inside folder
  csv_files <- list.files(current_folder,
                          pattern='\\.csv$',
                          full.names=TRUE)
  
  #create empty matrix to store decisions
  decision_matrix <- matrix(NA,
                            nrow=length(pred_internal),
                            ncol=length(csv_files))
  
  #add predictor names as row names
  rownames(decision_matrix) <- pred_internal
  
  #loop through Boruta repetitions
  for(j in 1:length(csv_files))
  {
    #read current Boruta result
    x <- read.csv(csv_files[j], row.names=1)
    
    #extract decision column
    d <- x$decision
    
    #add predictor names to decision vector
    names(d) <- rownames(x)
    
    #rename number-of-records predictor to unified name
    names(d)[grepl('^nOcc', names(d))] <- 'nOcc_EQ_log'
    
    #convert Boruta decisions to numeric values
    d[d=='Confirmed'] <- 1
    d[d=='Rejected'] <- -1
    d[d=='Tentative'] <- 0
    d <- as.numeric(d)
    
    #add predictor names again after numeric conversion
    names(d) <- rownames(x)
    
    #rename number-of-records predictor to unified name again
    names(d)[grepl('^nOcc', names(d))] <- 'nOcc_EQ_log'
    
    #store decisions in matrix
    decision_matrix[,j] <- d[pred_internal]
  }
  
  #sum decisions across repetitions
  summed_decisions <- rowSums(decision_matrix)
  
  #convert folder name to table row name
  row_name <- gsub('_EQ_relPol',' relpol EQ',folder_name)
  row_name <- gsub('_POL_relPol',' relpol POL',row_name)
  row_name <- gsub('_distEdge',' distEdge',row_name)
  
  #insert results into temperature table
  boruta_temp_table[row_name,] <- summed_decisions
}


#loop through precipitation Boruta folders
for(i in 1:length(folders))
{
  #select current folder
  current_folder <- folders[i]
  
  #get folder name only
  folder_name <- basename(current_folder)
  
  #skip non-precipitation folders
  if(!grepl('^(minPPT|meanPPT|maxPPT)_', folder_name)) next
  
  #list csv files inside folder
  csv_files <- list.files(current_folder,
                          pattern='\\.csv$',
                          full.names=TRUE)
  
  #create empty matrix to store decisions
  decision_matrix <- matrix(NA,
                            nrow=length(pred_internal),
                            ncol=length(csv_files))
  
  #add predictor names as row names
  rownames(decision_matrix) <- pred_internal
  
  #loop through Boruta repetitions
  for(j in 1:length(csv_files))
  {
    #read current Boruta result
    x <- read.csv(csv_files[j], row.names=1)
    
    #extract decision column
    d <- x$decision
    
    #add predictor names to decision vector
    names(d) <- rownames(x)
    
    #rename number-of-records predictor to unified name
    names(d)[grepl('^nOcc', names(d))] <- 'nOcc_EQ_log'
    
    #convert Boruta decisions to numeric values
    d[d=='Confirmed'] <- 1
    d[d=='Rejected'] <- -1
    d[d=='Tentative'] <- 0
    d <- as.numeric(d)
    
    #add predictor names again after numeric conversion
    names(d) <- rownames(x)
    
    #rename number-of-records predictor to unified name again
    names(d)[grepl('^nOcc', names(d))] <- 'nOcc_EQ_log'
    
    #store decisions in matrix
    decision_matrix[,j] <- d[pred_internal]
  }
  
  #sum decisions across repetitions
  summed_decisions <- rowSums(decision_matrix)
  
  #convert folder name to table row name
  row_name <- gsub('_EQ_relPol',' relpol EQ',folder_name)
  row_name <- gsub('_POL_relPol',' relpol POL',row_name)
  row_name <- gsub('_distEdge',' distEdge',row_name)
  
  #insert results into precipitation table
  boruta_ppt_table[row_name,] <- summed_decisions
}


#save tables
setwd(WD_base)

#save temperature Boruta summary table as csv file
write.csv(boruta_temp_table,
          file='boruta_temperature_summary.csv',
          row.names=TRUE)

#save precipitation Boruta summary table as csv file
write.csv(boruta_ppt_table,
          file='boruta_precipitation_summary.csv',
          row.names=TRUE)
