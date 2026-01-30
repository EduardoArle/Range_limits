####  This script deals with the taxonomic harmonisation of the names in IUCN
####  range maps and in the body mass database (Smith_etal_2003_Ecology), using
####  the GBIF backbone as reference

#load packages 
library(rgbif)

#list WDs
wd_ranges <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps'
wd_mass <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Mammal_trait_data/Smith_etal_2003_Ecology'
wd_out <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review'

########### IUCN #############

#make a species list
setwd(wd_ranges)
iucn_names <- gsub('.shp', '', list.files(pattern = '.shp$'))

#build crosswalk table to be populated in the for loop
xref <- data.frame(iucn_name = iucn_names, gbif_name = NA_character_,
                   canonicalName = NA_character_, status = NA_character_,
                   matchType = NA_character_, confidence = NA_real_,
                   usageKey = NA_integer_, order = NA_character_,
                   family = NA_character_, stringsAsFactors = FALSE)

#loop though species getting the harmonised names
for(i in 1:nrow(xref))
{
  #get IUCN species name
  iucn_name <- xref$iucn_name[i]
  
  #get GBIF backbone info (do not crash on timeouts/errors)
  bb <- tryCatch(name_backbone(name = iucn_name),
                 error = function(e) NULL)
  
  #defaults (so missing fields never crash)
  xref$gbif_name[i] <- NA_character_
  xref$canonicalName[i] <- NA_character_
  xref$status[i] <- NA_character_
  xref$matchType[i] <- NA_character_
  xref$confidence[i] <- NA_real_
  xref$usageKey[i] <- NA_integer_
  xref$order[i] <- NA_character_
  xref$family[i] <- NA_character_
  
  #if query failed, move on
  if(is.null(bb)){
    print(i)
    next
  }
  
  #populate GBIF name (species if available, otherwise canonicalName)
  if(!is.null(bb$species) && length(bb$species) > 0 &&
     !is.na(bb$species) && nzchar(bb$species)){
    xref$gbif_name[i] <- bb$species
  } else if(!is.null(bb$canonicalName) && length(bb$canonicalName) > 0 &&
            !is.na(bb$canonicalName) && nzchar(bb$canonicalName)){
    xref$gbif_name[i] <- bb$canonicalName
  }
  
  #populate canonicalName
  if(!is.null(bb$canonicalName) && length(bb$canonicalName) > 0 &&
     !is.na(bb$canonicalName) && nzchar(bb$canonicalName)){
    xref$canonicalName[i] <- bb$canonicalName
  }
  
  #populate status
  if(!is.null(bb$status) && length(bb$status) > 0 &&
     !is.na(bb$status) && nzchar(bb$status)){
    xref$status[i] <- bb$status
  }
  
  #populate matchType
  if(!is.null(bb$matchType) && length(bb$matchType) > 0 &&
     !is.na(bb$matchType) && nzchar(bb$matchType)){
    xref$matchType[i] <- bb$matchType
  }
  
  #populate confidence
  if(!is.null(bb$confidence) && length(bb$confidence) > 0 &&
     !is.na(bb$confidence)){
    xref$confidence[i] <- bb$confidence
  }
  
  #populate usageKey
  if(!is.null(bb$usageKey) && length(bb$usageKey) > 0 &&
     !is.na(bb$usageKey)){
    xref$usageKey[i] <- bb$usageKey
  }
  
  #populate order
  if(!is.null(bb$order) && length(bb$order) > 0 &&
     !is.na(bb$order) && nzchar(bb$order)){
    xref$order[i] <- bb$order
  }
  
  #populate family
  if(!is.null(bb$family) && length(bb$family) > 0 &&
     !is.na(bb$family) && nzchar(bb$family)){
    xref$family[i] <- bb$family
  }
  
  print(i)
}


#save harmonisation table
setwd(wd_out)
write.csv(xref, 'Harmonised_table.csv', row.names = F)


########### Smith et al. 2003 #############

#load table with species list
setwd(wd_mass)
mass <- read.csv('Mammals_bodymass_Smmith_2003.csv')

#build crosswalk table to be populated in the for loop
yref <- data.frame(smith_name = mass$Species, gbif_name = NA_character_,
                   canonicalName = NA_character_, status = NA_character_,
                   matchType = NA_character_, confidence = NA_real_,
                   usageKey = NA_integer_, stringsAsFactors = FALSE)

#loop though species getting the harmonised names
for(i in 1:nrow(yref))
{
  #get Smith species name
  smith_name <- yref$smith_name[i]
  
  #get GBIF backbone info (do not crash on timeouts/errors)
  bb <- tryCatch(name_backbone(name = smith_name),
                 error = function(e) NULL)
  
  #if query failed, populate NAs and move on
  if(is.null(bb)){
    yref$gbif_name[i] <- NA_character_
    yref$canonicalName[i] <- NA_character_
    yref$status[i] <- NA_character_
    yref$matchType[i] <- NA_character_
    yref$confidence[i] <- NA_real_
    yref$usageKey[i] <- NA_integer_
    print(i)
    next
  }
  
  #populate GBIF name (species if available, otherwise canonicalName)
  gbif_name <- NA_character_
  if(!is.null(bb$species) && !is.na(bb$species) && nzchar(bb$species)){
    gbif_name <- bb$species
  } else if(!is.null(bb$canonicalName) &&
            !is.na(bb$canonicalName) &&
            nzchar(bb$canonicalName)){
    gbif_name <- bb$canonicalName
  }
  yref$gbif_name[i] <- gbif_name
  
  #populate canonicalName
  if(!is.null(bb$canonicalName) &&
     !is.na(bb$canonicalName) &&
     nzchar(bb$canonicalName)){
    yref$canonicalName[i] <- bb$canonicalName
  } else {
    yref$canonicalName[i] <- NA_character_
  }
  
  #populate status
  if(!is.null(bb$status) && !is.na(bb$status)){
    yref$status[i] <- bb$status
  } else {
    yref$status[i] <- NA_character_
  }
  
  #populate matchType
  if(!is.null(bb$matchType) && !is.na(bb$matchType)){
    yref$matchType[i] <- bb$matchType
  } else {
    yref$matchType[i] <- NA_character_
  }
  
  #populate confidence
  if(!is.null(bb$confidence) && !is.na(bb$confidence)){
    yref$confidence[i] <- bb$confidence
  } else {
    yref$confidence[i] <- NA_real_
  }
  
  #populate usageKey
  if(!is.null(bb$usageKey) && !is.na(bb$usageKey)){
    yref$usageKey[i] <- bb$usageKey
  } else {
    yref$usageKey[i] <- NA_integer_
  }
  
  print(i)
}

#join tables
mass2 <- merge(mass, yref, by.x = "Species", by.y = "smith_name", all.x = TRUE)

#save table
setwd(wd_mass)
write.csv(mass2, 'Harmonised_Mammals_bodymass_Smmith_2003.csv')
