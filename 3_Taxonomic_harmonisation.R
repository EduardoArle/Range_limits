####  This script deals with the taxonomic harmonisation of the names in IUCN
####  range maps and in the body mass database (Smith_etal_2003_Ecology), using
####  the GBIF backbone as reference

#load packages 
library(rgbif)

#list WDs
wd_ranges <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps'
wd_mass <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Body mass/COMBINE_archives'
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


########### COMBINE #############

#setwd to trait table
setwd(wd_mass)
traits <- read.csv('trait_data_imputed.csv')

#keep only the columns you need
traits_mass <- traits[, c('iucn2020_binomial', 'genus', 'adult_mass_g')]

#drop duplicated species if present
traits_mass <- traits_mass[!duplicated(traits_mass$iucn2020_binomial), ]

#build crosswalk table
zref <- data.frame(combine_name = traits_mass$iucn2020_binomial,
                   gbif_name = NA_character_,
                   canonicalName = NA_character_,
                   status = NA_character_,
                   matchType = NA_character_,
                   confidence = NA_real_,
                   usageKey = NA_integer_,
                   stringsAsFactors = FALSE)

#loop through COMBINE names
for(i in 1:nrow(zref))
{
  combine_name <- zref$combine_name[i]
  
  bb <- tryCatch(name_backbone(name = combine_name),
                 error = function(e) NULL)
  
  if(is.null(bb)){
    zref$gbif_name[i] <- NA_character_
    zref$canonicalName[i] <- NA_character_
    zref$status[i] <- NA_character_
    zref$matchType[i] <- NA_character_
    zref$confidence[i] <- NA_real_
    zref$usageKey[i] <- NA_integer_
    print(i)
    next
  }
  
  if(!is.null(bb$species) && !is.na(bb$species) && nzchar(bb$species)){
    zref$gbif_name[i] <- bb$species
  } else if(!is.null(bb$canonicalName) && !is.na(bb$canonicalName) &&
            nzchar(bb$canonicalName)){
    zref$gbif_name[i] <- bb$canonicalName
  } else {
    zref$gbif_name[i] <- NA_character_
  }
  
  if(!is.null(bb$canonicalName) && !is.na(bb$canonicalName) &&
     nzchar(bb$canonicalName)){
    zref$canonicalName[i] <- bb$canonicalName
  } else {
    zref$canonicalName[i] <- NA_character_
  }
  
  if(!is.null(bb$status) && !is.na(bb$status)){
    zref$status[i] <- bb$status
  } else {
    zref$status[i] <- NA_character_
  }
  
  if(!is.null(bb$matchType) && !is.na(bb$matchType)){
    zref$matchType[i] <- bb$matchType
  } else {
    zref$matchType[i] <- NA_character_
  }
  
  if(!is.null(bb$confidence) && !is.na(bb$confidence)){
    zref$confidence[i] <- bb$confidence
  } else {
    zref$confidence[i] <- NA_real_
  }
  
  if(!is.null(bb$usageKey) && !is.na(bb$usageKey)){
    zref$usageKey[i] <- bb$usageKey
  } else {
    zref$usageKey[i] <- NA_integer_
  }
  
  print(i)
}


#join body mass table with harmonisation
traits_mass2 <- merge(traits_mass, zref, by.x = 'iucn2020_binomial',
                      by.y = 'combine_name', all.x = TRUE)

#create provenance column
traits_mass2$mass_source <- NA_character_

#loop through species
for(i in 1:nrow(traits_mass2))
{
  #if species already has mass, keep it and flag as direct
  if(!is.na(traits_mass2$adult_mass_g[i])){
    traits_mass2$mass_source[i] <- "species"
    
  } else {
    
    #get genus of focal species
    gen <- traits_mass2$genus[i]
    
    #get masses of other species in same genus
    gen_masses <- traits_mass2$adult_mass_g[traits_mass2$genus == gen]
    gen_masses <- gen_masses[!is.na(gen_masses)]
    
    #if genus has available masses, fill with genus mean
    if(length(gen_masses) > 0){
      traits_mass2$adult_mass_g[i] <- mean(gen_masses)
      traits_mass2$mass_source[i] <- "genus_mean"
      
    } else {
      #if no mass available in genus, keep NA
      traits_mass2$mass_source[i] <- "no_info"
    }
  }
}

#save
setwd(wd_mass)
write.csv(traits_mass2, 'Harmonised_COMBINE_bodymass.csv', row.names = FALSE)
