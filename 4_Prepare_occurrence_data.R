#load packages 
library(terra); library(sf); library(data.table); library(raster)

#list WDs
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review'
wd_raw <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Raw data'
wd_variables <- '/Users/carloseduardoaribeiro/Documents/CNA/Data/Variables/wc2-5'
wd_PA_thinned <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Thinned Occ'

#check which orders we have represented in the selected sps list
setwd(wd_tables)
sps_list <- read.csv('Harmonised_table.csv')
orders <- unique(sps_list$order)
sort(orders)

#manually download each order occurrences from GBIF

#load one variable to thin the records by the same resolution
setwd(wd_variables)
var <- raster('bio2.bil')

#make ID raster to thin points by variable resolution
ID_raster <- var
ID_raster[] <- c(1:length(ID_raster)) 


##########################
#       RODENTIA         #
##########################

#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Rodentia.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Rodentia.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Rodentia_thin.csv', row.names = F) 
 

##########################
#      CHIROPTERA        #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Chiroptera.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Chiroptera.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Chiroptera_thin.csv', row.names = F) 



##########################
#       CARNIVORA        #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Carnivora.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Carnivora.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Carnivora_thin.csv', row.names = F) 



##########################
#       ARTIODACTYLA     #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Artiodactyla.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Artiodactyla.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Artiodactyla_thin.csv', row.names = F) 



##########################
#      DIPROTODONTIA     #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Diprotodontia.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Diprotodontia.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Diprotodontia_thin.csv', row.names = F) 



##########################
#        LAGOMORPHA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Lagomorpha.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Lagomorpha.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Lagomorpha_thin.csv', row.names = F) 



##########################
#      SORICOMORPHA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Soricomorpha.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Soricomorpha.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Soricomorpha_thin.csv', row.names = F) 



##########################
#    ERINACEPMORPHA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Erinaceomorpha.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Erinaceomorpha.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Erinaceomorpha_thin.csv', row.names = F) 



##########################
#        PRIMATES        #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Primates.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Primates.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Primates_thin.csv', row.names = F) 



##########################
#     PERAMELEMORPHIA    #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Peramelemorphia.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Peramelemorphia.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Peramelemorphia_thin.csv', row.names = F) 



##########################
#     DASYUROMORPHIA     #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Dasyuromorphia.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Dasyuromorphia.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Dasyuromorphia_thin.csv', row.names = F) 




##########################
#      PERISSODACTLA     #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Perissodactyla.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Perissodactyla.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Perissodactyla_thin.csv', row.names = F) 



##########################
#       MONOTREMATA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Monotremata.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Monotremata.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Monotremata_thin.csv', row.names = F) 



##########################
#     DIDELPHIMORPHIA    #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Didelphimorphia.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Didelphimorphia.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Didelphimorphia_thin.csv', row.names = F) 



##########################
#        CINGULATA       #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Cingulata.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Cingulata.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Cingulata_thin.csv', row.names = F) 



##########################
#      PROBOSCIDEA       #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Proboscidea.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Proboscidea.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Proboscidea_thin.csv', row.names = F) 



##########################
#         PILOSA         #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Pilosa.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Pilosa.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Pilosa_thin.csv', row.names = F) 



##########################
#      HYRACOIDEA       #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Hyracoidea.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Hyracoidea.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Hyracoidea_thin.csv', row.names = F) 



##########################
#      AFROSORICIDA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Afrosoricida.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Afrosoricida.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Afrosoricida_thin.csv', row.names = F) 



##########################
#      MACROSCELIDEA     #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Macroscelidea.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Macroscelidea.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(order_occ_sf),
                         by = c('cellID', 'species'))


#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(no_reps, 'Macroscelidea_thin.csv', row.names = F) 



##########################
#        SCANDENTIA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Scandentia.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Scandentia.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Scandentia_thin.csv', row.names = F) 



##########################
#        DERMOPTERA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Dermoptera.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Dermoptera.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Dermoptera_thin.csv', row.names = F) 



##########################
#    MICROBIOTHERIA      #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Microbiotheria.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Microbiotheria.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))

#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(order_occ_thin, 'Microbiotheria_thin.csv', row.names = F) 



##########################
#        PHOLIDOTA       #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Pholidota.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Pholidota.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))


#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(no_reps, 'Pholidota_thin.csv', row.names = F) 



##########################
#      TUBULIDENTATA     #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Tubulidentata.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Tubulidentata.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))


#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(no_reps, 'Tubulidentata_thin.csv', row.names = F) 



##########################
#    PAUCITUBERCULATA    #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Paucituberculata.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Paucituberculata.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))


#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(no_reps, 'Paucituberculata_thin.csv', row.names = F) 



##########################
#    NOTORYCTEMORPHIA    #
##########################



#inspect the headers of the table
setwd(wd_raw)
cols <- names(read.delim("Notoryctemorphia.csv", sep = "\t", nrows = 1))
cols

#define columns to keep
keep <- c("gbifID", "scientificName", "species", "order", "decimalLatitude",  
          "decimalLongitude", "basisOfRecord")

#build colClasses
colClasses <- rep("NULL", length(cols))

colClasses[match("gbifID", cols)]           <- "character"
colClasses[match("scientificName", cols)]   <- "character"
colClasses[match("species", cols)]          <- "character"
colClasses[match("order", cols)]            <- "character"
colClasses[match("decimalLatitude", cols)]  <- "numeric"
colClasses[match("decimalLongitude", cols)] <- "numeric"
colClasses[match("basisOfRecord", cols)]    <- "character"

#load full table
setwd(wd_raw)

all_ord <- read.delim("Notoryctemorphia.csv", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      colClasses = colClasses)


#eliminate rows with same "species", "decimalLatitude", and "decimalLongitude"
no_reps <- unique(as.data.table(all_ord),
                  by = c('species', 'decimalLatitude', 'decimalLongitude'))

#create spatial object
order_occ_sf <- st_as_sf(no_reps, 
                         coords = c('decimalLongitude', 'decimalLatitude'),
                         crs = st_crs(var))

#get cellID value to thin the records to max one per sps per grid cel
cellIDs <- extract(ID_raster, order_occ_sf)

#include cell value and coordinates into species data
no_reps$cellID <- cellIDs

#keep only one entry per cell
order_occ_thin <- unique(as.data.table(no_reps),
                         by = c('cellID', 'species'))


#save just in case it crashes
setwd(wd_PA_thinned)
write.csv(no_reps, 'Notoryctemorphia_thin.csv', row.names = F) 
