####  This script plots the ranges of selected species showing geographical and taxonomic distribution os our sample

#load libraries
library(sf); library(terra); library(plyr); library(rnaturalearth)

#list wds
wd_tables <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Tables'
wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"
wd_thinned_occ <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Occurrences/Thinned Occ'
wd_vars <- '/Users/carloseduardoaribeiro/Documents/Post-doc/Variable layes/BioClim_layers'
wd_map_stuff <- '/Users/carloseduardoaribeiro/Documents/Map stuff/Carstens_world_proj'
wd_tab_SI <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/SI/Tables'
wd_mdd <- '/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Manuscript/Submission NEE/Review/Figures/Fig2/Mammal Diversity Database'


#disable s2 geometry engine to avoid topology errors in range polygons
sf::sf_use_s2(FALSE)

#load table with info on sps ranges
setwd(wd_tables)
sps_table <- read.csv('Selected_species.csv')

#load elevation raster (just for mapping)
setwd(wd_vars)
elev <- rast('wc2.1_2.5m_elev.tif')

#list species
setwd(wd_tab_SI)
tab <- read.csv('Table_S1.csv')
sps_list <- tab$sps

#count species number in orders represented in the study
ord_tab <- ddply(tab, .(order), function(x)
  {data.frame(n_study = length(unique(x$sps)))})

#loas mammal diversity database table
setwd(wd_mdd)
mdd <- read.csv("MDD_v2.4_6871species.csv")

#keep only wild extant mammals
mdd_wild <- mdd[which(mdd$extinct == 0 & mdd$domestic == 0),]

#remove pinnipeds
pinniped_families <- c("Phocidae","Otariidae","Odobenidae")
mdd_terr <- mdd_wild[!(mdd_wild$family %in% pinniped_families),]

#count species per order
mdd_counts <- as.data.frame(table(mdd_terr$order), stringsAsFactors = FALSE)
names(mdd_counts) <- c("order_mdd","n_total")
mdd_counts$n_total <- as.numeric(mdd_counts$n_total)
mdd_counts[order(mdd_counts$n_total, decreasing = TRUE),]

#add a column to fill later with total richness per order
ord_tab$n_total <- NA_real_

#fill total richness per order using mdd counts
ord_tab$n_total <- mdd_counts$n_total[match(ord_tab$order, mdd_counts$order_mdd)]

#map gbif Soricomorpha to mdd Eulipotyphla
ord_tab$n_total[ord_tab$order == "Soricomorpha"] <-
  mdd_counts$n_total[mdd_counts$order_mdd == "Eulipotyphla"]

#sort by richness in the study
ord_tab <- ord_tab[order(ord_tab$n_study, decreasing = TRUE),]

#calculate percentage of species represented in the study
ord_tab$perc_study <- round((ord_tab$n_study / ord_tab$n_total) * 100, 0)

#create colour palette for orders
colours <- c("#6A3D9A","#FB9A99","#FFED6F","#33A02C","#FF7F00","#B15928",
             "#A6CEE3","#B2DF8A","#FF00EE","#FDBF6F","#CAB2D6","#FFFF99",
             "#8DD3C7","#80B1D3","#FDB462","#B3DE69","#1F78B4","#BC80BD",
             "#E31A1C","#CCEBC5")

#order for plotting (smallest to largest so largest appears on top)
ord <- order(ord_tab$n_study, decreasing = FALSE)
counts <- ord_tab$n_study[ord]
orders <- ord_tab$order[ord]
perc <- ord_tab$perc_study[ord]

#assign colours to orders
col_map <- setNames(colours, ord_tab$order)

#retrieve colours in plotting order
cols <- col_map[orders]



##### PLOT SPECIES PER TAXON #####



#set margins for plot
par(mar = c(5,12,3,2), cex = 1.2)

#bar plot species per order (study sample only)
bp <- barplot(counts, names.arg = orders, horiz = TRUE, las = 1, col = cols,
              xlim = c(0,560), xaxt = 'n', xlab = 'Number of species',
              space = 0.25, ylab = NA, main = 'Species per order',
              cex.names = 1.1)

axis(1, at = c(0,100,200,300,400,500))

#plot counts (right aligned)
text(480, bp, labels = counts, adj = 1, cex = 1)

#plot percentages (right aligned)
text(530, bp, labels = paste0('(',perc,'%)'), adj = 1, cex = 1)

#save 1000

##### PLOT LARGE MAP #####



#make table relating order to colour code (following barplot scheme)
ord_col <- data.frame(order = ord_tab$order, colour = colours)

#include colour code for each species based on order
tab_2 <- merge(tab, ord_col, by = 'order', sort = FALSE)

#load map stuff
setwd(wd_map_stuff)
world <- readRDS('wrld.rds')
worldmapframe <- readRDS('Worldmapframe.rds')

#transform SpatialPolygons into sf object
worldmapframe <- st_as_sf(worldmapframe)

#change the crs of the frame and the elev raster
worldmapframe <- st_transform(worldmapframe, crs = st_crs(world))
elev <- project(elev, world)

#set margins for plot
oldpar <- par(mar = c(5,6,3,3))

#plot world
plot(elev, box = NA, legend = NA, axes = NA, col = 'grey90')

#add frame
plot(worldmapframe, add = TRUE)

#randomise plotting order so colours do not appear in blocks
sps_list <- sample(tab$sps)

for(i in 1:length(sps_list))
{
  #select species
  sps <- sps_list[i]
  
  #load species range map
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps)
  
  #select the range representing only the native range of the species
  range <- range[!range$legend %in% c('Introduced','Vagrant'), ]
  
  #get area of main fragment of the range
  dom_frag_km2 <- as.numeric(strsplit(sps_table$range_frag_km2[
    match(sps, sps_table$species)], ";")[[1]][
      sps_table$dom_frag[match(sps, sps_table$species)]])
  
  #split into fragments (one POLYGON per fragment)
  frags <- st_cast(range, 'POLYGON')
  
  #project to a local equal-area CRS so areas are meaningful
  cen <- st_coordinates(st_centroid(st_union(frags)))[1, ]
  crs_laea <- paste0('+proj=laea +lat_0=', cen[2], ' +lon_0=', cen[1],
                     ' +datum=WGS84 +units=m +no_defs')
  frags_m <- st_transform(frags, crs_laea)
  
  #areas in km2
  a_km2 <- as.numeric(st_area(frags_m)) / 1e6
  
  #pick the fragment whose area is closest to dom_frag_km2
  frag_pos <- which.min(abs(a_km2 - dom_frag_km2))
  
  #keep only the dominant fragment
  range_dom <- frags[frag_pos, ]
  
  #transform to the map CRS for plotting
  range_dom <- st_transform(range_dom, st_crs(world))
  
  #plot species range
  col_i <- tab_2$colour[match(sps, tab_2$sps)]
  plot(range_dom, add = TRUE,
       col = paste0(col_i, '20'),
       border = paste0(col_i, '30'))
  
  print(i)
}



##### PLOT SMALL MAP WITH SPECIES PER CONTINENT #####



#load world polygons
world_cont <- ne_countries(scale = 'medium', returnclass = 'sf')

#filter continents where we have species
keep <- c('North America','South America','Europe','Africa','Asia','Oceania')
world_cont <- world_cont[world_cont$continent %in% keep, ]

#dissolve countries into continent polygons:
continents <- aggregate(world_cont['continent'],
                        by = list(world_cont$continent), FUN = length)

names(continents)[1] <- 'continent'

#make CRS matche my ranges
continents <- st_transform(continents, st_crs(world))

#prepare vector to store results
sps_cont <- vector('list', length(tab$sps))
names(sps_cont) <- tab$sps

#loop through species to find continent(s) of each
for(i in 1:length(tab$sps))
{
  #get species name
  sps <- tab$sps[i]
  
  #load species range
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps, quiet = TRUE)
  
  #remove introduced and vagrant occurrences
  range <- range[!range$legend %in% c('Introduced','Vagrant'), ]
  
  #get area of main fragment of the range
  dom_frag_km2 <- as.numeric(strsplit(sps_table$range_frag_km2[
    match(sps, sps_table$species)], ';')[[1]][
      sps_table$dom_frag[match(sps, sps_table$species)]])
  
  #split range into fragments
  frags <- st_cast(range, 'POLYGON')
  
  #get centroid to build local equal-area projection
  cen <- st_coordinates(st_centroid(st_union(frags)))[1, ]
  
  #define local equal-area projection
  crs_laea <- paste0('+proj=laea +lat_0=', cen[2], ' +lon_0=', cen[1],
                     ' +datum=WGS84 +units=m +no_defs')
  
  #project fragments to equal-area CRS
  frags_m <- st_transform(frags, crs_laea)
  
  #calculate fragment areas in km2
  a_km2 <- as.numeric(st_area(frags_m)) / 1e6
  
  #identify fragment closest to dominant fragment area
  frag_pos <- which.min(abs(a_km2 - dom_frag_km2))
  
  #keep only dominant fragment
  range_dom <- frags[frag_pos, ]
  
  #transform dominant fragment to continent CRS
  range_dom <- st_transform(range_dom, st_crs(continents))
  
  #identify continent intersection
  inter <- st_intersects(range_dom, continents)
  
  #store continent names
  sps_cont[[i]] <- continents$continent[inter[[1]]]
  
  #print progress
  print(i)
}



## Convert the list sps_cont into a table.

#create empty table
cont_table <- data.frame()

#expand species–continent combinations
for(i in 1:length(sps_cont))
{
  if(length(sps_cont[[i]]) > 0)
  {
    tmp <- data.frame(species = tab$sps[i],
                      continent = sps_cont[[i]])
    cont_table <- rbind(cont_table, tmp)
  }
}

#count unique species per continent
cont_counts <- aggregate(cont_table$species,
                         by = list(cont_table$continent),
                         FUN = function(x) length(unique(x)))

names(cont_counts) <- c('continent','n_species')


#add species counts to continent polygons
continents$n_species <- cont_counts$n_species[
  match(continents$continent, cont_counts$continent)]

#assign colours to continents
cont_cols <- c('North America' = '#1F78B4', 'South America' = '#33A02C',
               'Europe' = '#6A3D9A', 'Africa' = '#FF7F00',
               'Asia' = '#E31A1C', 'Oceania' = '#B15928')

#get colours in polygon order
cols_cont <- cont_cols[continents$continent]

#plot continents with chosen colours
plot(st_geometry(continents),
     col = cols_cont,
     border = 'NA',
     axes = FALSE)

#add frame to the map
plot(worldmapframe, add = TRUE)

#get label positions inside polygons
lab_pts <- st_point_on_surface(continents)

#extract coordinates
coords <- st_coordinates(lab_pts)

#North America (+10° north)
coords[continents$continent == "North America",2] <-
  coords[continents$continent == "North America",2] + 1110000

#South America (+10° north)
coords[continents$continent == "South America",2] <-
  coords[continents$continent == "South America",2] + 1110000

#Africa (+15° north, 10° west)
coords[continents$continent == "Africa",2] <-
  coords[continents$continent == "Africa",2] + 1665000

coords[continents$continent == "Africa",1] <-
  coords[continents$continent == "Africa",1] - 1110000

#Asia (+15° north)
coords[continents$continent == "Asia",2] <-
  coords[continents$continent == "Asia",2] + 1110000

#add species numbers
text(coords,
     labels = continents$n_species,
     col = "white",
     font = 2,
     cex = 1.4)




##### PLOT HISTOGRAM WITH RANGE SIZES #####


#vector to store dominant range sizes
range_sizes <- numeric(length(sps_list))

for(i in 1:length(sps_list))
{
  sps <- sps_list[i]
  
  #load range
  setwd(wd_ranges)
  range <- st_read(dsn = wd_ranges, layer = sps, quiet = TRUE)
  
  #remove introduced and vagrant
  range <- range[!range$legend %in% c('Introduced','Vagrant'), ]
  
  #split fragments
  frags <- st_cast(range, "POLYGON")
  
  #centroid for equal-area projection
  cen <- st_coordinates(st_centroid(st_union(frags)))[1, ]
  
  crs_laea <- paste0("+proj=laea +lat_0=", cen[2],
                     " +lon_0=", cen[1],
                     " +datum=WGS84 +units=m +no_defs")
  
  #project fragments
  frags_m <- st_transform(frags, crs_laea)
  
  #areas in km²
  a_km2 <- as.numeric(st_area(frags_m)) / 1e6
  
  #dominant fragment
  dom_frag_km2 <- as.numeric(strsplit(
    sps_table$range_frag_km2[match(sps, sps_table$species)], ";")[[1]][
      sps_table$dom_frag[match(sps, sps_table$species)]])
  
  frag_pos <- which.min(abs(a_km2 - dom_frag_km2))
  
  #store area
  range_sizes[i] <- a_km2[frag_pos]
  
  print(i)
}

#range size bins (km2)
breaks <- c(0, 1e3, 1e4, 1e5, 1e6, 1e7,
            max(range_sizes, na.rm = TRUE))

#define labels
labels <- c("<10³", "10³–10⁴", "10⁴–10⁵", "10⁵–10⁶", "10⁶–10⁷", ">10⁷")

#assign species to bins
range_bins <- cut(range_sizes, breaks = breaks, labels = labels,
                  include.lowest = TRUE)

#count species per bin
bin_counts <- table(range_bins)

#set margins
par(mar = c(5,5,3,2))

#plot the histogram
bp <- barplot(bin_counts,
              col = "grey55",
              border = "white",
              space = 0.15,
              ylab = "Number of species",
              xlab = "Range size (km²)",
              ylim = c(0, max(bin_counts) * 1.15),
              cex.names = 1.1)


#add numbers
text(bp,
     bin_counts,
     labels = bin_counts,
     pos = 3,
     cex = 1)


box(bty = "l")

#save 800


