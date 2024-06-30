#load packages
library(sf)

## Shapefile contaioning ranges of all terrestrial mammals
#source: https://www.iucnredlist.org/resources/spatial-data-download

wd_ranges <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/All ranges"
wd_out <- "/Users/carloseduardoaribeiro/Documents/Post-doc/SHAP/Mammals/Range_maps"

ranges <- st_read(dsn = wd_ranges, layer = "MAMMALS_TERRESTRIAL_ONLY")

sps_list <- unique(ranges$sci_name)

for(i in 1:length(sps_list))
{
  sps_range <- ranges[which(ranges$sci_name == sps_list[i]),]
  st_write(sps_range,layer = sps_list[i], driver = "ESRI Shapefile",
           dsn = wd_out)
  print(i)
}

ranges_mar_terrestrial <- st_read(dsn = wd_ranges,
                                  layer = "MAMMALS_MARINE_AND_TERRESTRIAL")

sps_list <- unique(ranges_mar_terrestrial$sci_name)

for(i in 1:length(sps_list))
{
  sps_range <- ranges_mar_terrestrial[which(
    ranges_mar_terrestrial$sci_name == sps_list[i]),]
  st_write(sps_range,layer = sps_list[i], driver = "ESRI Shapefile",
           dsn = wd_out)
  print(i)
}
