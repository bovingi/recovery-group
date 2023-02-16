
setwd("/Users/megancattau/Dropbox/1_Conferences_workshops/2023_ForestResilience_EarthLab/Recovery")
getwd()

# Make shapefile of inside/outside fire
# 100 m buffer bt fire and 
# further 500m?
# matched by elevation, aspect/slope


# Import MTBS data
# Available here: https://www.epa.gov/eco-research/ecoregions-north-america
# will be downloaded directly to Data folder
MTBS_download <- file.path('data/MTBS', 'mtbs_perims_DD.shp')
if (!file.exists(MTBS_download)) {
  from <- "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/burned_area_extent_shapefile/mtbs_perimeter_data.zip"
  to <- paste0('MTBS', ".zip")
  download.file(from, to)
  unzip(to, exdir = 'data/MTBS')
  unlink(to)
  assert_that(file.exists(MTBS_download))
}
# Import MTBS polygons
MTBS_fires<- readOGR("data/Ecoregion","us_eco_l3")
crs(ecoregion)
class(ecoregion)
srockies<-ecoregion[ecoregion@data$US_L3NAME=="Southern Rockies",] 


# Look at VCF as function of distance from forest edge to inform when stops increasing and if levels off which will inform:
# make buffer inside fire
# make buffer outside fire
# make further 500m? buffer of comparable size to fire area

# Import DEM
# Calculate mean slope, aspect, elevation 
# make sure not too different from fire
