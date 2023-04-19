
setwd("/Users/megancattau/Dropbox/1_Conferences_workshops/2023_ForestResilience_EarthLab/Recovery")
getwd()


# This code determines where the 'recovering' samples should be extracted from (i.e., how far from the perimeter edge within a fire boundary) and where the 'control' samples should be extracted from (i.e., how far from the perimeter edge outside of a fire boundary)


# Required packages:
library(rgdal)
library(sp)
library(raster)
library(tidyr)
library(dplyr)
library(assertthat)
library(sf)
library(terra)
library(lubridate)
library(stars)


### IMPORT DATA AND PROCESSING ###
# You can skip this step if you just import:
# samples<-read.csv("data/processed_data/samples.csv")


## create a data folder 

## Import EPA Level III Ecoregions 
# Available here: https://www.epa.gov/eco-research/ecoregions-north-america
# will be downloaded directly to data folder
Ecoregion_download <- file.path('data/Ecoregion', 'us_eco_l3.shp')
if (!file.exists(Ecoregion_download)) {
  from <- "https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip"
  to <- paste0('Ecoregion', ".zip")
  download.file(from, to)
  unzip(to, exdir = 'data/Ecoregion')
  unlink(to)
  assert_that(file.exists(Ecoregion_download))
}
# Import EPA Level III Ecoregion polygons, subset just Southern Rockies, and remove US-wide ecoregions
ecoregion_sf<- st_read("data/Ecoregion","us_eco_l3")
srockies_sf<-ecoregion_sf[ecoregion_sf$US_L3NAME=="Southern Rockies",] 
rm(ecoregion_sf)


## Import MTBS data
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
# Import MTBS polygons and preprocess
MTBS_fires_sf<- st_read("data/MTBS","mtbs_perims_DD")
MTBS_fires_sf$ign.date<-ymd(MTBS_fires_sf$Ig_Date)													# change to date
MTBS_fires_sf$FireYear<-year(MTBS_fires_sf$ign.date)												# new Year col
MTBS_fires_sf$FireDay<-day(MTBS_fires_sf$ign.date)													# new Day col
MTBS_fires_sf$StartMonth<-month(MTBS_fires_sf$ign.date)										# new Month col


## Import RCMap tree cover
# Available here: https://www.mrlc.gov/data?f%5B0%5D=category%3ARCMAP%20-%20Time-Series%20-%20Cover
# will be downloaded directly to Data folder
RCMAP_download <- file.path('data/RCMAP', 'rcmap_tree_2003.tif')
if (!file.exists(RCMAP_download)) {
  from <- "https://s3-us-west-2.amazonaws.com/mrlc/Tree_1997_2008.zip"
  to <- paste0('RCMAP', ".zip")
  download.file(from, to)
  unzip(to, exdir = 'data/RCMAP')
  unlink(to)
  assert_that(file.exists(RCMAP_download))
}
# Import just one year's treecover (selected 2003 bc 2002 fires most replicates in S Rockies - determined below)
tree<- raster("data/RCMAP/rcmap_tree_2003.tif")


## Align data and crop to S Rockies

# Reproject fires and S rockies boundary 
compareCRS(tree, MTBS_fires_sf)
MTBS_fires_proj<-st_transform(MTBS_fires_sf, crs(tree)) 
compareCRS(MTBS_fires_proj, tree)
rm(MTBS_fires_sf)

compareCRS(srockies_sf, tree)
srockies_proj<-st_transform(srockies_sf, crs(tree)) 
compareCRS(srockies_proj, tree)
rm(srockies_sf)
st_write(obj=srockies_proj, dsn = 'data/processed_data', layer = 'srockies_proj', driver = "ESRI Shapefile")

# Crop tree and fires to srockies
tree_sf<-st_as_stars(tree)
tree_clip<-st_crop(tree_sf, srockies_proj)
rm(tree)
rm(tree_sf)
fires_clip<-st_intersection(MTBS_fires_proj, srockies_proj)
rm(MTBS_fires_proj)
st_write(obj=fires_clip, dsn = 'data/processed_data', layer = 'fires_clip', driver = "ESRI Shapefile")


## Subset to fire year with most fires, buffer by 5km, export both layers, and crop tree data
# Get fire year with most fires just to increse replicates
data.frame(fires_clip) %>% count(as.factor(FireYear))

fire_2002<-fires_clip[fires_clip$FireYear==2002,]
st_write(obj=fire_2002, dsn = 'data/processed_data', layer = 'fire_2002', driver = "ESRI Shapefile")
fire_2002_buff<-st_buffer(fire_2002, dist=5000, dissolve=FALSE)
st_write(obj=fire_2002_buff, dsn = 'data/processed_data', layer = 'fire_2002_buff', driver = "ESRI Shapefile")

# Crop tree data to buffered 2002 fire, create sample point locations and export
tree_clip_2002<-st_crop(tree_clip, fire_2002_buff)
tree_clip_2002_ras<- rast(tree_clip_2002) # takes forever, and values get all messed up here, so don't use these for tree percent
tree_clip_2002_points<-as.points(tree_clip_2002_ras, values=TRUE, na.rm=TRUE)
writeVector(x= tree_clip_2002_points, filename="data/processed_data/tree_clip_2002_points", filetype="ESRI Shapefile")
# If need to import back in: tree_clip_2002_points<-vect("data/processed_data/tree_clip_2002_points")

rm(tree_clip)
rm(tree_clip_2001)
rm(tree_clip_2001_ras)


## Create distance to fire layers in QGIS, and sample
# In QGIS:
# create distance to fire perimeter layer (going out)
# Add FID to fire 2002 attribute table, rasterize (10 georef units, extent of buffered fire and -9999 as ND value), proximity (max dist 6000, 1 as target value)
fire_2002_dist<-terra::rast("data/processed_data/fire_2002_dist.tif")

# create distance to fire perimeter layer (going in) 
# fires 2002 to lines, rasterize (same param as above), proximity (max dist 20000 bc 6000 didnt go all the way into interior of fires, 1 as target value, fire_2002_dist_in), clip to fires
fire_2002_dist_in_clip<-terra::rast("data/processed_data/fire_2002_dist_in_clip.tif")

# Sample raster layers at tree points within buffered fire: distance to fire, inside or outside fire perimeter, tree value 
# extract values at sample points
samples1<-terra::extract(x= fire_2002_dist, y= tree_clip_2002_points, xy=TRUE)
samples2<-terra::extract(x= fire_2002_dist_in_clip, y= tree_clip_2002_points)
rm(fire_2002_dist)
rm(fire_2002_dist_in_clip)
tree2<- rast("data/RCMAP/rcmap_tree_2003.tif")
samples3<-terra::extract(x= tree2, y= tree_clip_2002_points)
samples<-cbind(samples1, samples2[,2], samples3[,2])
names(samples)<-c("ID", "dist", "x", "y", "dist_in", "tree")
write.csv(samples, "data/processed_data/samples.csv")



### MODEL ###

# We use two approaches to determine where the 'recovering' samples should be extracted from (i.e., how far from the perimeter edge within a fire boundary) and where the 'control' samples should be extracted from (i.e., how far from the perimeter edge outside of a fire boundary)
# 1) t tests and ks tests to determine at what distance outside of / into the fire perimeter tree cover percent values stop differing in their mean and distribution, respectively. That threshold would be the point beyond which we no longer see edge effects. We test this at 50m increments
# 2) Segmentation analysis to examine tree cover percent value as a function of distance to fire perimeter. The point at which that slope changes determines at what distance outside of / into the fire perimeter would be the point beyond which we no longer see edge effects.

# If need to import back in:
samples<-read.csv("data/processed_data/samples.csv")
samples<-samples[,-1]

# OUTSIDE OF PERIMETER - "unburned' pixels
# limit to 3km
samples_df_out<-data.frame(samples$tree, samples$dist)
samples_df_out<-samples_df_out[samples_df_out$samples.dist<3000,]
names(samples_df_out)<-c("tree", "dist")
dist_vec<-seq(50, 2950, 50)

for (i in 1:length(dist_vec)){
	samples_df_out[,i+2]<-ifelse(samples_df_out$dist<dist_vec[i], "close", "far")
}

samples_df_out_names<-character(length(dist_vec))
for (i in 1:length(dist_vec)){
	samples_df_out_names[i]<-paste0('m', dist_vec[i])
}	
names(samples_df_out)<-c("tree", "dist", samples_df_out_names)


## t test
t_results<-vector("list", length(dist_vec)-1)		# empty list
t_p<-vector("numeric", length(dist_vec)-1)		# empty vector

for (i in 1:(length(dist_vec)-1)){
	t_results[[i]]<-t.test(samples_df_out[samples_df_out[,i+2]=="close", 1], samples_df_out[samples_df_out[,i+2]=="far", 1])
}
t_results

# Just look at p values
for (i in 1:length(t_results)){
	t_p[i]<-t_results[[i]]$p.value
}
t_p
# means all significantly different from one another (e.g., 50m from edge diff tree values from those farther away), with difference between groups highest at 50m threshold and decreasing as threshold increases


## ks test
ks_results<-vector("list", length(dist_vec)-1)		# empty list
k_p<-vector("list", length(dist_vec)-1)		# empty list

for (i in 1:(length(dist_vec)-1)){
	ks_results[[i]]<-ks.test(samples_df_out[samples_df_out[,i+2]=="close", 1], samples_df_out[samples_df_out[,i+2]=="far", 1])
}
ks_results

for (i in 1:length(ks_results)){
	k_p[[i]]<-ks_results[[i]][2]
}
k_p
# distributions all significantly different from one another (e.g., 50m from edge diff tree values from those farther away), with d value (magnitude of difference between groups) highest at 50m threshold and decreasing as threshold increases


# INSIDE PERIMETER - "Recovering" pixels
samples_df_in<-data.frame(samples$tree, samples$dist_in)
samples_df_in<-samples_df_in[samples_df_in$samples.dist_in<3000,]

for (i in 1:length(dist_vec)){
	samples_df_in[,i+2]<-ifelse(samples_df_in$samples.dist_in<dist_vec[i], "close", "far")
}

names(samples_df_in)<-c("tree", "dist_in", samples_df_out_names)


## t test
t_results_in<-vector("list", length(dist_vec)-1)		# empty list
t_p_in<-vector("numeric", length(dist_vec)-1)		# empty list
t_diff<-vector("numeric", length(dist_vec)-1)		# empty vector

for (i in 1:(length(dist_vec)-1)){
	t_results_in[[i]]<-t.test(samples_df_in[samples_df_in[,i+2]=="close", 1], samples_df_in[samples_df_in[,i+2]=="far", 1])
}

t_results_in

for (i in 1:(length(dist_vec)-1)){
	t_p_in[i]<-t_results_in[[i]]$p.value
}
t_p_in

for (i in 1:(length(dist_vec)-1)){
	t_diff[i]<-t_results_in[[i]]$estimate[1]-t_results_in[[i]]$estimate[2]
}
t_diff

# means all significantly different from one another (e.g., 50m from edge diff tree values from those farther away), with difference between groups highest at 100m threshold


## ks test
ks_results_in<-vector("list", length(dist_vec)-1)		# empty list
k_p_in<-vector("list", length(dist_vec)-1)		# empty list

for (i in 1:(length(dist_vec)-1)){
	ks_results_in[[i]]<-ks.test(samples_df_in[samples_df_in[,i+2]=="close", 1], samples_df_in[samples_df_in[,i+2]=="far", 1])
}

ks_results_in

for (i in 1:length(ks_results_in)){
	k_p_in[[i]]<-ks_results_in[[i]][2]
}
k_p_in

# distributions all significantly different from one another (e.g., 50m from edge diff tree values from those farther away) and d value highest 100m in from edge

	
## linear model and segmented analysis	
reg_out<-lm(tree~dist, data=samples_df_out)

# segmented for breakpoint analysis
library(segmented)
davies_out<-davies.test(obj=reg_out, seg.Z=~dist, k = 10)	# Davies' test on regression model object. If significant, there's a breakpoint
davies_out						# Davies' test indicates that there is a breakpoint w the p value ("best" at 222)
seg_out<-segmented(obj=reg_out, seg.Z=~dist, psi=200)	# estimate the best breakpoint
# Estimated breakpoint at 16 m

	
reg_in<-lm(tree~dist_in, data=samples_df_in)
davies_in<-davies.test(obj=reg_in, seg.Z=~dist_in, k = 10)
davies_in
seg_in<-segmented(obj=reg_in, seg.Z=~dist_in, psi=100)	
# Estimated breakpoint at 66 m



### All evidence points to threshold 50m out from fire perimeter and 100m into fire perimeter, which we use as liberal thresholds below ###


### CREATE FIRE BUFFERS ###

## Conservative threshold of 375m in both directions:

# "Unburned" shapefile
fires_clip_unburned_375m<-st_buffer(fires_clip, 375)
st_write(obj=fires_clip_unburned_375m, dsn = 'data/processed_data', layer = 'fires_clip_unburned_375m', driver = "ESRI Shapefile")
fires_clip_unburned_5km<-st_buffer(fires_clip, 5000)
st_write(obj=fires_clip_unburned_5km, dsn = 'data/processed_data', layer = 'fires_clip_unburned_5km', driver = "ESRI Shapefile")
# Difference tool in QGIS to get unburned samples (called unburned.shp), meaning all areas 375m from fire and within 5km

# "Recovery" shapefile
fires_clip_recovery_375<-st_buffer(fires_clip, -375)
st_write(obj=fires_clip_recovery_375, dsn = 'data/processed_data', layer = 'recovery', driver = "ESRI Shapefile")


## Liberal threshold of 50m outside of fire and 100m into:

# "Unburned" shapefile
fires_clip_unburned_50<-st_buffer(fires_clip, 50)
st_write(obj=fires_clip_unburned_50, dsn = 'data/processed_data', layer = 'fires_clip_unburned_50', driver = "ESRI Shapefile")
# Difference tool in QGIS to get unburned samples (called unburned50.shp), meaning all areas 50m from fire and within 5km

# "Recovery" shapefile

fires_clip_recovery_100<-st_buffer(fires_clip, -100)
st_write(obj=fires_clip_recovery_100, dsn = 'data/processed_data', layer = 'recovery_100', driver = "ESRI Shapefile")


