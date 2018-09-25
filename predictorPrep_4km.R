library(dismo)
library(rgdal)
library(spatstat)
library(maptools)

rasterOptions(tmpdir="/Volumes/RasterTemp/")
downloadFolder = "/Volumes/RasterTemp/"

# projection to use
sr_cea = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sr_geog = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# establish folders
inputFolder = "/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_input/"
dataPrep_Folder = "/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_output/"

# load raw data
studyArea = readOGR(paste0(inputFolder,"tapir_studyArea_new.shp"),"tapir_studyArea_new")
studyArea_proj = spTransform(studyArea,sr_cea)
water_bodies_proj = readOGR(paste0(inputFolder,"water_bodies_proj.shp"),"water_bodies_proj")
water_bodies = spTransform(water_bodies_proj, sr_geog)
evi = raster(paste0(inputFolder,"EVI.tif"))/10000
vcf = raster(paste0(inputFolder,"VCF.tif"))
pa_polys = readOGR(dsn=paste0(inputFolder,"WDPA_selection.shp"),layer="WDPA_selection")
treecover2000 = raster(paste0(downloadFolder,"treecover2000.tif"))
forestloss = raster(paste0(downloadFolder,"forestloss.tif"))
datamask = raster(paste0(downloadFolder,"datamask.tif"))
elevation = raster(paste0(downloadFolder,"elevation.tif"))
slope = raster(paste0(downloadFolder,"slope.tif"))
globcover_forest = raster(paste0(downloadFolder,"globcover_forest.tif"))
fires = read.csv(paste0(downloadFolder,"DL_FIRE_M6_22345/fire_archive_M6_22345.csv"))
rLength = raster(paste0(inputFolder,"rLength.tif"))
names(rLength) = "road_length"
roads_dist = raster(paste0(downloadFolder,"road_distance.tif"))
water_occurrence = raster(paste0(downloadFolder,"water_occurrence.tif"))

chelsa.files = list.files(paste0(downloadFolder),full.names=T,pattern="CHELSA_")
chelsa.list = lapply(chelsa.files, raster)
chelsa.stack = stack(chelsa.list)
# http://chelsa-climate.org/bioclim/
chelsa.names = c("annual_mean_temp","mean_temp_warmest_quarter","precip_wettest_month","precip_driest_month",
                 "mean_temp_coldest_quarter","annual_precip","precip_seasonality","precip_wettest_quarter",
                 "precip_driest_quarter","precip_warmest_quarter","precip_coldest_quarter","mean_diurnal_range",
                 "isothermality","temp_seasonality","max_temp_warmest_month","min_temp_coldest_month",
                 "temp_annual_range","mean_temp_wettest_quarter","mean_temp_driest_quarter")
names(chelsa.stack) = paste0("chelsa_",chelsa.names)

worldclim.files = list.files(paste0(downloadFolder),full.names=T,pattern="worldclim_")
worldclim.list = lapply(worldclim.files, raster)
worldclim.stack = stack(worldclim.list)

### aggregating 

mean_elevation = aggregate(elevation,(11*4),mean,na.rm=TRUE)
names(mean_elevation) = "mean_elevation"
max_elevation = aggregate(elevation,(11*4),max,na.rm=TRUE)
min_elevation = aggregate(elevation,(11*4),min,na.rm=TRUE)
range_elevation = max_elevation - min_elevation # this is like roughness?
std_elevation = aggregate(elevation,(11*4),sd,na.rm=TRUE) # this is like roughness?
names(std_elevation) = "std_elevation"
max_slope = aggregate(slope,(11*4),max,na.rm=TRUE)
names(max_slope) = "max_slope"
mean_slope = aggregate(slope,(11*4),mean,na.rm=TRUE)
names(mean_slope) = "mean_slope"

# globcover grain exactly 1/3 of 1 km, use the mode
globcover_forest_mask = mask(globcover_forest,studyArea)
globcover_forest_agg = aggregate(globcover_forest_mask,(3*4),modal,na.rm=TRUE)
names(globcover_forest_agg) = "forest"

# EVI grain 1/4 of 1 km
evi_crop = crop(evi,extent(studyArea))
evi_mask = mask(evi_crop, studyArea)
evi_agg = aggregate(evi_mask,(3*4),mean,na.rm=TRUE)

# VCF grain 1/4 of 1 km
vcf_crop = crop(vcf,extent(studyArea))
vcf_mask = mask(vcf_crop, studyArea)
vcf_agg = aggregate(vcf_mask,(4*4),mean,na.rm=TRUE)
names(vcf_agg) = "treecover2010"

# water occurrence grain is about
water_occurrence_agg = aggregate(water_occurrence,(36*4),sum,na.rm=TRUE)
names(water_occurrence_agg) = "water_occurence"
water_occurrence_mask = mask(water_occurrence_agg,studyArea)

# Treecover grain is about 1/36 of 1 km
treecover2000_agg = aggregate(treecover2000,(36*4),mean,na.rm=TRUE)
names(treecover2000_agg) = "treecover2000"
treecover2000_mask = mask(treecover2000_agg,studyArea)

forestloss_agg = aggregate(forestloss,(36*4),mean,na.rm=TRUE) # changed sum to mean
names(forestloss_agg) = "forestloss"
forestloss_mask = mask(forestloss_agg,studyArea)

datamask_agg = aggregate(datamask,(36*4),sum,na.rm=TRUE)
names(datamask_agg) = "datamask"
datamask_mask = mask(datamask_agg,studyArea)

########################################
### all below depends on snap raster ###
########################################

# project worldclim, crop it, etc...use as snap raster below
snap_raster = projectRaster(worldclim.stack[[1]], res=4000, method = "bilinear", crs = sr_cea) # problem is I save worldclim as a raster with bands
snap_raster = mask(snap_raster, water_bodies_proj, inverse=TRUE)

# load here if coming in manually
#outputFolder = paste0(dataPrep_Folder,"4km_chp4/")
#snap_raster = raster(outputFolder,"mean_slope.tif")

# fires
coordinates(fires) = c("longitude","latitude")
fires$fire = 1
fire.density = rasterize(fires,worldclim.stack[[1]],field='fire',fun='count')
names(fire.density) = "fire_density"
fire.density[is.na(fire.density)] = 0
fire.density = mask(fire.density,worldclim.stack[[1]])
fire.density = mask(fire.density, water_bodies, inverse=TRUE)

# PA stuff
pa_polys = crop(pa_polys,studyArea)
pa_polys_proj = spTransform(pa_polys,sr_cea)
pa_polys_proj$protected = 1
pa_raster = raster(snap_raster)
pa_raster = rasterize(pa_polys_proj,pa_raster,field="protected")

# calculate distance to PA
pa_dist1 = distance(pa_raster)
pa_raster_inverse = pa_raster
pa_raster_inverse[is.na(pa_raster_inverse)] = 2
pa_raster_inverse[pa_raster_inverse == 1] = NA
pa_dist2 = distance(pa_raster_inverse)
pa_dist2 = pa_dist2*-1
distancePA = mask(pa_dist1 + pa_dist2, snap_raster)
names(distancePA) = "distancePA"

pa_raster[is.na(pa_raster)] = 0
pa_raster = mask(pa_raster, snap_raster)
names(pa_raster) = "protected"

# project and resample at same time
rasters_to_project = list(worldclim.stack, chelsa.stack, mean_elevation, mean_slope, std_elevation, max_slope,
                          evi_agg, vcf_agg, treecover2000_mask, forestloss_mask, 
                          roads_dist, rLength, fire.density)
projected_rasters = lapply(rasters_to_project, FUN = projectRaster, to = snap_raster, method = "bilinear", crs = sr_cea)

#projected_raster = projectRaster(water_occurrence_mask,to = snap_raster, method = "bilinear", crs = sr_cea)
#projected_raster_temp = projected_raster
#projected_raster_temp[projected_raster_temp==0]=NA
#water_distance = distance(projected_raster_temp)
#water_distance = mask(water_distance,projected_raster)
#names(water_distance) = "water_distance"
#continuous_layers_final = stack(water_distance,projected_raster)

continuous_layers = stack(projected_rasters)

rasters_to_project = list(pa_raster, globcover_forest_agg)
projected_rasters = lapply(rasters_to_project, FUN = projectRaster, to = snap_raster, method = "ngb", crs = sr_cea)

binary_layers = stack(projected_rasters)

continuous_layers_final = mask(continuous_layers, calc(continuous_layers,fun = sum)) # mask through NAs just to make sure
plot(continuous_layers_final$max_slope) # test to see if we have NAs in the right place, i.e. mostly the ocean, before scaling

# calculate quadratic terms, scale
continuous_layers_scale = scale(continuous_layers_final)
continuous_layers_sq = continuous_layers_scale^2
names(continuous_layers_sq) <- paste(names(continuous_layers_final),"_sq",sep="")

# do distancePA separately since we want to preserve the meaning of 0
distancePA_proj = projectRaster(distancePA, to = snap_raster, method = "bilinear", crs = sr_cea)
distancePA_mask = mask(distancePA_proj,water_bodies_proj,inverse=TRUE)
distancePA_scale = scale(distancePA_mask, center=FALSE)

### create interaction variables here, and add to the all_layers
treecover2000_EVI = continuous_layers_scale$treecover2000 * continuous_layers_scale$EVI
names(treecover2000_EVI) = "treecover2000_EVI"

# all layers
all_layers = stack(continuous_layers_scale,distancePA_scale,continuous_layers_sq,binary_layers,treecover2000_EVI)
#all_layers = stack(continuous_layers_scale,continuous_layers_sq)

# write layers
outputFolder = paste0(dataPrep_Folder,"4km_chp4/")
dir.create(file.path(outputFolder), showWarnings = FALSE)
unstack_layers = unstack(all_layers)
for(i in seq_along(unstack_layers)){
  writeRaster(unstack_layers[[i]], file=paste0(outputFolder,names(unstack_layers[[i]]), ".tif"), overwrite=TRUE)}

# to add interactions

distancePA = raster(paste0(outputFolder,"distancePA.tif"))
protected = raster(paste0(outputFolder,"protected.tif"))
treecover2000 = raster(paste0(outputFolder,"treecover2000.tif"))
distancePA_treecover = distancePA * treecover2000
names(distancePA_treecover) = "distancePA_treecover"
writeRaster(distancePA_treecover, file=paste0(outputFolder,"distancePA_treecover.tif"), overwrite=TRUE)
protected_treecover = protected * treecover2000
names(protected_treecover) = "protected_treecover"
writeRaster(protected_treecover, file=paste0(outputFolder,"protected_treecover.tif"), overwrite=TRUE)

chelsa_precip_driest_quarter = raster(paste0(outputFolder,"chelsa_precip_driest_quarter.tif")) 
chelsa_mean_temp_driest_quarter = raster(paste0(outputFolder,"chelsa_mean_temp_driest_quarter.tif")) 
precip_temp = chelsa_precip_driest_quarter*chelsa_mean_temp_driest_quarter
names(precip_temp) = "precip_temp"
writeRaster(precip_temp, file=paste0(outputFolder,"precip_temp.tif"), overwrite=TRUE)

twi = raster(paste0(outputFolder,"twi.tif")) 
precip_temp_twi = precip_temp*twi
names(precip_temp_twi) = "precip_temp_twi"
writeRaster(precip_temp_twi, file=paste0(outputFolder,"precip_temp_twi.tif"), overwrite=TRUE)
