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

# project worldclim, crop it, etc...use as snap raster below
snap_raster = projectRaster(worldclim.stack[[1]], res=1000, method = "bilinear", crs = sr_cea) # problem is I save worldclim as a raster with bands
snap_raster = mask(snap_raster, water_bodies_proj, inverse=TRUE)

# fires
coordinates(fires) = c("longitude","latitude")
fires$fire = 1
fire.density = rasterize(fires,worldclim.stack[[1]],field='fire',fun='count')
names(fire.density) = "fire_density"
fire.density[is.na(fire.density)] = 0
fire.density = mask(fire.density,worldclim.stack[[1]])
fire.density = mask(fire.density, water_bodies, inverse=TRUE)
  
### aggregating 

mean_elevation = aggregate(elevation,11,mean,na.rm=TRUE)
names(mean_elevation) = "mean_elevation"
max_elevation = aggregate(elevation,11,max,na.rm=TRUE)
min_elevation = aggregate(elevation,11,min,na.rm=TRUE)
range_elevation = max_elevation - min_elevation # this is like roughness?
std_elevation = aggregate(elevation,11,sd,na.rm=TRUE) # this is like roughness?
names(std_elevation) = "std_elevation"
max_slope = aggregate(slope,11,max,na.rm=TRUE)
names(max_slope) = "max_slope"
mean_slope = aggregate(slope,11,mean,na.rm=TRUE)
names(mean_slope) = "mean_slope"

# globcover grain exactly 1/3 of 1 km, use the mode
globcover_forest_mask = mask(globcover_forest,studyArea)
globcover_forest_agg = aggregate(globcover_forest_mask,3,modal,na.rm=TRUE)
names(globcover_forest_agg) = "forest"

# EVI grain 1/4 of 1 km
evi_crop = crop(evi,extent(studyArea))
evi_mask = mask(evi_crop, studyArea)
evi_agg = aggregate(evi_mask,4,mean,na.rm=TRUE)

# VCF grain 1/4 of 1 km
vcf_crop = crop(vcf,extent(studyArea))
vcf_mask = mask(vcf_crop, studyArea)
vcf_agg = aggregate(vcf_mask,4,mean,na.rm=TRUE)
names(vcf_agg) = "treecover2010"

# Treecover grain is about 1/36 of 1 km
treecover2000_agg = aggregate(treecover2000,36,mean,na.rm=TRUE)
names(treecover2000_agg) = "treecover2000"
treecover2000_mask = mask(treecover2000_agg,studyArea)

forestloss_agg = aggregate(forestloss,36,mean,na.rm=TRUE) # changed sum to mean
names(forestloss_agg) = "forestloss"
forestloss_mask = mask(forestloss_agg,studyArea)

datamask_agg = aggregate(datamask,36,sum,na.rm=TRUE)
names(datamask_agg) = "datamask"
datamask_mask = mask(datamask_agg,studyArea)

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

continuous_layers = stack(projected_rasters)

rasters_to_project = list(pa_raster, globcover_forest_agg)
projected_rasters = lapply(rasters_to_project, FUN = projectRaster, to = snap_raster, method = "ngb", crs = sr_cea)

binary_layers = stack(projected_rasters)

fw = focalWeight(continuous_layers$forestloss, 10000, type='circle') 
fw[fw>0] = 1 
forestloss_focal = focal(continuous_layers$forestloss, w=fw, fun=mean, na.rm=TRUE)
forestloss_focal = mask(forestloss_focal,water_bodies_proj,inverse=TRUE)
names(forestloss_focal) = "forestloss_focal"

fw = focalWeight(continuous_layers$road_length, 10000, type='circle') 
fw[fw>0] = 1 
rLength_focal = focal(continuous_layers$road_length, w=fw, fun=mean, na.rm=TRUE)
rLength_focal = mask(rLength_focal,water_bodies_proj,inverse=TRUE)
names(rLength_focal) = "road_length_focal"

fw = focalWeight(continuous_layers$fire_density, 10000, type='circle') 
fw[fw>0] = 1 
fire_density_focal = focal(continuous_layers$fire_density, w=fw, fun=mean, na.rm=TRUE)
fire_density_focal = mask(fire_density_focal,water_bodies_proj,inverse=TRUE)
names(fire_density_focal) = "fire_density_focal"

continuous_layers_final = stack(continuous_layers,forestloss_focal, rLength_focal, fire_density_focal)
continuous_layers_final = mask(continuous_layers_final, calc(continuous_layers_final,fun = sum)) # mask through NAs just to make sure
plot(continuous_layers_final$max_slope) # test to see if we have NAs in the right place, i.e. mostly the ocean, before scaling

# calculate quadratic terms, scale
continuous_layers_scale = scale(continuous_layers_final)
continuous_layers_sq = continuous_layers_scale^2
names(continuous_layers_sq) <- paste(names(continuous_layers_final),"_sq",sep="")

# do distancePA separately since we want to preserve the meaning of 0
distancePA_proj = projectRaster(distancePA, to = snap_raster, method = "bilinear", crs = sr_cea)
distancePA_mask = mask(distancePA_proj,water_bodies_proj,inverse=TRUE)
distancePA_scale = scale(distancePA_mask, center=FALSE)

# all layers
all_layers = stack(continuous_layers_scale,distancePA_scale,continuous_layers_sq,binary_layers)

# write 1 km layers
outputFolder = paste0(dataPrep_Folder,"1km/")
dir.create(file.path(outputFolder), showWarnings = FALSE)
unstack_layers = unstack(all_layers)
for(i in seq_along(unstack_layers)){
  writeRaster(unstack_layers[[i]], file=paste0(outputFolder,names(unstack_layers[[i]]), ".tif"), overwrite=TRUE)}

########################################
### Aggregate to coarser resolutions ###
########################################

library(raster)
dataPrep_Folder = "/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_output/"
outputFolder = paste0(dataPrep_Folder,"1km/")
rasterFiles = list.files(outputFolder, full.names = TRUE)
# remove quadratic terms, so we can do them again, after aggregating
rasterFiles = rasterFiles[!grepl("_sq",rasterFiles)]
all_layers = stack(rasterFiles)
distancePA_scale = subset(all_layers,"distancePA")
binary_layer_names = c("protected","forest")
binary_layers = subset(all_layers,binary_layer_names)
continuous_layer_names = names(all_layers)[!(names(all_layers) %in% c("protected","forest","distancePA"))]
continuous_layers_final = subset(all_layers,continuous_layer_names)

resolutions = c(2,4,8,16,32,64,128)

for(resolution in resolutions){
  
  continuous_layers_agg = aggregate(continuous_layers_final, resolution, fun=mean, na.rm=TRUE)
  continuous_layers_scale = scale(continuous_layers_agg)
  continuous_layers_sq = continuous_layers_scale^2
  names(continuous_layers_sq) <- paste(names(continuous_layers_agg),"_sq",sep="")
  #continuous_layers_sq = scale(continuous_layers_sq)
  
  # binary layers
  binary_layers_agg  = aggregate(binary_layers, resolution, fun=modal, na.rm=TRUE)
  
  distancePA_scale_agg = aggregate(distancePA_scale, resolution, fun=mean, na.rm=TRUE)
  distancePA_scale_agg = scale(distancePA_scale_agg, center=FALSE)
  
  # all layers
  all_layers_agg = stack(continuous_layers_scale,distancePA_scale_agg,continuous_layers_sq,binary_layers_agg)
  
  # write layers
  outputFolder = paste0(dataPrep_Folder,resolution,"km/")
  dir.create(file.path(outputFolder), showWarnings = FALSE)
  unstack_layers = unstack(all_layers_agg)
  for(i in seq_along(unstack_layers)){
    writeRaster(unstack_layers[[i]], file=paste0(outputFolder,names(unstack_layers[[i]]), ".tif"), overwrite=TRUE)}
  
}
