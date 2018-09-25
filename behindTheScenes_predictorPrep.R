library(dismo)
library(rgdal)
library(spatstat)
library(maptools)
library(osmar) # osmosis must be installed, used homebrew

rasterOptions(tmpdir="/Volumes/RasterTemp/")
downloadFolder = "/Volumes/RasterTemp/"
inputFolder = "/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_input/"

# load study area
studyArea = readOGR(paste0(inputFolder,"tapir_studyArea_new.shp"),"tapir_studyArea_new")

# load water bodies data
glwd_1 = readOGR(paste0(inputFolder,"glwd_1.shp"),"glwd_1")
glwd_2 = readOGR(paste0(inputFolder,"glwd_2.shp"),"glwd_2")
glwd_1_crop = crop(glwd_1,studyArea)
glwd_2_crop = crop(glwd_2,studyArea)
proj4string(glwd_1_crop) = sr_geog
proj4string(glwd_2_crop) = sr_geog
glwd_1_proj = spTransform(glwd_1_crop,sr_cea)
glwd_1_proj = as(glwd_1_proj,"SpatialPolygons")
glwd_2_proj = spTransform(glwd_2_crop,sr_cea)
glwd_2_proj = as(glwd_2_proj,"SpatialPolygons")
water_bodies_proj = rbind(glwd_1_proj,glwd_2_proj)
water_bodies_proj = as(water_bodies_proj,"SpatialPolygonsDataFrame")
writeOGR(water_bodies_proj,dsn=paste0(inputFolder,"water_bodies_proj.shp"), layer="water_bodies_proj", driver="ESRI Shapefile")

# roads
# http://download.geofabrik.de/south-america/colombia-latest-free.shp.zip
colombia_roads = readOGR(paste0(inputFolder,"colombia_roads.shp"),"colombia_roads")
colombia_roads_crop = crop(colombia_roads,studyArea)
colombia_roads_proj = spTransform(colombia_roads_crop,sr_cea)
writeOGR(colombia_roads_proj,dsn=paste0(inputFolder,"colombia_roads_proj.shp"), layer="colombia_roads_proj", driver="ESRI Shapefile")
colombia_roads_proj = readOGR(paste0(inputFolder,"colombia_roads_proj.shp"),"colombia_roads_proj")

# http://download.geofabrik.de/north-america/mexico-latest-free.shp.zip
mexico_roads = readOGR(paste0(inputFolder,"mexico_roads.shp"),"mexico_roads")
mexico_roads_crop = crop(mexico_roads,studyArea)
mexico_roads_proj = spTransform(mexico_roads_crop,sr_cea)
writeOGR(mexico_roads_proj,dsn=paste0(inputFolder,"mexico_roads_proj.shp"), layer="mexico_roads_proj", driver="ESRI Shapefile")
mexico_roads_proj = readOGR(paste0(inputFolder,"mexico_roads_proj.shp"),"mexico_roads_proj")

ca_roads = readOGR(paste0(inputFolder,"lines_hwy_clip.shp"),"lines_hwy_clip")
ca_roads_proj = spTransform(ca_roads,sr_cea)
writeOGR(ca_roads_proj,dsn=paste0(inputFolder,"central_america_roads_proj.shp"), 
         layer="central_america_roads_proj", driver="ESRI Shapefile")
ca_roads_proj = readOGR(paste0(inputFolder,"central_america_roads_proj.shp"),"central_america_roads_proj")

colombia_roads_proj = as(colombia_roads_proj,"SpatialLines")
mexico_roads_proj = as(mexico_roads_proj,"SpatialLines")
ca_roads_proj = as(ca_roads_proj,"SpatialLines")
roads_merge = rbind(colombia_roads_proj,mexico_roads_proj,ca_roads_proj, makeUniqueIDs = TRUE)

roads_merge$id = row.names(roads_merge)
writeOGR(roads_merge,dsn=paste0(inputFolder,"roads_merge.shp"), 
         layer="roads_merge", driver="ESRI Shapefile")

roads_merge = readOGR(paste0(inputFolder,"roads_merge.shp"),"roads_merge")

# Hansen data
tiles = c("30N_100W","30N_090W","20N_100W","20N_090W","20N_080W","10N_090W","10N_080W")
for(tile in tiles){
  download.file(paste0("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2016-v1.4/Hansen_GFC-2016-v1.4_treecover2000_",tile,".tif"),
                paste0(downloadFolder,"Hansen_GFC2016_treecover2000_",tile,".tif"),mode="wb")
  treecover2000 = raster(paste0(downloadFolder,"Hansen_GFC2016_treecover2000_",tile,".tif"))
  download.file(paste0("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2016-v1.4/Hansen_GFC-2016-v1.4_lossyear_",tile,".tif"),
                paste0(downloadFolder,"Hansen_GFC2016_lossyear_",tile,".tif"),mode="wb")
  forestloss = raster(paste0(downloadFolder,"Hansen_GFC2016_lossyear_",tile,".tif"))
  download.file(paste0("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2016-v1.4/Hansen_GFC-2016-v1.4_datamask_",tile,".tif"),
                paste0(downloadFolder,"Hansen_GFC2016_datamask_",tile,".tif"),mode="wb")
  datamask = raster(paste0(downloadFolder,"Hansen_GFC2016_datamask_",tile,".tif"))
  
  datamask[datamask==2] = NA
  forestloss[forestloss>0] = 1
  
  datamask = crop(datamask,file=paste0(downloadFolder,"Hansen_GFC2016_datamask_crop_",tile,".tif"),extent(studyArea))
  forestloss = crop(forestloss,file=paste0(downloadFolder,"Hansen_GFC2016_forestloss_crop_",tile,".tif"),extent(studyArea))
  treecover2000 = crop(treecover2000,file=paste0(downloadFolder,"Hansen_GFC2016_treecover2000_crop_",tile,".tif"),extent(studyArea))
  
  datamask = mask(datamask,file=paste0(downloadFolder,"Hansen_GFC2016_datamask_mask_",tile,".tif"),studyArea)
  
  forestloss = mask(forestloss,file=paste0(downloadFolder,"Hansen_GFC2016_forestloss_mask_",tile,".tif"),datamask)
  treecover2000 = mask(treecover2000,file=paste0(downloadFolder,"Hansen_GFC2016_treecover2000_mask_",tile,".tif"),datamask)
  
}


treecover_files = list.files(downloadFolder,pattern="treecover",full.names=TRUE)
treecover_list = lapply(treecover_files, raster)
rasters.mosaicargs = treecover_list
rasters.mosaicargs$fun = mean
rasters.mosaicargs$filename = paste0(downloadFolder,"treecover2000.tif")
treecover2000 = do.call(mosaic, rasters.mosaicargs)

forestloss_files = list.files(downloadFolder,pattern="loss",full.names=TRUE)
forestloss_list = lapply(forestloss_files, raster)
rasters.mosaicargs = forestloss_list
rasters.mosaicargs$fun = max
rasters.mosaicargs$filename = paste0(downloadFolder,"forestloss.tif")
forestloss = do.call(mosaic, rasters.mosaicargs)

datamask_files = list.files(downloadFolder,pattern="datamask",full.names=TRUE)
datamask_list = lapply(datamask_files, raster)
rasters.mosaicargs = datamask_list
rasters.mosaicargs$fun = max
rasters.mosaicargs$filename = paste0(downloadFolder,"datamask.tif")
datamask = do.call(mosaic, rasters.mosaicargs)

#download elevation data for all countries in study area, mosaic, and crop/mask
### I actually had to download these manually because the getData did not work in this instance
srtm_list = list()
i = 1
tiles = list(c(76,1), c(-76,6),c(-81,6),c(-86,11),c(-81,11),c(-91,16),c(-86,16),c(-86,21),c(-91,11),c(-76,1),c(-96,16),c(-81,16),c(-86,6),c(-91,21),c(-101,16))
for(tile in tiles){
  srtm_list[[i]] = raster::getData('SRTM', lon=tile[1], lat=tile[2])
  i = i + 1
}

srtm_files = list.files(downloadFolder,pattern="srtm",full.names=TRUE)
srtm_list = lapply(srtm_files, raster)
rasters.mosaicargs = srtm_list
rasters.mosaicargs$fun = mean
alt_mosaic = do.call(mosaic, rasters.mosaicargs)
alt_mosaic = crop(alt_mosaic,extent(studyArea))
alt_mosaic = mask(alt_mosaic, water_bodies, inverse=TRUE) 
alt_mosaic = mask(alt_mosaic, studyArea)
alt_mosaic_proj = projectRaster(alt_mosaic,crs=sr_cea,method="bilinear")
names(alt_mosaic_proj) = "elevation"
writeRaster(alt_mosaic_proj,paste0(downloadFolder,"elevation.tif"),overwrite=T)
slope = terrain(alt_mosaic_proj,'slope',unit="degrees",filename=paste0(downloadFolder,"slope.tif"),overwrite=T)

# download two worldclim tiles, mosaic, rename variables, and crop/mask
worldclim1 = raster::getData('worldclim', var='bio', lon=-75, lat=15, res=0.5, path = downloadFolder)
worldclim2 = raster::getData('worldclim',var='bio', lon=-105, lat=15, res=0.5, path = downloadFolder)
worldclim = mosaic(worldclim1,worldclim2,fun=mean) # save this to file?
worldclim = crop(worldclim,extent(studyArea))
worldclim = mask(worldclim,studyArea,filename=paste0(downloadFolder,"worldclim.tif"))
#rename variables using this list http://www.worldclim.org/bioclim
worldclim.names = c("annual_mean_temp","mean_diurnal_range","isothermality","temp_seasonality",
                     "max_temp_warmest_month","min_temp_coldest_month","temp_annual_range",
                     "mean_temp_wettest_quarter","mean_temp_driest_quarter","mean_temp_warmest_quarter",
                     "mean_temp_coldest_quarter","annual_precip","precip_wettest_month","precip_driest_month",
                     "precip_seasonality","precip_wettest_quarter","precip_driest_quarter","precip_warmest_quarter",
                     "precip_coldest_quarter")
names(worldclim) = worldclim.names
worldclim = unstack(worldclim)
outputnames = paste(downloadFolder,"worldclim_",worldclim.names, ".tif",sep="") #write with worldclim at the front
for(i in seq_along(worldclim)){writeRaster(worldclim[[i]], file=outputnames[i], overwrite=T)}

# globcover landcover data
download.file("http://due.esrin.esa.int/files/Globcover2009_V2.3_Global_.zip",
              paste0(downloadFolder,"Globcover2009_V2.3_Global_.zip"))
unzip(paste0(downloadFolder,"Globcover2009_V2.3_Global_.zip"),exdir=substr(downloadFolder, 0, nchar(downloadFolder)-1))
globcover_tif = raster(paste0(downloadFolder,"GLOBCOVER_L4_200901_200912_V2.3.tif"))
globcover_tif_crop = crop(globcover_tif,studyArea)
forest_values = c(40, 50, 60, 70, 90, 100, 160, 170)
globcover_forest = globcover_tif_crop
for(value in forest_values){
  globcover_forest[globcover_forest == value] = 1}
globcover_forest[globcover_forest == 230] = NA
globcover_forest[globcover_forest > 1] = 0
writeRaster(globcover_forest,paste0(downloadFolder,"globcover_forest.tif"))

chelsa.files = list.files(paste0(downloadFolder,"CHELSA"),full.names=T,pattern="tif")
for(chelsa.file in chelsa.files){
  chelsa.raster = raster(chelsa.file)
  chelsa.crop = crop(chelsa.raster,studyArea)
  chelsa.mask = mask(chelsa.crop,studyArea,filename=paste0(downloadFolder,basename(chelsa.file)))
}

pspSl = as.psp(roads_merge)
px = pixellate(pspSl, eps=1000)
rLength = raster(px)
crs(rLength) = crs(roads_merge)
writeRaster(rLength,paste0(downloadFolder,"rLength.tif"))

roads_merge$roads = 1
roads_raster = raster(snap_raster)
roads_raster = rasterize(roads_merge,roads_raster,field="roads") #write this to file, takes a long time
roads_dist = distance(roads_raster)
roads_dist = mask(roads_dist, snap_raster)
names(roads_dist) = "road_distance"
writeRaster(roads_dist,paste0(downloadFolder,"road_distance.tif"),overwrite=T)

# downloaded manually
water_files = list.files("/Users/codyschank/Downloads/",pattern="occurrence",full.names=T)
water_files_list = lapply(water_files, raster)
rasters.mosaicargs = water_files_list
rasters.mosaicargs$fun = max
rasters.mosaicargs$filename = paste0(downloadFolder,"water_occurrence.tif")
water_occurrence = do.call(mosaic, rasters.mosaicargs)

### OSM testing

# http://download.geofabrik.de/
# download osm.bz2, unzip
src = osmsource_osmosis(file = "/Users/codyschank/Downloads/central-america-latest.osm")
box <- corner_bbox(-85,11,-84,12) # 2 minutes?
muc <- get_osm(box, src)
summary(muc$nodes)
summary(muc$ways)
summary(muc$relations)

find_ids = find(muc, node(tags(v == "hamlet")))
o1 = subset(muc, node_ids = find_ids)
hamlet_points = as_sp(o1, "points")
plot(hamlet_points)

find_ids = find(muc, relation(tags(k == "boundary")))
boundaries = lapply(find_ids,
                    function(i){
                      raw <- get_osm(relation(i), full = TRUE)
                      as_sp(raw, "lines")
                    })
for(i in seq(1,27)){
  plot(boundaries[[14]])
}

o1 = subset(muc, relation_ids = find_ids)
boundary_polygons = as_sp(o1, "polygons")
plot(hamlet_points)

find_ids = find(muc, node(tags(v == "protected_area")))
find_ids = find(muc, way(tags(v == "protected_area")))
find_ids = find(muc, node(tags(v %agrep% "protected"))) # park, reserve
find_ids

# nodes, ways, and relations
# key-value = k-v

