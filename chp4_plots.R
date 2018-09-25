###########################
### load libraries, etc ###
###########################

libs = c('ggplot2','reshape2','raster','rgdal','xtable','rgeos','igraph','gdistance',"classInt","RColorBrewer","dpylr")
lapply(libs, require, character.only = TRUE)

shadowtext <- function(x, y=NULL, labels, bg="white", col="black", theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}

sr_cea = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sr_geo = " +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#################
### data prep ###
#################

output.folder = "/Users/codyschank/Dropbox/Research/Dissertation/outputs_chp4_4km_final2/"

foldernames = list.dirs(path = output.folder,full.names = TRUE)

predictors.folder = "/Users/codyschank/Dropbox/Research/Dissertation/DataPrep_output/4km_chp4/"
plots.folder = paste0(output.folder,"plots/")
dir.create(file.path(plots.folder), showWarnings = FALSE)

csvs = list.files(foldernames,full.names=TRUE,pattern="outputs_df")
merged = do.call("rbind", lapply(csvs, read.csv, header=TRUE,stringsAsFactors=FALSE))
merged$candidate_model = paste("model",merged$candidate_model,sep="")

allData.df = data.frame(merged)

resolution = unique(allData.df$resolution)
sampling.radius = round(unique(allData.df$r.po),0)
n.sample = unique(allData.df$n.samples)
model = unique(allData.df$model)

xCovs = paste0("x.",na.omit(unique(unlist(strsplit(allData.df$x.po.covs,",")))))
beta0.xCovs = c("beta0",xCovs)

wCovs.po = paste0("w.po.",na.omit(unique(unlist(strsplit(allData.df$w.po.covs,",")))))
wCovs.pa = paste0("w.pa.",na.omit(unique(unlist(strsplit(allData.df$w.pa.covs,",")))))
wCovs.pa.sampling = paste0("w.pa.sampling.",na.omit(unique(unlist(strsplit(allData.df$w.pa.sampling.covs,",")))))
alpha0.wCovs = c("alpha0.po","alpha0.pa",wCovs.po,wCovs.pa,wCovs.pa.sampling)

allCovs = c(beta0.xCovs,alpha0.wCovs)
allCovs.se = paste0(allCovs,".se")

model.averaging.final = unique(allData.df[1,c("model","resolution","n.samples","r.po")]) # hard coded this in for now
model.averaging.final[,allCovs] = NA
model.averaging.final[,allCovs.se] = NA

# check recip and total population to see if I need to fine tune recip
minrecipCondNum = 1e-6
View(allData.df[,c("total_population","recip")])
sum(allData.df$recip>1e-6,na.rm=T)
sum(allData.df$recip>1e-5,na.rm=T)
sum(allData.df$recip>1e-4,na.rm=T)

new_minrecipCondNum = 1e-4 # as opposed to 1e-6
allData.df$estimated.se = allData.df$recip > new_minrecipCondNum

# remove outlier models for total population, changing recipCondNum doesn't get all of them
#allData.df$estimated.se = (allData.df$recip > minrecipCondNum) & (allData.df$total_population < 100000) 

model.averaging.final[, allCovs] = colMeans(allData.df[allData.df$estimated.se, allCovs], na.rm=TRUE)
model.averaging.final[, allCovs.se] = colMeans(allData.df[allData.df$estimated.se, allCovs.se], na.rm=TRUE)

all.covs = na.omit(unique(unlist(strsplit(allData.df$x.po.covs,","))))
sgrid = raster::stack(paste0(predictors.folder,all.covs,".tif"))
xycov = na.omit(raster::rasterToPoints(sgrid))
xy = xycov[,1:2]
des = as.data.frame(xycov[,3:ncol(xycov)])
X.back = as.matrix(cbind(rep(1, dim(des)[1]), des))

#d_cor = as.matrix(cor(des))
#d_cor_melt = arrange(melt(d_cor), -abs(value))
#d_cor_melt # twi and slope are correlated

##########################
### Spatial prediction ###
##########################

linearPredictor.sdm = X.back %*% as.numeric(model.averaging.final[,beta0.xCovs])
response.sdm.lambda = exp(linearPredictor.sdm)
raster.sdm.lambda = raster(sgrid)
cells = cellFromXY(raster.sdm.lambda, xy)
raster.sdm.lambda[cells] = response.sdm.lambda
raster.sdm.psi = 1-exp(-1*(raster.sdm.lambda)*resolution^2)

file_name_base = paste0(resolution,"km_",sampling.radius,"m_",n.sample,"samples_",model)[1] # hard coded for now

maps.folder = paste0(plots.folder,"maps/")
dir.create(file.path(maps.folder), showWarnings = FALSE)

tifs.folder = paste0(maps.folder,"tifs/")
dir.create(file.path(tifs.folder), showWarnings = FALSE)

writeRaster(raster.sdm.lambda, paste0(tifs.folder,file_name_base,"_lambda.tif"), "GTiff", overwrite=TRUE)
writeRaster(raster.sdm.psi, paste0(tifs.folder,file_name_base,"_psi.tif"), "GTiff", overwrite=TRUE)

#raster.abundance = raster.sdm.lambda*16
#plot(raster.abundance >= 1)

##########################
### Read in Map layers ###
##########################

input.folder = "/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_input/"

forestloss = raster(paste0(input.folder,"forestloss1km.tif"))
datamask = raster(paste0(input.folder,"datamask1km.tif"))
pa_raster = raster("/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_output/1km/protected.tif") # need to fix this to match 4km
pa_raster_4km = raster("/Users/codyschank/Dropbox/Research/Dissertation/dataPrep_output/4km_chp4/protected.tif")

studyArea = readOGR(paste0(input.folder,"tapir_studyArea_new.shp"),"tapir_studyArea_new")
studyArea_proj = spTransform(studyArea,sr_cea)

presenceData = read.csv("/Users/codyschank/Dropbox/Research/Dissertation/Presence_Data/presence_only_proj.csv")
paData = read.csv(paste0("/Users/codyschank/Dropbox/Research/Dissertation/Presence_Data/detection_histories_proj_10samples.csv"))
paData = paData[paData$Y > extent(studyArea_proj)[3],]

coordinates(presenceData) = c("X","Y")
proj4string(presenceData) = sr_cea
coordinates(paData) = c("X","Y")
proj4string(paData) = sr_cea

protected_areas = readOGR(paste0(input.folder,"WDPA_selection.shp"),"WDPA_selection", stringsAsFactors = F, use_iconv = TRUE, encoding = "UTF-8")
protected_areas_proj = spTransform(protected_areas,sr_cea)
lake = readOGR(paste0(input.folder,"Lake_proj.shp"),"Lake_proj")
countries = readOGR(paste0(input.folder,"country1m_vmap.shp"),"country1m_vmap")
e.buffer = extent(extent(studyArea)[1]-1,extent(studyArea)[2]+1,
                  extent(studyArea)[3]-1,extent(studyArea)[4]+1)
countries.crop.geo = crop(countries,e.buffer)
countries.crop = spTransform(countries.crop.geo,sr_cea)
canals = readOGR(paste0(input.folder,"canal_paths.shp"), layer = "canal_paths")

# testing some additional layers

mx_community_areas = readOGR(paste0(input.folder,"COMMUNITY_PA.shp"), layer = "COMMUNITY_PA")
mx_community_areas_proj = spTransform(mx_community_areas,sr_cea)

indigenous_presence = readOGR(paste0(input.folder,"presindigw.shp"), layer = "presindigw")
indigenous_presence$prcnt_indigenous = indigenous_presence$pobindi / indigenous_presence$pobtot

iucn_indigenous_areas = readOGR(paste0(input.folder,"PI_UICN_ORMACC_2016.shp"), layer = "PI_UICN_ORMACC_2016")
iucn_indigenous_areas_proj = spTransform(iucn_indigenous_areas,sr_cea)

#spplot(indigenous_presence,"prcnt_indigenous")

################
### Hotspots ###
################

intensity = raster.sdm.lambda
ones = raster(intensity)
ones = setValues(ones, 1) 
ones = mask(ones,intensity)

global.mean = cellStats(intensity, "mean")
global.sd = cellStats(intensity, "sd")
cell.nums = Which(!is.na(intensity), cells = TRUE)
all.n = length(intensity)

fw = focalWeight(intensity, 10000, type='circle') 
fw[fw>0] = 1 
sum.weights = focal(ones, w=fw, fun=sum, na.rm=TRUE)
sum.values = focal(intensity, w=fw, fun=sum, na.rm=TRUE)

localG.numerator = sum.values - (sum.weights*global.mean)
localG.denominator = global.sd*(((all.n*sum.weights)-sum.weights^2)/(all.n-1))^0.5
localG = localG.numerator/localG.denominator
localG = mask(localG,intensity)

writeRaster(localG, paste0(maps.folder,"localG.tif"), overwrite = TRUE)

cores = localG > 1.96
cores[cores == 0] = NA

cores.clump = clump(cores, directions = 4)

source("/Users/codyschank/Dropbox/Research/Dissertation/scripts/gdal_polygonizeR.R")
# requires python, and gdal to be correctly installed
cores.poly = gdal_polygonizeR(cores.clump)
cores.poly = gBuffer(cores.poly, byid=TRUE, width = 0) # repairs geometry for later
writeOGR(cores.poly, dsn=paste0(maps.folder,"cores.poly.shp"), layer="cores.poly", driver="ESRI Shapefile", overwrite_layer = TRUE)

plot(cores.poly) 

###############
### overlay ###
###############

# get population in each core
cores.raster = raster(intensity)
cores.raster = rasterize(cores.poly,cores.raster,field="DN")
cores.zonal = data.frame(zonal(intensity*resolution^2, cores.raster, fun = 'sum', na.rm = TRUE)) # multiply intensity by area of site
names(cores.zonal)[2] = "sum"
cores.poly = merge(cores.poly,cores.zonal,by.x="DN",by.y="zone")

cores.poly$area_sqkm = area(cores.poly) / 1000^2
cores.poly = cores.poly[cores.poly$sum>10,] # remove any cores with less than 10 individuals predicted

# redo indexes and DN after removing cores with less than 10 individuals
cores.poly$DN = seq(1,length(cores.poly))
for(i in seq(1,length(cores.poly))){cores.poly@polygons[[i]]@ID = as.character(cores.poly@data$DN[i])} # set ID to DN for each polygon

# calculate confirmed presence in cores
paData.presence = paData[paData$Presence == 1,]
cores.poly[cores.poly$DN %in% intersect(cores.poly,presenceData)$DN,"PO"] = 1
cores.poly[is.na(cores.poly$PO), "PO"] = 0
cores.poly[cores.poly$DN %in% intersect(cores.poly,paData.presence)$DN,"PA"] = 1
cores.poly[is.na(cores.poly$PA), "PA"] = 0
cores.poly[(cores.poly$PO == 1) | (cores.poly$PA == 1),"Presence"] = 1
cores.poly[is.na(cores.poly$Presence), "Presence"] = 0

# calculate forest loss
forestloss_pts = rasterToPoints(forestloss,spatial=TRUE)
forestloss_proj = spTransform(forestloss_pts,sr_cea)
forestloss_intersect = intersect(forestloss_proj,cores.poly)
forestloss_summary = aggregate(forestloss_intersect$forestloss1km,by=list(DN=forestloss_intersect$DN), FUN=sum)

datamask_pts = rasterToPoints(datamask,spatial=TRUE)
datamask_proj = spTransform(datamask_pts,sr_cea)
datamask_intersect = intersect(datamask_proj,cores.poly)
datamask_summary = aggregate(datamask_intersect$datamask1km,by=list(DN=datamask_intersect$DN), FUN=sum)

pa_pts = rasterToPoints(pa_raster,spatial=TRUE)
pa_intersect = intersect(pa_pts,cores.poly)
pa_summary = aggregate(pa_intersect$protected,by=list(DN=pa_intersect$DN), FUN=sum)

cores.poly$prcnt_FL = forestloss_summary$x / datamask_summary$x
cores.poly$FL_index = cores.poly$prcnt_FL/max(cores.poly$prcnt_FL)

cores.poly$pa_sqkm = pa_summary$x
cores.poly$prcnt_prot = cores.poly$pa_sqkm/ cores.poly$area_sqkm
cores.poly$prcnt_unprot = 1 - cores.poly$prcnt_prot

cores.poly$popindex = cores.poly$sum/max(cores.poly$sum)
cores.poly$inv_popindex = 1-cores.poly$popindex

######################################################################
### forest loss in each protected area that intersects with a core ###
######################################################################

pas_cores_intersect = intersect(protected_areas_proj,cores.poly)
pas_cores_select = protected_areas_proj[protected_areas_proj$WDPAID %in% unique(pas_cores_intersect$WDPAID),]
keeps = c("WDPAID")
pas_cores_select = pas_cores_select[ , (names(pas_cores_select) %in% keeps)]

wdpaids = unique(pas_cores_select$WDPAID)
FL_pas = data.frame(WDPAID = wdpaids)
FL_pas$FL = NA
for(wdpaid in wdpaids){
  pas_cores_select_temp = pas_cores_select[pas_cores_select$WDPAID == wdpaid,]
  
  possibleError.occu = tryCatch(
    expr = (forestloss_intersect_pas = intersect(forestloss_proj,pas_cores_select_temp)),
    error=function(e) e
  )
  
  if(!inherits(possibleError.occu, "error")){
    datamask_pas_intersect = intersect(datamask_proj,forestloss_intersect_pas)
    FL_pa = sum(forestloss_intersect_pas$forestloss1km,na.rm=TRUE)/sum(datamask_pas_intersect$datamask1km,na.rm=TRUE)
    FL_pas[FL_pas$WDPAID==wdpaid,"FL"] = FL_pa
  }else{FL_pas[FL_pas$WDPAID==wdpaid,"FL"] = NA}
  
}

FL_pas = FL_pas[!is.na(FL_pas$FL),]

pas_cores_join = merge(protected_areas_proj,FL_pas,by="WDPAID", all=FALSE)
pas_cores_join = merge(pas_cores_join,aggregate(DN ~ WDPAID, pas_cores_intersect, paste, collapse = ", "), by="WDPAID")
pas_cores_join = pas_cores_join[pas_cores_join$FL > 0.05,]
pas_cores_join$area_sqkm = area(pas_cores_join) / 1000^2
pas_cores_join = pas_cores_join[order(pas_cores_join$FL, decreasing = T),]

#############################
### Connectivity analysis ###
#############################

# manually do different dispersal distances by changing the first number in the calc below, and saving flow.matrix as flow.matrix.1000, for example

area.home.range = 16*1000000 
linear.dimension.home.range = area.home.range^0.5
max.dispersal.dist = 40*linear.dimension.home.range # from Bowman et al 2002
median.dispersal.dist = 7*linear.dimension.home.range  # from Bowman et al 2002
mean.dispersal = median.dispersal.dist/1000 # in km, for some reason this seems better for dispersal kernel models below

# dispersal kernels from Nathan et al 2012
a.neg.exp = mean.dispersal/2
neg.exp.dispersal = function(x) (1/(2*pi*a.neg.exp^2))*exp(-1*(x/a.neg.exp))
png(paste0(plots.folder,"negative_exponential_curve_",linear.dimension.home.range,".png"))
curve(neg.exp.dispersal,from=0,to=mean.dispersal*2) 
dev.off()

a.guassian = (mean.dispersal*2)/pi^0.5
gaussian.dispersal = function(x) (1/(pi*a.guassian^2))*exp(-1*(x^2/a.guassian^2))
curve(gaussian.dispersal,from=0,to=mean.dispersal*2) 

d.all.cores.countries = gDistance(cores.poly, byid=T) # distances are in meters
flow.matrix = (d.all.cores.countries<max.dispersal.dist)*1
diag(flow.matrix) = 0

transition.layer = transition(intensity, transitionFunction=mean, directions=16, symm=TRUE)

abundance = intensity*resolution^2
abundance.pts = rasterToPoints(abundance,spatial=TRUE)

for(DN1 in cores.poly$DN){
  for(DN2 in cores.poly$DN){
    if(flow.matrix[DN1,DN2] == 1){
      print(paste0(DN1,"->",DN2," start"))
      core.select1 = cores.poly[cores.poly$DN==DN1,]
      core.select1.intersect = intersect(abundance.pts,core.select1)
      core.select2 = cores.poly[cores.poly$DN==DN2,]
      core.select2.intersect = intersect(abundance.pts,core.select2)
      nearestPoints = gNearestPoints(core.select1, core.select2)
      nearest.buffer = gBuffer(nearestPoints[2,], width=max.dispersal.dist*1.10, byid=TRUE)
      #core.select1.crop = intersect(core.select1.intersect,nearest.buffer) # subset only those within max dispersal 
      core.select1.crop = core.select1.intersect[which(over(core.select1.intersect,nearest.buffer) == 1),] # bug in raster package, https://stackoverflow.com/questions/45614227/very-likely-bug-in-the-r-raster-package-intersect-function-when-intersecting-1-p
      flux = 0
      if(nrow(core.select1.crop)>0){
        for(i in seq(1,nrow(core.select1.crop))){
          
          nearest.buffer = gBuffer(core.select1.crop[i,], width=max.dispersal.dist*1.25, byid=TRUE)
          core.select2.crop = intersect(core.select2.intersect,nearest.buffer) # subset only those within max dispersal 

          shortestPath.lines = shortestPath(transition.layer, core.select1.crop[i,], core.select2.crop, "SpatialLines")
          shortestPath.distances = SpatialLinesLengths(shortestPath.lines)
          if(sum(shortestPath.distances<max.dispersal.dist)>0){
            shortestPath.distances.select = shortestPath.distances[shortestPath.distances<max.dispersal.dist]
            #flux = flux + sum((core.select1.crop@data[i,"layer"]^0.5)*(neg.exp.dispersal(shortestPath.distances.select/1000)/length(shortestPath.distances.select))) # added normalization by dividing by number of LCPs
            flux = flux + sum((core.select1.crop@data[i,"layer"]^0.5)*(neg.exp.dispersal(shortestPath.distances.select/1000))) # removed normalization, on second thought doesn't make sense
          }
        }
      }
      flow.matrix[DN1,DN2] = flux
      print(paste0(DN1,"->",DN2," complete"))
    }
  }
}

#flow.matrix.1000 = flow.matrix

write.csv(flow.matrix.1000,paste0(plots.folder,"flow.matrix.1000.csv"))
write.csv(flow.matrix.2000,paste0(plots.folder,"flow.matrix.2000.csv"))
write.csv(flow.matrix.4000,paste0(plots.folder,"flow.matrix.4000.csv"))

####################
### Map of cores ###
####################

core.centroids = gCentroid(cores.poly, byid=TRUE, id=cores.poly$DN)

# where are the points of highest intensity in Integrated Model
intensity.matrix = values(intensity)
max.index = which.max(intensity.matrix)
max.location = xyFromCell(intensity, max.index)

guanacaste = coordinates(gCentroid(protected_areas_proj[protected_areas_proj$NAME=="Guanacaste",]))
la_muralla = coordinates(gCentroid(protected_areas_proj[protected_areas_proj$NAME=="La Muralla",]))
volcan_barva = coordinates(gCentroid(protected_areas_proj[protected_areas_proj$NAME=="Braulio Carrillo",]))
la_sepultura = coordinates(gCentroid(protected_areas_proj[protected_areas_proj$NAME=="La Sepultura",]))
el_triunfo = coordinates(gCentroid(protected_areas_proj[protected_areas_proj$NAME=="El Triunfo",]))

# convert pa raster layer to polygons so I can use instead of protect_areas_proj
pa_raster_copy = pa_raster
pa_raster_copy[pa_raster_copy==0] = NA
pa_raster_polygonize = gdal_polygonizeR(pa_raster_copy)

pdf(paste0(maps.folder,"cores_map.pdf"), width = 5, height = 5)
par(mai=c(0.0,0.0,0.0,0.0))
plot(countries.crop, border = FALSE)
e = extent(countries.crop)
rect(e[1], e[3], e[2], e[4], col = alpha("light blue", 0.3), border = FALSE)
plot(countries.crop, col = alpha("tan", 0.3), border = FALSE, add = TRUE)
plot(studyArea_proj,col = alpha("tan", 0.6), border = "black", lwd = 0.25, add = TRUE)
plot(pa_raster_polygonize,col=alpha("brown", 0.6),border=FALSE, add=TRUE) # could use pas_cores_select instead
plot(cores.poly, col=alpha("dark green", 0.6), border = FALSE, add=TRUE)
plot(lake, col=alpha("light blue", 0.6), border = FALSE, add = TRUE)
plot(countries.crop, border = "black", lwd = 0.25,  add = TRUE)
plot(canals, col= "purple", lwd = 2, add = TRUE)
shadowtext(core.centroids@coords[,1], core.centroids@coords[,2], labels=cores.poly$DN, cex=0.5, r=0.05)
shadowtext(guanacaste, labels="(a)", cex=0.75, r=0.05)
shadowtext(la_sepultura, labels="(b)", cex=0.75, r=0.05)
shadowtext(el_triunfo, labels="(c)", cex=0.75, r=0.05)
shadowtext(la_muralla, labels="(d)", cex=0.75, r=0.05)
#shadowtext(volcan_barva, labels="(d)", cex=0.75, r=0.05)
dev.off()

png(paste0(maps.folder,"cores_map.png"), width=5, height=5, units="in", res=300) 
par(mai=c(0.0,0.0,0.0,0.0))
plot(countries.crop, border = FALSE)
e = extent(countries.crop)
rect(e[1], e[3], e[2], e[4], col = alpha("light blue", 0.3), border = FALSE)
plot(countries.crop, col = alpha("tan", 0.3), border = FALSE, add = TRUE)
plot(studyArea_proj,col = alpha("tan", 0.6), border = "black", lwd = 0.25, add = TRUE)
plot(pa_raster_polygonize,col=alpha("brown", 0.6),border=FALSE, add=TRUE) # could use pas_cores_select instead
plot(cores.poly, col=alpha("dark green", 0.6), border = FALSE, add=TRUE)
plot(lake, col=alpha("light blue", 0.6), border = FALSE, add = TRUE)
plot(countries.crop, border = "black", lwd = 0.25,  add = TRUE)
plot(canals, col= "purple", lwd = 2, add = TRUE)
shadowtext(core.centroids@coords[,1], core.centroids@coords[,2], labels=cores.poly$DN, cex=0.5, r=0.05)
shadowtext(guanacaste, labels="(a)", cex=0.75, r=0.05)
shadowtext(la_sepultura, labels="(b)", cex=0.75, r=0.05)
shadowtext(el_triunfo, labels="(c)", cex=0.75, r=0.05)
shadowtext(la_muralla, labels="(d)", cex=0.75, r=0.05)
#shadowtext(volcan_barva, labels="(d)", cex=0.75, r=0.05)
dev.off()

# zoomed in map
DN1 = 23
DN2 = 21
i = 1
core.select1 = cores.poly[cores.poly$DN==DN1,]
core.select1.intersect = intersect(abundance.pts,core.select1)
core.select2 = cores.poly[cores.poly$DN==DN2,]
core.select2.intersect = intersect(abundance.pts,core.select2)
nearestPoints = gNearestPoints(core.select1, core.select2)
nearest.buffer = gBuffer(nearestPoints[2,], width=max.dispersal.dist*1.10, byid=TRUE)
core.select1.crop = intersect(core.select1.intersect,nearest.buffer) # subset only those within max dispersal
nearest.buffer = gBuffer(core.select1.crop[i,], width=max.dispersal.dist*1.25, byid=TRUE)
core.select2.crop = intersect(core.select2.intersect,nearest.buffer) # subset only those within max dispersal 
shortestPath.lines = shortestPath(transition.layer, core.select1.crop[i,], core.select2.crop, "SpatialLines")
shortestPath.distances = SpatialLinesLengths(shortestPath.lines)

min(shortestPath.distances) # minimum distance to get between these two cores
core.select1.crop@data[i,"layer"] # estimated N in the origina cell, not lambda
cores.poly.select = cores.poly[cores.poly$DN %in% c(DN1,DN2),]

e = as(extent(-9250064,-8780445,789027.4,1099974), 'SpatialPolygons')  
proj4string(e) = proj4string(cores.poly.select)
pa_raster_polygonize_crop = crop(pa_raster_polygonize,e)
cores.poly.crop = crop(cores.poly,e)
countries.crop.temp = crop(countries.crop,e)
canal.crop.temp = crop(canals,e)
intensity.crop = crop(intensity,e)
shortestPath.lines.crop = crop(shortestPath.lines,e)
DN1.label.loc = coordinates(gCentroid(crop(core.select1,e)))
DN2.label.loc = coordinates(gCentroid(crop(core.select2,e)))

breaks = c(0,0.01,0.02,0.04,0.08,0.16,0.32)
breaks = round(breaks,2)
rng = breaks[c(1,4,5,6,7)]
arg = list(at=rng, labels=rng, cex.axis=0.75)
palette = brewer.pal(n = 6, name = "YlOrBr")

pdf(paste0(maps.folder,"map_zoom_lcps.pdf"), width = 6, height = 4)
par(mai=c(0.0,0.0,0.0,0.0))
plot(e, col="light blue", border=F)
plot(intensity.crop, breaks=breaks, col=palette, axis.args=arg,  axes=FALSE, box=F, legend =F, add=T)
plot(cores.poly.crop, lwd = 1, add=TRUE) #col=alpha("dark green", 0.6)
plot(countries.crop.temp, border = "black", lwd = 0.25,  add = TRUE)
plot(canal.crop.temp, col= "purple", lwd = 2, add = TRUE)
plot(core.select1.crop, add=TRUE, pch=19, cex=0.3)
plot(shortestPath.lines.crop, col="blue", add=TRUE)
points(core.select1.crop[i,],col="blue",pch=19)
shadowtext(DN1.label.loc, labels=DN1, cex=1, r=0.1)
shadowtext(DN2.label.loc, labels=DN2, cex=1, r=0.1)
par(cex = 0.7)
# x1, x2, y1, y2
plot(intensity.crop, breaks=breaks, col=palette, axis.args=arg, axes=FALSE, box=F, legend.only=T, 
     smallplot=c(.85,.88, .1,.3))
dev.off()

png(paste0(maps.folder,"map_zoom_lcps.png"), width = 6, height = 4, units="in", res=300) 
par(mai=c(0.0,0.0,0.0,0.0))
plot(e, col="light blue", border=F)
plot(intensity.crop, breaks=breaks, col=palette, axis.args=arg,  axes=FALSE, box=F, legend =F, add=T)
plot(cores.poly.crop, lwd = 1, add=TRUE) #col=alpha("dark green", 0.6)
plot(countries.crop.temp, border = "black", lwd = 0.25,  add = TRUE)
plot(canal.crop.temp, col= "purple", lwd = 2, add = TRUE)
plot(core.select1.crop, add=TRUE, pch=19, cex=0.3)
plot(shortestPath.lines.crop, col="blue", add=TRUE)
points(core.select1.crop[i,],col="blue",pch=19)
shadowtext(DN1.label.loc, labels=DN1, cex=1, r=0.1)
shadowtext(DN2.label.loc, labels=DN2, cex=1, r=0.1)
par(cex = 0.7)
# x1, x2, y1, y2
plot(intensity.crop, breaks=breaks, col=palette, axis.args=arg, axes=FALSE, box=F, legend.only=T, 
     smallplot=c(.85,.88, .1,.3))
dev.off()

#########################
### Calculate Indices ###
#########################

flow.matrices = list(flow.matrix.1000,flow.matrix.2000,flow.matrix.4000)
linear.dimensions = c(1000,2000,4000)

for(i in seq(1,3)){
  
  flow.matrix = flow.matrices[[i]]
  linear.dimension = linear.dimensions[i]
  
  flow.matrix[flow.matrix==0] = NA
  diag(flow.matrix) = NA
  
  my.graph = graph_from_adjacency_matrix(flow.matrix, weighted=TRUE, mode="directed")
  my.graph = delete.edges(my.graph, which(is.na(E(my.graph)$weight))) 
  V(my.graph)$name = unique(cores.poly$DN)
  cores.poly@data[,paste0("between",linear.dimension)] = betweenness(my.graph)/max(betweenness(my.graph))
  
  # calculate farness by hand
  for(v in V(my.graph)$name){
    ds = distances(my.graph,v,mode="in")
    cores.poly[cores.poly$DN==v,paste0("farness",linear.dimension)] = sum(ds*is.finite(ds),na.rm=TRUE)
  }
  cores.poly@data[,paste0("farness",linear.dimension)] = cores.poly@data[,paste0("farness",linear.dimension)] / max(cores.poly@data[,paste0("farness",linear.dimension)])
  
  #test resilience of network by removing vertices one at a time to see which is most important to overall connectivity
  for(core in cores.poly$DN){
    my.graph.temp = delete.vertices(my.graph,as.character(core))
    cores.poly[cores.poly$DN == core, paste0("diameter",linear.dimension)] = diameter(my.graph.temp)
    cores.poly[cores.poly$DN == core, paste0("components",linear.dimension)] = components(my.graph.temp)$no
    cores.poly[cores.poly$DN == core, paste0("mean_distance",linear.dimension)] = mean_distance(my.graph.temp)
  }
  
  cores.poly@data[,paste0("diameter",linear.dimension)] = cores.poly@data[,paste0("diameter",linear.dimension)] / max(cores.poly@data[,paste0("diameter",linear.dimension)])
  cores.poly@data[,paste0("mean_distance",linear.dimension)] = cores.poly@data[,paste0("mean_distance",linear.dimension)] / max(cores.poly@data[,paste0("mean_distance",linear.dimension)])
  
  cores.poly@data[,paste0("components_created",linear.dimension)] = cores.poly@data[,paste0("components",linear.dimension)]  - components(my.graph)$no
  cores.poly@data[cores.poly@data[,paste0("components_created",linear.dimension)] < 0, paste0("components_created",linear.dimension)]=0
  cores.poly@data[,paste0("components_rank",linear.dimension)] = cores.poly@data[,paste0("components_created",linear.dimension)] / max(cores.poly@data[,paste0("components_created",linear.dimension)])
  
  cores.poly@data[,paste0("import",linear.dimension)] = cores.poly$popindex + cores.poly$Presence + cores.poly@data[,paste0("between",linear.dimension)] + cores.poly@data[,paste0("components_rank",linear.dimension)] 
  cores.poly@data[,paste0("vuln",linear.dimension)] = cores.poly$FL_index + cores.poly$prcnt_unprot + cores.poly$inv_popindex + cores.poly@data[,paste0("farness",linear.dimension)]
  
  # rescale importance and vulnerability to be from 0 to 1
  cores.poly@data[,paste0("import",linear.dimension)] = cores.poly@data[,paste0("import",linear.dimension)]/max(cores.poly@data[,paste0("import",linear.dimension)])
  cores.poly@data[,paste0("vuln",linear.dimension)] = cores.poly@data[,paste0("vuln",linear.dimension)]/max(cores.poly@data[,paste0("vuln",linear.dimension)])

  cores.poly@data[,paste0("comb",linear.dimension)]= cores.poly@data[,paste0("import",linear.dimension)] + cores.poly@data[,paste0("vuln",linear.dimension)]  
  
}

View(cores.poly@data)

######################
### Map of network ###
######################

for(i in seq(1,3)){
  
  flow.matrix = flow.matrices[[i]]
  
  my.graph = graph_from_adjacency_matrix(flow.matrix, weighted=TRUE, mode="directed")
  my.graph = delete.edges(my.graph, which(is.na(E(my.graph)$weight))) 
  
  vertex.weight = cores.poly$sum
  vertex.weight[vertex.weight<300]=300
  vertex.weight = round(vertex.weight,1)
  vertex.weight = vertex.weight*5000
  
  edge.weight = E(my.graph)$weight/max(E(my.graph)$weight)
  edge.weight = (edge.weight/max(edge.weight))*5
  # rescale to between 0.5 and 5
  edge.weight[edge.weight<0.5]=0.5
  edge.weight[edge.weight>10]=10
  
  L= layout_with_fr(my.graph)
  L[,1]=coordinates(core.centroids)[,1]
  L[,2]=coordinates(core.centroids)[,2]
  
  linear.dimension = linear.dimensions[i]
  
  pdf(paste0(maps.folder,"network_map",linear.dimension,".pdf"), width = 5, height = 5)
  par(mai=c(0.0,0.0,0.0,0.0))
  plot(countries.crop, border = FALSE)
  e = extent(countries.crop)
  rect(e[1], e[3], e[2], e[4], col = alpha("light blue", 0.3), border = FALSE)
  plot(countries.crop, col = alpha("tan", 0.3), border = FALSE, add = TRUE)
  plot(studyArea_proj,col = alpha("tan", 0.6), border = "black", lwd = 0.25, add = TRUE)
  plot(cores.poly, col=alpha("dark green", 0.6), border = FALSE, add=TRUE)
  plot(lake, col=alpha("light blue", 0.6), border = FALSE, add = TRUE)
  plot(countries.crop, border = "black", lwd = 0.25,  add = TRUE)
  plot(canals, col= "purple", lwd = 2, add = TRUE)
  plot(my.graph, layout=L, vertex.label=NA, edge.label=NA, vertex.frame.color = NA, vertex.size=vertex.weight, edge.color="black",edge.width=edge.weight,
       add = TRUE, rescale = FALSE, edge.curved=TRUE, edge.arrow.mode=0 ) 
  dev.off()
  
  png(paste0(maps.folder,"network_map",linear.dimension,".png"), width = 5, height = 5, units="in", res=300) 
  par(mai=c(0.0,0.0,0.0,0.0))
  plot(countries.crop, border = FALSE)
  e = extent(countries.crop)
  rect(e[1], e[3], e[2], e[4], col = alpha("light blue", 0.3), border = FALSE)
  plot(countries.crop, col = alpha("tan", 0.3), border = FALSE, add = TRUE)
  plot(studyArea_proj,col = alpha("tan", 0.6), border = "black", lwd = 0.25, add = TRUE)
  plot(cores.poly, col=alpha("dark green", 0.6), border = FALSE, add=TRUE)
  plot(lake, col=alpha("light blue", 0.6), border = FALSE, add = TRUE)
  plot(countries.crop, border = "black", lwd = 0.25,  add = TRUE)
  plot(canals, col= "purple", lwd = 2, add = TRUE)
  plot(my.graph, layout=L, vertex.label=NA, edge.label=NA, vertex.frame.color = NA, vertex.size=vertex.weight, edge.color="black",edge.width=edge.weight,
       add = TRUE, rescale = FALSE, edge.curved=TRUE, edge.arrow.mode=0 ) 
  dev.off()
  
  i = i + 1
  
}


##############################
### Make latex figure file ###
##############################

figurefile = file(paste0(plots.folder,"figures.tex"), open = "wt")
write("\\documentclass{article}",figurefile)
write("\\usepackage{graphicx}",figurefile)
write("\\usepackage{subcaption}",figurefile)
write("\\usepackage{nopageno}",figurefile)
write("\\usepackage[font=small,labelfont=bf]{caption}",figurefile)
write("\\usepackage[letterpaper, total={8in, 8in}]{geometry}",figurefile)
write("\\begin{document}",figurefile)

write("\\begin{figure}",figurefile)
write(paste0("\\includegraphics[width=\\textwidth,height=\\textheight,keepaspectratio]{",maps.folder,"map_zoom_lcps","}"),figurefile)
caption = paste0("Example of flow calculation for one source pixel going from patch ",DN1," to patch ",DN2,". The source pixel (blue dot) contains ", round(core.select1.crop@data[i,"layer"],3), " individuals. Black dots are additional pixels that fall within range of the destination patch, and will subsequently have least cost paths (LCPs) and flow calculated. Blue lines are LCPs, purple line is the Panama Canal, estimated tapir density (individuals/km2) shown in shades of brown. ")
write(paste0("\\caption{",caption,"}"),figurefile)
write("\\end{figure}",figurefile)

write("\\begin{figure}",figurefile)
write(paste0("\\includegraphics[width=\\textwidth,height=\\textheight,keepaspectratio]{",plots.folder,"EVI_treecover_response","}"),figurefile)
caption = "3D response plot of intensity to EVI and forest cover."
write(paste0("\\caption{",caption,"}"),figurefile)
write("\\end{figure}",figurefile)

write("\\begin{figure}",figurefile)
write(paste0("\\includegraphics[width=\\textwidth,height=\\textheight,keepaspectratio]{",maps.folder,"cores_map","}"),figurefile)
caption = "Map of tapir habitat patches, across their geographic range, with an estimated population more than 10 individuals, identified using an Integrated Species Distribution Model (ISDM). Habitat patches are numbered 1-28 and shown in green, protected areas shown in brown, purple lines represent the Panama Canal and planned path of the Nicaragua Canal. Specific locations mentioned in the Discussion: (a) - Santa Rosa, (b) - La Sepultura, (c) - El Triunfo, (d) - La Muralla."
#caption = "Map of habitat patches. Green areas are habitat patches. Brown areas are protected areas. Purple lines represent the Panama Canal and the planned path of the Nicaragua Canal. (a) - Santa Rosa, (b) - La Sepultura, (c) - El Triunfo, (d) - La Muralla."
write(paste0("\\caption{",caption,"}"),figurefile)
write("\\end{figure}",figurefile)

write("\\begin{figure}",figurefile)
write("\\begin{subfigure}{.33\\textwidth}",figurefile)
write(paste0("\\includegraphics[width=7cm,keepaspectratio]{",maps.folder,"network_map1000","}"),figurefile)
write("\\caption{Max dispersal = 40 km}",figurefile)
write("\\end{subfigure}",figurefile)
write("\\begin{subfigure}{.33\\textwidth}",figurefile)
write(paste0("\\includegraphics[width=7cm,keepaspectratio]{",maps.folder,"network_map2000","}"),figurefile)
write("\\caption{Max dispersal = 80 km}",figurefile)
write("\\end{subfigure}",figurefile)
write("\\begin{subfigure}{.33\\textwidth}",figurefile)
write(paste0("\\includegraphics[width=7cm,keepaspectratio]{",maps.folder,"network_map4000","}"),figurefile)
write("\\caption{Max dispersal = 160 km}",figurefile)
write("\\end{subfigure}",figurefile)
caption = "Maps of network connectivity. Thickness of edges indicates strength of flow between vertices, size of nodes indicates estimated population size."
write(paste0("\\caption{",caption,"}"),figurefile)
write("\\end{figure}",figurefile)

write("\\end{document}",figurefile)
closeAllConnections()


###################
### Make tables ###
###################

cores.poly@data$sum = round(cores.poly@data$sum,0)
final.indices = cores.poly@data[,c("DN","Presence","sum","between4000","components_rank4000","import4000","prcnt_unprot","FL_index","farness4000","vuln4000","comb4000")]
final.indices = final.indices[order(final.indices$comb4000,decreasing=TRUE),]
names(final.indices) = c("ID","Occup.","Pop.","Between","Components","Import.","% Unprot.","Forest loss","Farness","Vuln.","Combined")

tablefile = file(paste0(plots.folder,"tables.tex"), open = "wt")
write("\\documentclass{article}\n",tablefile)
write("\\usepackage[margin=1in]{geometry}\n",tablefile)
write("\\usepackage[width=.75\\textwidth]{caption}\n",tablefile)
write("\\usepackage[utf8]{inputenc}\n",tablefile)
write("\\begin{document}\n",tablefile)
sink(tablefile)

coefficients.table.ests = model.averaging.final[,allCovs]
coefficients.table.se = model.averaging.final[,allCovs.se]
coefficients.table = cbind(t(coefficients.table.ests),t(coefficients.table.se))
coefficients.table = data.frame(coefficients.table)
coefficients.table[,3] = abs(coefficients.table[,1]/coefficients.table[,2])
coefficients.table[,4] = 2*pnorm(-abs(coefficients.table[,3]))
names(coefficients.table) = c("Estimate","SE","z-score","p-value")
row.names(coefficients.table) = gsub(".pa",".so",row.names(coefficients.table))

table.caption = "Coefficient estimtes, standard errors, z-scores, and p-values"
coefficients.table.xtable = xtable(coefficients.table, caption = table.caption)

sink(tablefile)
print(coefficients.table.xtable, type="latex", hline.after=c(0),include.rownames=T, comment=F, caption.placement = "bottom")
write.csv(coefficients.table,paste0(plots.folder,"coefficients.table.csv"))

table.caption = "Vulernability, importance, and combined indices for conservation prioritization in the maximum dispersal = 160 km scenario."
indices.table.xtable = xtable(final.indices, caption = table.caption)

sink(tablefile)
print(indices.table.xtable, type="latex", hline.after=c(0), include.rownames=F, comment=F, caption.placement = "bottom")
write.csv(final.indices,paste0(plots.folder,"final.indices.csv"))

PA_FL_table = pas_cores_join@data[,c("DN","NAME","FL","area_sqkm","SUB_LOC")] # ,"IUCN_CAT"
names(PA_FL_table) = c("Patch ID","Name","Forest Loss","Area (sq km)","Location") # ,"IUCN category"
table.caption = "Forest loss rates (2000-2014), protected status, area, and location of protected areas that intersect tapir habitat patches."
PA_FL_table.xtable = xtable(PA_FL_table, caption = table.caption)

sink(tablefile)
print(PA_FL_table.xtable, type="latex", hline.after=c(0), include.rownames=F, comment=F, caption.placement = "bottom")
write.csv(PA_FL_table,paste0(plots.folder,"PA_FL_table.csv"))

write("\\end{document}",tablefile)
closeAllConnections() 

########################
### 3d response plot ###
########################

# x is EVI
min.x = min(des$EVI,na.rm=T)
max.x = max(des$EVI,na.rm=T)
all.x = seq(min.x,max.x,0.1)

# y is Treecover
min.y = min(des$treecover2000,na.rm=T)
max.y = max(des$treecover2000,na.rm=T)
all.y = seq(min.y,max.y,0.1)
all.x = all.y # persp() only works when x and y are the same length?

all.z = data.frame(matrix(nrow=length(all.x), ncol=length(all.y)))
xCoefs = model.averaging.final[,c("beta0",xCovs)]

k = 1
for(x in seq(1,length(all.x))){
  for(y in seq(1,length(all.y))){
    newData = c(1,colMeans(des,na.rm=T))
    newData[grep("EVI",names(newData))[1]] = all.x[x]
    newData[grep(paste0("EVI_sq"),names(newData))[1]] = all.x[x]^2
    newData[grep("treecover2000",names(newData))[1]] = all.y[y]
    newData[grep("evi_treecover2000",names(newData))[1]] = all.x[x]*all.y[y]
    prediction = exp(as.matrix(xCoefs) %*% as.matrix(newData))
    all.z[x,y]=prediction
    print(k/(length(all.x)*length(all.y)))
    k=k+1
  }
}

z = as.matrix(all.z)
nrz <- nrow(z)
ncz <- ncol(z)
nb.col=100
color <- topo.colors(nb.col)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nb.col)
png(paste0(plots.folder,"EVI_treecover_response.png"), height=6, width=6, units="in",res=300)
persp(all.x, all.x, z=z, phi = 45, theta = 45, col = color[facetcol],
      xlab = "EVI", ylab="Treecover", zlab = "lambda")
dev.off()


# x is chelsa_precip_driest_quarter
min.x = min(des$chelsa_precip_driest_quarter,na.rm=T)
max.x = max(des$chelsa_precip_driest_quarter,na.rm=T)
all.x = seq(min.x,max.x,0.1)

# y is chelsa_mean_temp_driest_quarter
min.y = min(des$chelsa_mean_temp_driest_quarter,na.rm=T)
max.y = max(des$chelsa_mean_temp_driest_quarter,na.rm=T)
all.y = seq(min.y,max.y,0.1)
all.x = all.y # persp() only works when x and y are the same length?

all.z = data.frame(matrix(nrow=length(all.x), ncol=length(all.y)))
xCoefs = model.averaging.final[,c("beta0",xCovs)]

k = 1
for(x in seq(1,length(all.x))){
  for(y in seq(1,length(all.y))){
    newData = c(1,colMeans(des,na.rm=T))
    newData[grep("chelsa_precip_driest_quarter",names(newData))[1]] = all.x[x]
    newData[grep(paste0("chelsa_precip_driest_quarter_sq"),names(newData))[1]] = all.x[x]^2
    newData[grep("chelsa_mean_temp_driest_quarter",names(newData))[1]] = all.y[y]
    newData[grep(paste0("chelsa_mean_temp_driest_quarter"),names(newData))[1]] = all.y[y]^2
    newData[grep("precip_temp",names(newData))[1]] = all.x[x]*all.y[y]
    prediction = exp(as.matrix(xCoefs) %*% as.matrix(newData))
    all.z[x,y]=prediction
    print(k/(length(all.x)*length(all.y)))
    k=k+1
  }
}

z = as.matrix(all.z)
nrz <- nrow(z)
ncz <- ncol(z)
nb.col=100
color <- topo.colors(nb.col)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nb.col)
png(paste0(plots.folder,"precip_temp_response.png"), height=6, width=6, units="in",res=300)
persp(all.x, all.x, z=z, phi = 45, theta = 135, col = color[facetcol],
      xlab = "precip", ylab="temp", zlab = "lambda")
dev.off()


#######################
### important stats ###
#######################

### Estimated intensity and patches section

# mean intensity from model
cellStats(intensity, "mean") 
# max intensity from model
cellStats(intensity, "max") 
# percent of land covered by these cores
sum(cores.poly@data[cores.poly@data$sum>10,"area_sqkm"])/(cellStats(!is.na(intensity), "sum")*resolution^2) 
# total individuals in these cores
sum(cores.poly@data[cores.poly@data$sum>10,"sum"])
# % of total population from individuals in these cores
sum(cores.poly@data[cores.poly@data$sum>10,"sum"]) / (cellStats(intensity,"sum")*resolution^2)

### Network stats section

# number of patches causing breaks
sum(cores.poly@data$components_created > 0)

### Indices section

View(cores.poly)
# vulnerability components for most vulnerable core
cores.poly[which(cores.poly$vuln4000==max(cores.poly$vuln4000)),c("DN","sum","prcnt_prot","prcnt_FL","farness4000")]
# mean population in certain vulnerable but unimportant cores
mean(cores.poly.select.low1$sum)
cores.poly.select.low1$Presence
cores.poly.select.low1$between
# population in most important core
cores.poly[which(cores.poly$import==max(cores.poly$import)),c("DN","sum")]
# percent difference between top combined score and next one
(cores.poly@data[order(cores.poly$comb, decreasing=T),"comb"][2] - cores.poly@data[order(cores.poly$comb, decreasing=T),"comb"][1])/cores.poly@data[order(cores.poly$comb, decreasing=T),"comb"][2]

plot(cores.poly$import,cores.poly$vuln)

### Protection stats of patches

# percent of cores under some kind of protection
sum(cores.poly$pa_sqkm)/sum(cores.poly$area_sqkm)

# connectivity stuff

my.graph = graph_from_adjacency_matrix(flow.matrix.1000, weighted=TRUE, mode="directed")
components(my.graph)$no
my.graph = graph_from_adjacency_matrix(flow.matrix.4000, weighted=TRUE, mode="directed")
components(my.graph)$no

### Indigenous overlap stats

# calculate protected area coverage in Central America
central_america_countries_names = c("Belize","Guatemala","Honduras","El Salvador","Nicaragua","Costa Rica","Panama")
central_america_countries = countries.crop[countries.crop$NAME %in% central_america_countries_names,]
central_america_cores = intersect(cores.poly,central_america_countries)
central_america_cores$area_sqkm = area(central_america_cores) / 1000^2
total_ca_cores = sum(central_america_cores$area_sqkm) # total core are in Central America

ca_pa_intersect = intersect(pa_pts,central_america_cores)
ca_pa_summary = aggregate(ca_pa_intersect$protected,by=list(DN=ca_pa_intersect$DN), FUN=sum)
total_ca_cores_protected = sum(ca_pa_summary$x)
total_ca_unprotected = (total_ca_cores-total_ca_cores_protected)

total_ca_cores_protected/total_ca_cores # % of cores protected in Central America

# overall falls in indigenous area
iucn_indigenous_areas_proj$indigenous = 1
ia_raster = pa_raster
ia_raster = rasterize(iucn_indigenous_areas_proj,ia_raster,field="indigenous")
ia_pts = rasterToPoints(ia_raster,spatial=TRUE)
ia_intersect = intersect(ia_pts,central_america_cores)
ia_summary = aggregate(ia_intersect$layer,by=list(DN=ia_intersect$DN), FUN=sum)
sum(ia_summary$x)/total_ca_cores

# not protected, but falls in indigenous area
pa_raster_mask = pa_raster
pa_raster_mask[is.na(pa_raster_mask$protected)]=0
pa_raster_mask[pa_raster_mask==1]=NA
ia_raster_masked = mask(ia_raster,pa_raster_mask)
ia_masked_pts = rasterToPoints(ia_raster_masked,spatial=TRUE)
ia_masked_intersect = intersect(ia_masked_pts,central_america_cores)
ia_masked_summary = aggregate(ia_masked_intersect$layer,by=list(DN=ia_masked_intersect$DN), FUN=sum)
total_ia_not_pa = sum(ia_masked_summary$x)

total_ia_not_pa/total_ca_unprotected # % of unprotected core area that is in an indigenous area
(total_ia_not_pa+total_ca_cores_protected)/total_ca_cores # % of cores protected or indigenous in Central America 

# percent of core area that falls in a KBA
kbas_proj$kba = 1
kbas_raster = pa_raster
kbas_raster = rasterize(kbas_proj,kbas_raster,field="kba")
kbas_pts = rasterToPoints(kbas_raster,spatial=TRUE)
kbas_intersect = intersect(kbas_pts,cores.poly)
kbas_summary = aggregate(kbas_intersect$layer,by=list(DN=kbas_intersect$DN), FUN=sum)
sum(kbas_summary$x)/sum(cores.poly$area_sqkm)

# percent of core area that falls in a JCU
jcus_proj$jcu = 1
jcus_raster = pa_raster
jcus_raster = rasterize(jcus_proj,jcus_raster,field="jcu")
jcus_pts = rasterToPoints(jcus_raster,spatial=TRUE)
jcus_intersect = intersect(jcus_pts,cores.poly)
jcus_summary = aggregate(jcus_intersect$layer,by=list(DN=jcus_intersect$DN), FUN=sum)
sum(jcus_summary$x)/sum(cores.poly$area_sqkm)


### Discussion

length(unique(pas_cores_select@data$WDPAID)) # number of protected areas that intersect with cores
nrow(PA_FL_table) # number of those protected areas with greater than 5% forest loss
table(PA_FL_table$Location)

# how many patches do not have confirmed occupancy
sum(cores.poly$Presence == 0)
# what is their average population size
mean(cores.poly@data[cores.poly$Presence == 0, "sum"])
max(cores.poly@data[cores.poly$Presence == 0, "sum"])

cores.poly@data[cores.poly$DN == 30, "sum"] # population in La Amistad - Panama

####################################

# some stuff for the YIBS

hotspots = localG>1.96
hotspots.poly = gdal_polygonizeR(hotspots)
hotspots.poly = hotspots.poly[hotspots.poly@data$DN==1,]
coldspots = localG<(-1.96)
coldspots.poly = gdal_polygonizeR(coldspots)
coldspots.poly = coldspots.poly[coldspots.poly@data$DN==1,]

cellStats(intensity,"mean")

breaks = c(cellStats(intensity,"min"),0.01,0.02,0.04,0.08,0.16,cellStats(intensity,"max"))
breaks = round(breaks,2)
rng = breaks[c(1,4,5,6,7)]
arg = list(at=rng, labels=rng, cex.axis=0.75)
palette = brewer.pal(n = 6, name = "YlOrBr")

png(paste0(maps.folder,"intensity.png"), width=5, height=5, units="in", res=300) 
par(mar=c(0,0,0,0))
plot(countries.crop, border = FALSE)
e = extent(countries.crop)
rect(e[1], e[3], e[2], e[4], col = alpha("light blue", 0.3), border = FALSE)
plot(countries.crop, col = alpha("tan", 0.3), border = FALSE, add = TRUE)
plot(intensity, breaks=breaks, col=palette, axis.args=arg,  axes=FALSE, box=F, legend =F, add=T)
plot(lake, col=alpha("light blue", 0.6), border = FALSE, add = T)
plot(intensity, breaks=breaks, col=palette, axis.args=arg, axes=FALSE, box=F, legend.only=T, 
     smallplot=c(.05,.1, .15,.5))
dev.off()

breaks = c(cellStats(localG,"min"),-3,-2,-1,0,3,7,14,cellStats(localG,"max"))
rng = c(cellStats(localG,"min"),0,7,14,cellStats(localG,"max"))
arg = list(at=rng, labels=round(rng, 0), cex.axis=0.75)
palette = rev(brewer.pal(n = 8, name = "RdYlBu"))

png(paste0(maps.folder,"localG.png"), width=5, height=5, units="in", res=300) 
par(mar=c(0,0,0,0))
plot(countries.crop, border = FALSE)
e = extent(countries.crop)
rect(e[1], e[3], e[2], e[4], col = alpha("light blue", 0.3), border = FALSE)
plot(countries.crop, col = alpha("tan", 0.3), border = FALSE, add = TRUE)
plot(localG, breaks=breaks, col=palette, axis.args=arg,  axes=FALSE, box=F, legend =F, add=T)
plot(hotspots.poly,col=NA, border = 1,add=T)
plot(lake, col=alpha("light blue", 0.6), border = FALSE, add = T)
plot(localG, breaks=breaks, col=palette, axis.args=arg, axes=FALSE, box=F, legend.only=T, 
     smallplot=c(.05,.1, .15,.5))
dev.off()

# end stuff for YIBS

# Named patches

patch.names = c("Maya","Sian Ka'an","Oaxaca","Ocote","Chimalapas","Lacandona","Pico Bonito","Lacandona","9",
                "Karatasca","Montana de Botaderos","Maya","13","Santa B?rbara","Azul Me?mbar",
                "16","Monta?a Verde","Nicaragua","19","20","21","22","Talamanca-La Amistad to Golfo de los Mosquitos",
                "24","Osa","26","27","28")

patch.names.table = data.frame(cbind(cores.poly.countries$DN,patch.names))
names(patch.names.table) = c("Patch ID","Patch Names")
table.caption = "Patch names."
patch.names.table.xtable = xtable(patch.names.table, caption = table.caption)

sink(tablefile)
print(patch.names.table.xtable, type="latex", hline.after=c(0), include.rownames=F, comment=F, caption.placement = "bottom")

