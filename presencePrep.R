library(sp)
library(plyr)
library(rgdal)

# projection to use
sr_geo = " +init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
sr_cea = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

outputFolder = "/Users/Cody/Dropbox/Research/Dissertation/Presence_Data/"

# Prep PO data
inputFolder = "/Users/Cody/Dropbox/Research/Dissertation/Presence_Data/presence-only/"
poFiles = list.files(inputFolder, pattern = ".csv", full.names = TRUE)
merged_po = do.call("rbind", lapply(poFiles, read.csv, header = TRUE, stringsAsFactors=FALSE))
coordinates(merged_po) = c("Longitude","Latitude")
proj4string(merged_po) = sr_geo
merged_po_proj = spTransform(merged_po,sr_cea)
po_final = data.frame(cbind(coordinates(merged_po),coordinates(merged_po_proj)))
names(po_final) = c("Longitude","Latitude","X","Y")
po_final$species = "Tapirus bairdii"
#remove duplicates, if they are not already removed
write.csv(unique(po_final),paste0(outputFolder,"presence_only_proj.csv"), row.names=FALSE)
writeOGR(merged_po_proj,dsn=paste0(outputFolder,"presence_only_proj.shp"),layer="presence_only_proj",driver="ESRI Shapefile", overwrite_layer=T)

# Prep PA data
inputFolder = "/Users/Cody/Dropbox/Research/Dissertation/Presence_Data/detection_histories"
siteFiles = list.files(inputFolder, pattern = "sites.csv", full.names = TRUE)

# check input for mismatches
for(site in siteFiles){
  
  stations = read.csv(site,stringsAsFactors = FALSE)
  detections = read.csv(gsub("sites","detections",site), stringsAsFactors = FALSE)
  
  test_join = merge(detections,stations,by="Site",all.x=TRUE)
  
  if(!("start_date" %in% names(test_join))){print(paste0("start_name not in ",site)); break}
  
  missing_stations = unique(test_join[is.na(test_join$start_date),"Site"])
  
  if(length(missing_stations) > 0){print(paste0(paste(missing_stations, collapse = ", "), " not in ",site)); break }

}

# combine files into data.frames
merged_sites = do.call("rbind", lapply(siteFiles, read.csv, header = TRUE, stringsAsFactors=FALSE))
merged_sites$start_date = as.Date(merged_sites$start_date,format="%m/%d/%Y")
merged_sites$end_date = as.Date(merged_sites$end_date,format="%m/%d/%Y")

# this should get rid of NA rows
merged_sites = merged_sites[!is.na(merged_sites$Latitude),]

# check if there are duplicated names in sites
sum(duplicated(merged_sites$Site)) > 0

detectionFiles = list.files(inputFolder, pattern = "detections.csv", full.names = TRUE)
merged_detections = do.call("rbind", lapply(detectionFiles, read.csv, header = TRUE, stringsAsFactors=FALSE))

max_samples = 10
sample_occasion = 10

det_history_final = data.frame(matrix(NA, nrow=nrow(merged_sites), ncol=max_samples))
names(det_history_final) = paste0("sample_", seq(1,max_samples))

det_history_final = cbind(merged_sites,det_history_final)

det_history_final$Days = NA

test = 0

for(i in seq(1,nrow(merged_sites))){
  
  site = merged_sites[i,"Site"]
  start_date = merged_sites[i,"start_date"]
  end_date = merged_sites[i,"end_date"]
  trap_days = as.numeric(end_date - start_date)
  
  det_history_final[i,"Days"] = trap_days
  det_history = rep(0,trap_days)
  
  missing_dates = merged_sites[i,"missing_dates"]
  if(!is.na(missing_dates)){
    missing_dates_split = as.Date(unlist(strsplit(missing_dates,",")),format="%m/%d/%Y")
    days_sequence = seq(start_date, end_date, by="days")
    missing_dates_index = which(days_sequence %in% missing_dates_split)
    det_history[missing_dates_index] = NA
  }
  
  # this adds necessary NAs to the data
  if(length(det_history) < max_samples*sample_occasion){
    NAs = rep(NA,((max_samples*sample_occasion)-length(det_history)))
    det_history = c(det_history,NAs)
  }
  
  detections = merged_detections[merged_detections$Site == site,]
  
  if(nrow(detections) > 0){
    for(j in seq(1,nrow(detections))){
      capture_day = as.numeric(as.Date(detections[j,"Detection"],format="%m/%d/%Y") - start_date)
      det_history[(capture_day+1)] = 1
      #test for problems, usually with dates
      if((capture_day < 0)|(abs(capture_day) > trap_days)){test=1}
    }
  }

  #test for problems, usually with dates
  if(test==1){print(paste0("test for problems, usually with dates")); break}
  
  det_history_reduced = rep(NA,max_samples)
  
  breaks = c(0,seq(sample_occasion,max_samples*sample_occasion,sample_occasion))
  for (j in seq(1,(max_samples))){
      det_history_reduced[j] <- (sum(det_history[(breaks[j]+1):breaks[j+1]]) > 0) * 1
  }
  
  names(det_history_reduced) =  paste0("sample_", seq(1,max_samples))
  det_history_final[i,paste0("sample_", seq(1,max_samples))] = det_history_reduced
}

det_history_final = det_history_final[det_history_final$Days >= 10,]

coordinates(det_history_final) = c("Longitude","Latitude")
proj4string(det_history_final) = sr_geo

det_history_final_proj = spTransform(det_history_final,sr_cea)

det_history_final_coords = data.frame(coordinates(det_history_final))
det_history_final_proj_coords = data.frame(coordinates(det_history_final_proj))
names(det_history_final_proj_coords) = c("X","Y")

det_history_final_join = cbind(det_history_final_coords,det_history_final_proj_coords,det_history_final_proj@data)

det_history_final_join$Year = as.numeric(format(det_history_final_join$start_date,'%Y'))
det_history_final_join$Species = "Tapirus bairdii"

#det_history_final_join = det_history_final_join[,c("Site","Provider","Year","Days","Longitude",
#                                                   "Latitude",paste0("sample_", seq(1,max_samples)),
#                                                   "Species","X","Y")]

det_history_final_join = det_history_final_join[,c("Provider","Year","Days","Longitude",
                                                   "Latitude",paste0("sample_", seq(1,max_samples)),
                                                   "Species","X","Y","Tapir","Cat","OnTrail")]

det_history_final_join$Presence = (rowSums(det_history_final_join[,paste0("sample_", seq(1,max_samples))],na.rm = T)>0)*1


write.csv(det_history_final_join,paste0(outputFolder,"detection_histories_proj_",max_samples,"samples.csv"), row.names=FALSE)

det_history_final_join_shp = det_history_final_join
coordinates(det_history_final_join_shp) = c("X","Y")
proj4string(det_history_final_join_shp) = sr_cea
writeOGR(det_history_final_join_shp,dsn=paste0(outputFolder,"detection_histories_proj.shp"),layer="detection_histories_proj",driver="ESRI Shapefile", overwrite_layer=T)

### Appendix Table

inputFolder = "/Users/Cody/Dropbox/Research/Dissertation/Presence_Data/presence-only/"
poFiles = list.files(inputFolder, pattern = ".csv", full.names = TRUE)

inputFolder = "/Users/Cody/Dropbox/Research/Dissertation/Presence_Data/detection_histories"
siteFiles = list.files(inputFolder, pattern = "sites.csv", full.names = TRUE)

allFiles = c(poFiles,siteFiles)
appendixTable = data.frame(allFiles)
names(appendixTable) = "Filename"

inputFolder = "/Users/Cody/Dropbox/Research/Dissertation/dataPrep_input/"
countries = readOGR(dsn=paste0(inputFolder,"countries.crop.shp"),layer="countries.crop")
countries.geog = spTransform(countries,sr_geo)

for(dataFile in allFiles){
  
  inData = read.csv(dataFile,header=TRUE, stringsAsFactors = FALSE)
  numLocations = nrow(inData)
  coordinates(inData) = c("Longitude","Latitude")
  proj4string(inData) = sr_geo
  
  country = unique(na.omit(over(inData,countries.geog)["ISO_CODE"]))[1]
  
  if(grepl("presence-only",dataFile)){
    type = "PO"
    provider = gsub(".csv","",basename(dataFile))
  
  }else{
    provider = inData$Provider[1]
    type = "PA"
  }
  
  appendixTable[appendixTable$Filename == dataFile,"Type"] = type
  appendixTable[appendixTable$Filename == dataFile,"Provider"] = provider
  appendixTable[appendixTable$Filename == dataFile,"Locations"] = numLocations
  appendixTable[appendixTable$Filename == dataFile,"Country"] = country

}

write.csv(appendixTable, paste0(outputFolder,"data_sources.csv"),row.names=FALSE)

#######################################
##### Not sure what this stuff is #####
#######################################

#temp1 = apply(det_history_final[,paste0("sample_", seq(1,max_samples))] > 0 , MARGIN = 1, FUN=which)
#x1 = suppressWarnings(unlist(lapply(temp1, min)))
#x1 = x1[x1<500]
#plot(density(x1))
#x2 = merge(x1,det_history_final[,c("Site","Days")],by="row.names")
x2 = x2[,c("x","Days")]
x2 = x2[x2$Days>100,]
#day of first detection is conditional on number of days in detection history, 
#i.e. more days allows for the possibility of a later first detection

# 860 detections in whole dataset, at 337 sites...out of 1163 sites
# 62 sites out of 200 in the psi test data have detections (if I use 8 sampling occasions)
# 28 sites out of 93 in the psi test data hasve detections (if I use 10 sampling occasions)
# 11 of 200 sites have a detection after 80 days, but not before
# at some point the assumption of closure could be violated, and that could cause these sites to have detections after x days

#from model18 with 5x10 detection histories, alpha.pa = -0.8, thus p = 0.3100255
#from model18 with 8x10 detection histories, alpha.pa = -0.95, thus p = 0.28
#dbinom(0, 5, 0.28) = 0.05, it is present and not detected
#max days = 484

#sum(is.na(det_history_final_join[,6:13]))/(nrow(det_history_final)*8) #40% NA ratio
#sum(is.na(det_history_final_join[,6:10]))/(nrow(det_history_final)*5) #20% NA ratio