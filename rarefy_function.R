
#spatial sampling procedure
rarefy = function(pointSample,r){
  pointSampleCopy = pointSample
  pointSampleTrain = pointSample[0,]
  pointSampleTest = pointSample[0,]
  nsamples=1
  while(dim(pointSampleCopy)[1] > 0){
    rowID = sample(nrow(pointSampleCopy), 1)
    s1 = pointSampleCopy[rowID, ]
    pointSampleTrain[nsamples,] = s1
    coordinates(s1) = c("X","Y")
    proj4string(s1) = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
    coordinates(pointSampleCopy) = c("X","Y")
    proj4string(pointSampleCopy) = "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
    pointSampleCopy$distance = (spDistsN1(pointSampleCopy,s1)>r)
    pointSampleCopy = data.frame(pointSampleCopy)
    pointSampleCopy = pointSampleCopy[-rowID,]
    pointSampleTemp = pointSampleCopy[!pointSampleCopy$distance,]
    #rownames(pointSampleTemp) = NULL
    #since I reinstalled R, for some reason it adds an extra column during this function, I get rid of it below
    pointSampleTemp$distance = NULL
    pointSampleTemp = pointSampleTemp[,-dim(pointSampleTemp)[2]]
    pointSampleTest = rbind(pointSampleTest,pointSampleTemp)
    pointSampleCopy  = pointSampleCopy[pointSampleCopy$distance,]
    pointSampleCopy$distance = NULL
    #since I reinstalled R, for some reason it adds an extra column during this function, I get rid of it below
    pointSampleCopy = pointSampleCopy[,-dim(pointSampleCopy)[2]]
    nsamples=nsamples+1
    #rownames(pointSampleTrain) = NULL
  }
  return(list(pointSampleTrain,pointSampleTest,r))
}
