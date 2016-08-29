####################################################
# TVC Radar Processing
# J. King ECCC
# 29/08/16
####################################################
rm(list = ls()) #Clear workspace

### Paths
showPlots = TRUE
workingPath = "C:/Users/kingj/Documents/Projects/2016-2017/010116_TVC_SnowSAR"
codePath = paste0(workingPath, "/Code/R/")
setwd(codePath)
setwd(workingPath)
sigma0Path = paste0(workingPath, "/Radar_data/GeoTIFF_sig0/")
elocalPath = paste0(workingPath, "/Radar_data/GeoTIFF_elocal/")
elangPath = paste0(workingPath, "/Radar_data/GeoTIFF_elang/")
outputPath = paste0(workingPath, "/Radar_data/Output/")
dir.create(outputPath, showWarnings = FALSE)

#require(sp)
require(rgdal)
#require(rgeos)
require(raster)
#require(ggplot2)
#require(parallel)
#require(gstat) 

Mode <- function(x,na.rm=FALSE) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

radarNativeRes = 2 #in m, this should not change
scaleData.radar = c(10/radarNativeRes, 50/radarNativeRes, 100/radarNativeRes)

cat("Processing TVC", mission, "...\n")
GCS = "+proj=longlat +ellps=WGS84" #CRS geographic coordniate system WGS84
PCS = "+proj=utm +zone=8 +north +units=m +ellps=WGS84" #CRS projected coordniate system UTM8N/WGS84

###Load Radar Data 
filePaths = list.files(sigma0Path, pattern="*.tif$", full.names=TRUE)
fileNames = list.files(sigma0Path, pattern="*.tif$", full.names=FALSE)
fileDates = unique(substring(fileNames ,9,22))
swathNum = length(fileDates)
cat("Radar swaths to be processed:",swathNum,"\n")
	
for (i in 1:swathNum) {
	print(paste0("Processing swath ", fileDates[i],"..."))
	radarData.s0Files = list.files(sigma0Path, pattern=paste0("*", fileDates[i],"*"), full.names=TRUE)
	radarData.elocalFiles = list.files(elocalPath , pattern=paste0("*", fileDates[i],"*"), full.names=TRUE)
	radarData.elangFiles = list.files(elangPath , pattern=paste0("*", fileDates[i],"*"), full.names=TRUE)
	
	#Load all 4 radar channels
	radarData.s0 <- mclapply(radarData.s0Files, raster)

	#Resample to a common grid
	radarData.s0 <- mclapply (1:length(radarData.s0), function(x) {resample(radarData.s0[[x]],radarData.s0[[1]], method="bilinear")})

	#Filtering (Here we apply a 3x3 low pass filter); uncomment the next two lines to use a gausian instead
	radarData.s0 <- mclapply(1:length(radarData.s0), function(x) {focal(radarData.s0[[x]], w=matrix(1/9,nrow=3,ncol=3))})
	#gf= focalWeight(radarData.s0[[1]], 3, "Gauss")
	#radarData.s0.filter <- lapply(1:length(radarData.s0), function(x) {focal(radarData.s0[[x]],w=gf)})
	
	#Open local radar geometry files
	radarData.elocal <- lapply(radarData.elocalFiles, raster)
	radarData.elocal <- mclapply(1:length(radarData.elocal), function(x) {reclassify(radarData.elocal[[x]], c(-360, 0, NA))}) #Remove weird angles
	radarData.elocal <- mclapply(1:length(radarData.elocal), function(x) {reclassify(radarData.elocal[[x]], c(70, 360, NA))}) #Remove weird angles
	radarData.elocal <- mclapply(1:length(radarData.elocal), function(x) {resample(radarData.elocal[[x]],radarData.s0[[x]], method="bilinear")})
		
	#Open local radar geometry files
	radarData.elang <- lapply(radarData.elangFiles, raster)
	radarData.elang <- mclapply(1:length(radarData.elang), function(x) {reclassify(radarData.elang[[x]], c(-360, 0, NA))}) #Remove weird angles
	radarData.elang <- mclapply(1:length(radarData.elang), function(x) {reclassify(radarData.elang[[x]], c(70, 360, NA))}) #Remove weird angles
	radarData.elang <- mclapply(1:length(radarData.elang), function(x) {resample(radarData.elang[[x]],radarData.s0[[x]], method="bilinear")})
	print(paste0("Swath ", fileDates[i]," filtered."))

	if(showPlots){
		range.KVH <- c(round(minValue(radarData.s0[[1]]),2), round(maxValue(radarData.s0[[1]]),2))
		range.KVV <- c(round(minValue(radarData.s0[[2]]),2), round(maxValue(radarData.s0[[2]]),2))
		range.XVH <- c(round(minValue(radarData.s0[[3]]),2), round(maxValue(radarData.s0[[3]]),2))
		range.XVV <- c(round(minValue(radarData.s0[[4]]),2), round(maxValue(radarData.s0[[4]]),2))

		graphics.off()
		par(mfrow=c(2,2))

		plot(radarData.s0[[1]], maxpixels=50000, col = grey.colors(255, alpha = 1),
			main=paste("Ku VH", fileDates[i]),
			axis.args=list(at=seq(range.KVH[1], range.KVH[2], 0.02),
            labels=seq( range.KVH[1],  range.KVH[2], 0.02), cex.axis=0.6),
			legend.args=list(text=expression(sigma^0), side=4, font=2, line=2.5, cex=0.8))

		plot(radarData.s0[[2]], maxpixels=50000, col = grey.colors(255, alpha = 1),
			main=paste("Ku VV", fileDates[i]),
			axis.args=list(at=seq( range.KVV[1], range.KVV[2], 0.02),
            labels=seq( range.KVV[1],  range.KVV[2], 0.02), cex.axis=0.6),
			legend.args=list(text=expression(sigma^0), side=4, font=2, line=2.5, cex=0.8))

		plot(radarData.s0[[3]], maxpixels=50000, col = grey.colors(255, alpha = 1),
			main=paste("X VH", fileDates[i]),
			axis.args=list(at=seq( range.XVH[1],  range.XVH[2], 0.02),
            labels=seq( range.XVH[1],  range.XVH[2], 0.02), cex.axis=0.6),
			legend.args=list(text=expression(sigma^0), side=4, font=2, line=2.5, cex=0.8))

		plot(radarData.s0[[4]], maxpixels=50000, col = grey.colors(255, alpha = 1),
			main=paste("X VV", fileDates[i]),
			axis.args=list(at=seq(range.XVV[1], range.XVV[2], 0.02),
            labels=seq(range.XVV[1], range.XVV[2], 0.02), cex.axis=0.6),
			legend.args=list(text=expression(sigma^0), side=4, font=2, line=2.5, cex=0.8))
		}#END PLOTS
	
	print(paste0("Scaling ", fileDates[i]," ...."))
	#Scale data to desired resolution(s); uses the mean as a default, but you can add whatever function is desired
	#I only aggregate the geometry for one fo the bands (KVV) because it should all be equal
	sigFileNames = list.files(sigma0Path, pattern=paste0("*", fileDates[i],"*"), full.names=FALSE)
	elangFileNames = list.files(elangPath, pattern=paste0("*", fileDates[i],"*"), full.names=FALSE)
	elocalFileNames = list.files(elangPath, pattern=paste0("*", fileDates[i],"*"), full.names=FALSE)

	if(scaleData.radar[1]>0){
		for (s in 1:length(scaleData.radar)){
			sigOutFile = paste0(outputPath,substring(sigFileNames ,1,32),scaleData.radar[s]*radarNativeRes,"m.tif")
			elangOutFile = paste0(outputPath,substring(elangFileNames[[1]] ,1,23),"elang_",scaleData.radar[s]*radarNativeRes,"m.tif")
			elocalOutFile = paste0(outputPath,substring(elocalFileNames[[1]] ,1,23),"elocal_",scaleData.radar[s]*radarNativeRes,"m.tif")
			radarData.s0.agg <- lapply(1:length(radarData.s0), function(x) {aggregate(radarData.s0[[x]], scaleData.radar[s], fun=mean, na.rm=TRUE, expand=FALSE)})
			radarData.elocal.agg <- aggregate(radarData.elocal[[2]], scaleData.radar[s], fun=mean, na.rm=TRUE, expand=FALSE)
			radarData.elang.agg <- aggregate(radarData.elang[[2]], scaleData.radar[s], fun=mean, na.rm=TRUE, expand=FALSE)
			lapply(1:length(radarData.s0.agg), function(x) {writeRaster(radarData.s0.agg[[x]], filename=sigOutFile[x], format="GTiff", overwrite=TRUE)})
			writeRaster(radarData.elocal.agg, filename=elocalOutFile, format="GTiff", overwrite=TRUE)
			writeRaster(radarData.elang.agg, filename=elangOutFile, format="GTiff", overwrite=TRUE)
			}
		}
	}
	print(paste0("Finished swath ", fileDates[i], format="GTiff", overwrite=TRUE)
}