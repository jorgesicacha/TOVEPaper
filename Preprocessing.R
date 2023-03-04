###### Data preprocessing #####

library(rgdal)
library(raster)
library(rgeos)
library(dismo)
library(stringr)
library(RColorBrewer)
library(mapview)
library(viridis)
library(geosphere)
library(stringr)
library(sjmisc)
library(remotes)
library(threadr)
library(gdalUtils)
library(sqldf)
library(gstat)
library(plyr)
library(reshape2)
library(maptools)
library(ggplot2)


## Species list ##

specieslist <- read.csv("All_spp_lines_and_points.csv",header=T,sep=";",dec=".")
speciesinters <- specieslist[which(specieslist$Is_in_point==1  & specieslist$Is_in_transect==1),]
speciesonlypt <- specieslist[which(specieslist$Is_in_point==1  & is.na(specieslist$Is_in_transect)),]
speciesonlyline <- specieslist[which(is.na(specieslist$Is_in_point)  & specieslist$Is_in_transect==1),]

## Data from Norway ##

## Point counts ##
pointcount <- read.csv("Point_counts_converted_flock_to_pairs_SUMMED_BY_SPECIES.csv",header=T,sep=";",dec=".")
names(pointcount)[c(1,5)] <- c("ID","Yearx")
pointcount$Total <- pointcount$TOTAL_COUNT 

checkdata0 <- pointcount[!is.na(pointcount$Yearx),] ##Removing year=NA
checkdata0 <- checkdata0[which(checkdata0$Yearx!=2020),] ##Removing year=2020
checkdata <- dcast(checkdata0,ID+EURINGCode~Yearx,value.var = "Total",fun.aggregate = sum,fill=NA_real_) ##Creating the table as we want it!
checkdata$ID <- as.factor(as.character(checkdata$ID))
# ## Changing NA's for 0's where there has been sampling ##

for(i in 1:length(unique(checkdata$ID))){
  for(j in 1:(length(names(checkdata))-2)){
    if(nrow(pointcount[which(as.character(pointcount$ID)==as.character(unique(checkdata$ID)[i]) & pointcount$Yearx==2005+j),])!=0){checkdata[which(checkdata$ID==unique(checkdata$ID)[i]),j+2][is.na(checkdata[which(checkdata$ID==unique(checkdata$ID)[i]),j+2])]=0}
   }
}

## Aggregating temporally and by polygon ##

checkdata$LineID <- sub("_.*", "", checkdata$ID)

checkdataf <- checkdata
checkdata1 <- stats::aggregate.data.frame(checkdataf[,-c(1,2,ncol(checkdataf))],by=list(checkdataf$LineID),FUN=sum)
names(checkdata1)[1] <- "ID"

summat <- do.call("rbind",apply(checkdata1[,-1],1,summary))
total <- rowSums(checkdata1[,-1],na.rm = T)
num.pers <- apply(checkdata1,1,function(x) length(x[-1][!is.na(x[-1])]))
which.keep <-which(num.pers>0)

pointcount_sum <- data.frame(ID=checkdata1$ID[which.keep],Total=total[which.keep],Years=num.pers[which.keep],
                             Average = round(summat[which.keep,4],0))

pointcount_sum_final <- pointcount_sum

## Line counts ##

transcount <- read.csv("Line_counts_converted_flock_to_pairs.csv",header=T,sep=";",dec=".")
transcountinters <- transcount[which(transcount$EURINGCode %in% speciesinters$EURINGCode),] ## Counts for species observed in lines and points
others <- transcount[-which(transcount$EURINGCode %in% speciesinters$EURINGCode),] ## Counts for species only observed in lines

names(transcount)[c(4,5)] <- c("ID","Yearx")
transcount$Total <- transcount$Total_line_count 

checkdata0 <- transcount[!is.na(transcount$Yearx),] ##Removing year=NA
checkdata0 <- checkdata0[which(checkdata0$Yearx!=2020 & checkdata0$Yearx!=2005),] ##Removing year=2020 & 2005

checkdata_line <- dcast(checkdata0,ID+EURINGCode~Yearx,value.var = "Total",fun.aggregate = sum,fill=NA_real_)
checkdata_line$ID <- as.factor(as.character(checkdata_line$ID))

# ## Changing NA's for 0's where there has been sampling ##

for(i in 1:length(unique(checkdata_line$ID))){
  for(j in 1:(length(names(checkdata_line))-2)){
    if(nrow(transcount[which(as.character(transcount$ID)==as.character(unique(checkdata_line$ID)[i]) & transcount$Yearx==2005+j),])!=0){checkdata_line[which(checkdata_line$ID==unique(checkdata_line$ID)[i]),j+2][is.na(checkdata_line[which(checkdata_line$ID==unique(checkdata_line$ID)[i]),j+2])]=0}
  }}

checkdataf <- checkdata_line
checkdata1_line <- stats::aggregate.data.frame(checkdataf[,-c(1,2)],by=list(checkdataf$ID),FUN=sum)
names(checkdata1_line)[1] <- "ID"

summat <- do.call("rbind",apply(checkdata1_line[,-1],1,summary))
total <- rowSums(checkdata1_line[,-1],na.rm = T)
num.pers <- apply(checkdata1_line,1,function(x) length(x[-1][!is.na(x[-1])]))
which.keep <-which(num.pers>0)

linecount_sum <- data.frame(ID=checkdata1_line$ID[which.keep],Total=total[which.keep],Years=num.pers[which.keep],
                            Average = round(summat[which.keep,4],0))

## Putting all Norwegian data together ##
birdcensus <- shapefile("TOVERoutes_3006.shp")

birdcensus$x <- birdcensus@coords[,1]
birdcensus$y <- birdcensus@coords[,2]

listtmp <- lapply(unique(birdcensus$RuteID), function(x){
  Polygons(list(Polygon(birdcensus[which(as.character(birdcensus$RuteID)==x),c("x","y")])),ID=as.character(x))
})

spbirdspoly <- lapply(listtmp, function(x){
  a <- gBuffer(SpatialPolygons(list(x),proj4string = crs(crs_new)))
  a@polygons[[1]]@ID=x@ID
  b <- a@polygons[[1]]
  return(b)
})
spbirdspoly.buffer <- SpatialPolygons(Srl=spbirdspoly,proj4string = crs(crs_new))
mapview(spbirdspoly.buffer,viewer.suppress=T)

df0 <- data.frame(ID=as.character(unique(birdcensus$RuteID)))
linecount_poly_df <- plyr::join(df0,linecount_sum,"ID",type="left")
pointandlinecount_poly_df <- plyr::join(linecount_poly_df,pointcount_sum,"ID",type="left")

row.names(pointandlinecount_poly_df) <- as.character(unique(birdcensus$RuteID))
spbirdspoly.df <- SpatialPolygonsDataFrame(spbirdspoly.buffer,data = pointandlinecount_poly_df)

dmesh <- shapefile("dmeshgood.shp")
dmesh <- spTransform(dmesh,crs(crs_new))
pointandlinecount_sum_coord.sp <- raster::intersect(spbirdspoly.df,dmesh)

## Data from Sweden ##

f1 <- shapefile("stdruttPunkter_3006.shp")
f2 <- shapefile("stdruttPol_3006.shp")
f3 <- shapefile("stdruttLinjer_3006.shp")

## Point counts

swecounts <- read.csv("Bird_data_Sweden_2016_onwards.csv",header=T,sep=";",dec=".")
speciesid <- read.csv("speciescodes.csv",header=T,sep=";",dec=".")
speciesid <- speciesid[,c("art","euring")]

swecounts <- plyr::join(swecounts,speciesid,by="art",type="left")

##  We need to summarize this!! ##
swepointcounts <- swecounts[,c(1:13,22:24)] #Keep only information related to points
swepointcounts <- swepointcounts[which(swepointcounts$art!=0 & swepointcounts$art!=999),] ##Removing non-informative rows
swepointcounts <- swepointcounts[,c(2,4,6:13,16)]
swepointcounts <- swepointcounts[!is.na(swepointcounts$p1) | !is.na(swepointcounts$p2) |!is.na(swepointcounts$p3) | !is.na(swepointcounts$p4)|!is.na(swepointcounts$p5) | !is.na(swepointcounts$p6)|!is.na(swepointcounts$p7) | !is.na(swepointcounts$p8),] #Removing rows without observations at any of the eight points

library(reshape2)
swepoints0 <- melt(swepointcounts,id=c("karta","yr","euring"))
swepoints0$ID <- paste(swepoints0$karta,swepoints0$variable,sep="_")
swepointscheck <-  dcast(swepoints0,ID+euring~yr,value.var = "value",fun.aggregate = sum,fill=NA_real_)

swepointscheck$ID <- as.factor(as.character(swepointscheck$ID))
# ## Changing NA's for 0's where there has been sampling ##
for(i in 1:length(unique(swepointscheck$ID))){
  for(j in 1:(length(names(swepointscheck))-2)){
    if(nrow(swepointcounts[which(as.character(swepoints0$ID)==as.character(unique(swepointscheck$ID)[i]) & swepoints0$yr==2005+j),])!=0){swepointscheck[which(swepointscheck$ID==unique(swepointscheck$ID)[i]),j+2][is.na(swepointscheck[which(swepointscheck$ID==unique(swepointscheck$ID)[i]),j+2])]=0}
  }}

## Aggregating temporally and by polygon ##

swepointscheck$LineID <- substr(swepointscheck$ID,1,5) 
checkdataf <- swepointscheck
checkdata1_point_swe <- stats::aggregate.data.frame(checkdataf[,-c(1,2,ncol(checkdataf))],by=list(checkdataf$LineID),FUN=sum)
names(checkdata1_point_swe)[1] <- "ID"

summat <- apply(checkdata1_point_swe[,-1],1,summary)
total <- rowSums(checkdata1_point_swe[,-1],na.rm = T)
num.pers <- apply(checkdata1_point_swe,1,function(x) length(x[-1][!is.na(x[-1])]))
which.keep <-which(num.pers>0)

pointcount_sum_swe <- data.frame(ID=checkdata1_point_swe$ID[which.keep],Total=total[which.keep],Years=num.pers[which.keep],
                                 Average = round(do.call("c",lapply(summat[which.keep],function(x) x[4])),4))

pointcount_sum_swe_final <- pointcount_sum_swe

## Line counts

swelinecounts <- swecounts[,c(1:5,14:24)]
swelinecounts <- swelinecounts[which(swelinecounts$art!=0 & swelinecounts$art!=999),]
swelinecounts <- swelinecounts[,c(2,4,6:13,16)]
swelinecounts <- swelinecounts[!is.na(swelinecounts$l1) | !is.na(swelinecounts$l2) |                                                  !is.na(swelinecounts$l3) | !is.na(swelinecounts$l4)|                                                   !is.na(swelinecounts$l5) | !is.na(swelinecounts$l6)|                                                   !is.na(swelinecounts$l7) | !is.na(swelinecounts$l8),]

swelines0 <- melt(swelinecounts,id=c("karta","yr","euring"))
swelines0$ID <- paste(swelines0$karta,swelines0$variable,sep="_")
swelinescheck <-  dcast(swelines0,ID+euring~yr,value.var = "value",fun.aggregate = sum,fill=NA_real_)

swelinescheck$ID <- as.factor(as.character(swelinescheck$ID))
# ## Changing NA's for 0's where there has been sampling ##
for(i in 1:length(unique(swelinescheck$ID))){
  for(j in 1:(length(names(swelinescheck))-2)){
    if(nrow(swelinecounts[which(as.character(swelines0$ID)==as.character(unique(swelinescheck$ID)[i]) & swelines0$yr==2005+j),])!=0){swelinescheck[which(swelinescheck$ID==unique(swelinescheck$ID)[i]),j+2][is.na(swelinescheck[which(swelinescheck$ID==unique(swelinescheck$ID)[i]),j+2])]=0}
  }}

swelinescheck$LineID <- substr(swelinescheck$ID,1,5) 
checkdataf <- swelinescheck
checkdata1_line_swe <- stats::aggregate.data.frame(checkdataf[,-c(1,2,ncol(checkdataf))],by=list(checkdataf$LineID),FUN=sum)
names(checkdata1_line_swe)[1] <- "ID"

summat <- apply(checkdata1_line_swe[,-1],1,summary)
total <- rowSums(checkdata1_line_swe[,-1],na.rm = T)
num.pers <- apply(checkdata1_line_swe,1,function(x) length(x[-1][!is.na(x[-1])]))
which.keep <-which(num.pers>0)

linecount_sum_swe <- data.frame(ID=checkdata1_line_swe$ID[which.keep],Total=total[which.keep],Years=num.pers[which.keep],
                                Average = round(do.call("c",lapply(summat[which.keep],function(x) x[4])),4))

f2_data <- f2@data
names(f2_data)[1] <- "ID"

test <- plyr::join(f2_data,pointcount_sum_swe_final,by="ID",type="left")
test <- plyr::join(test,linecount_sum_swe,by="ID",type="left")
f2@data <- test
swepointandlinecount <- f2

pointandlinecount_sum_swe_coord.sp <- raster::intersect(swepointandlinecount,dmesh)

## Making points object ##
pointcount_sum_swe_coord.sp <- pointandlinecount_sum_swe_coord.sp[,c(1,6,7,8)]
names(pointcount_sum_swe_coord.sp) <- c("ID","Total","Years","Average")
pointcount_sum_coord.sp <- pointandlinecount_sum_coord.sp[,c(1,5,6,7)]
names(pointcount_sum_coord.sp) <- c("ID","Total","Years","Average")

linecount_sum_swe_coord.sp <- pointandlinecount_sum_swe_coord.sp[,c(1,9,10,11)]
names(linecount_sum_swe_coord.sp) <- c("ID","Total","Years","Average")
linecount_sum_coord.sp <- pointandlinecount_sum_coord.sp[,c(1,2,3,4)]
names(linecount_sum_coord.sp) <- c("ID","Total","Years","Average")

norswe_points <- rbind(pointcount_sum_coord.sp,pointcount_sum_swe_coord.sp)
norswe_lines <- rbind(linecount_sum_coord.sp,linecount_sum_swe_coord.sp)

centroids <- do.call("rbind",lapply(norswe_lines@polygons, function(x){colMeans(x@Polygons[[1]]@coords)}))

norswe_lines_aspoints <- SpatialPointsDataFrame(centroids,data= norswe_lines@data,proj4string = crs(crs_new))
norswe_points_aspoints <- SpatialPointsDataFrame(centroids,data= norswe_points@data,proj4string = crs(crs_new))
