### Exploratory Analysis of TOVE counts in lines and points ##

## Reading the data available ##


load("~/Documents/Research/DataIntegration/DeadBirds/Preprocessing.RData")
setwd("/Users/jorgespa/Documents/Research/DataIntegration/DeadBirds/")

library(raster)
names(norswe_lines_aspoints)[1] <- "ID"

## Creating the variable country ##
norswe_points_aspoints@data$country <- norswe_lines_aspoints@data$country <- NA ## Creating a variable called country

norswe_points_aspoints@data$country[!is.na(as.numeric(as.character(norswe_points_aspoints$ID)))] <-
  norswe_lines_aspoints@data$country[!is.na(as.numeric(as.character(norswe_lines_aspoints$ID)))] <-"Norway"

norswe_points_aspoints@data$country[is.na(as.numeric(as.character(norswe_points_aspoints$ID)))] <- norswe_lines_aspoints@data$country[is.na(as.numeric(as.character(norswe_lines_aspoints$ID)))] <-"Sweden"


## Creating the variable spprep ##
norswe_points_aspoints@data$spprep <- norswe_lines_aspoints@data$spprep <- NA

norswe_points_aspoints@data$spprep <- 219 ## This places have all species observed (the points)
norswe_lines_aspoints@data$spprep[which(norswe_lines_aspoints@data$country=="Norway")] <- 121 ## 121 species observed in lines in Norway
norswe_lines_aspoints@data$spprep[which(norswe_lines_aspoints@data$country=="Sweden")] <- 219 ## 219 species observed in lines in Norway

## Verifying counts for points and lines
boxplot(norswe_points_aspoints@data$Average[which(norswe_lines_aspoints@data$country=="Norway")],
        norswe_points_aspoints@data$Average[which(norswe_lines_aspoints@data$country=="Sweden")])

boxplot(norswe_lines_aspoints@data$Average[which(norswe_lines_aspoints@data$country=="Norway")],
        norswe_lines_aspoints@data$Average[which(norswe_lines_aspoints@data$country=="Sweden")])

norswe_lines_aspoints1 <- norswe_lines_aspoints[!is.na(norswe_lines_aspoints$Average) & !is.na(norswe_points_aspoints$Average),] ## Getting rid of NAs in averages in both points and counts
norswe_points_aspoints1 <- norswe_points_aspoints[!is.na(norswe_points_aspoints$Average) & !is.na(norswe_lines_aspoints$Average),] ## Getting rid of NAs in averages in both points and counts


## Putting all data together 

all_data <- cbind(norswe_lines_aspoints1@data[,c("ID","Average","country")],norswe_points_aspoints1@data[,"Average"])
names(all_data) <- c("ID","Average","country","Average2")
par(mfrow=c(1,2))
plot(all_data$Average[which(all_data$country=="Norway")],all_data$Average2[which(all_data$country=="Norway")],xlim=c(0,15),xlab = "Line Counts",ylab="Point Counts",main="Line vs Point Counts: Norway")

cor(all_data$Average[which(all_data$country=="Norway")],all_data$Average2[which(all_data$country=="Norway")])

plot(all_data$Average[which(all_data$country=="Sweden")],all_data$Average2[which(all_data$country=="Sweden")],xlab = "Line Counts",ylab="Point Counts",main="Line vs Point Counts: Sweden")
cor(all_data$Average[which(all_data$country=="Sweden")],all_data$Average2[which(all_data$country=="Sweden")])


all_data$ratio <- all_data$Average/all_data$Average2
all_data$ratio_inv <- all_data$Average2/all_data$Average
par(mfrow=c(1,2))
bp1.no <- boxplot(all_data$ratio[which(all_data$country=="Norway")],main="Line Count / Point Count : Norway")
bp1.se <-boxplot(all_data$ratio[which(all_data$country=="Sweden")],main="Line Count / Point Count : Sweden")

par(mfrow=c(1,2))
bp2.no <- boxplot(all_data$ratio_inv[which(all_data$country=="Norway")])
bp2.se <- boxplot(all_data$ratio_inv[which(all_data$country=="Sweden")])

## Getting the outliers in both ways in both countries ##

## Norway

no.ratio1 <- all_data$ratio[which(all_data$country=="Norway")]
no.out1 <- which(no.ratio1 %in% bp1.no$out)

no.ratio2 <- all_data$ratio_inv[which(all_data$country=="Norway")]
no.out2 <- which(no.ratio2 %in% bp2.no$out)

no.out <- c(no.out1,no.out2)
## Sweden

se.ratio1 <- all_data$ratio[which(all_data$country=="Sweden")]
se.out1 <- which(se.ratio1 %in% bp1.se$out)

se.ratio2 <- all_data$ratio_inv[which(all_data$country=="Sweden")]
se.out2 <- which(se.ratio2 %in% bp2.se$out)

se.out <- c(se.out1,se.out2)

par(mfrow=c(1,1))
plot(all_data$Average[which(all_data$country=="Norway")][-no.out],all_data$Average2[which(all_data$country=="Norway")][-no.out])
cor(all_data$Average[which(all_data$country=="Norway")][-no.out],all_data$Average2[which(all_data$country=="Norway")][-no.out])

boxplot(all_data$ratio[which(all_data$country=="Norway")][-no.out])
boxplot(all_data$ratio_inv[which(all_data$country=="Norway")][-no.out])


plot(all_data$Average[which(all_data$country=="Norway")],all_data$Average2[which(all_data$country=="Norway")],xlim=c(0,20))
points(all_data$Average[which(all_data$country=="Norway")][no.out],all_data$Average2[which(all_data$country=="Norway")][no.out],pch=19)

