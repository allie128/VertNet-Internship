## For outlier flagging based on climate extracts from range polygon

##Produces original dataset with flags if the point is too far from the range (distanceflag3sd), in a
##zoo (zoo flag), not located on land (landflag), or has very different climatic characteristics 
##(climflag)

##also produces dataset that excludes flagged pts

getwd()
setwd()

install.packages(c('XML','dismo','maptools', 'rgdal','sp', 'raster','rgeos','plyr','rvertnet','spatstat')); #if you have not already installed these packages

library(dismo)
library(rgdal)
library(maptools)
library(sp)
library(rgeos)
library(XML)
library(plyr)
library(raster)
library(spatstat)


#Input species info first
Genus <- "Pseudacris"
species <- "crucifer"

#Read in range map by choice - make sure it is projected to Geographic
Rangesp1 <- readShapeSpatial(file.choose())


##FROM HERE YOU CAN RUN THE REST OF THE SCRIPT WITHOUT INPUT##


##if range is in two separate shapefiles you must combine them  into one before running program
##if this step fails you may have to change the second entry, which specifies which field 
##to use to join the polygons. this just has to be a field where the values for all the polygons are the same
Rangesp1 <- unionSpatialPolygons(Rangesp1, Rangesp1$BINOMIAL, threshold=5) #Joins if more than one polygon

#Acquire species occurrence records from gbif 
species1<-gbif(Genus, species, geo=T); #returns records with Georefs only
species1 <-subset(species1,lat !=0 & lon !=0) #makes sure there are no zero values for coordinates
species1_xy <- species1
#Set coordinates to a Spatial data frame
coordinates(species1_xy) = c("lon", "lat")

#attach calculation to species1
species1["distfromrange"] <- gDistance(Rangesp1, species1_xy, byid=TRUE)
species1['index'] <- 1:nrow(species1)

#subset points outside and inside range
ptsoutsiderange <- subset(species1, distfromrange > 0)
ptsinsiderange <- subset(species1, distfromrange == 0)

##DISTANCE FLAG##

#flags pts that are more than 3sd (using dist from range for pts outside range) away from known range

species1['distanceflag3sd'] <- ifelse( 
  species1['distfromrange'] > (3*sd(ptsoutsiderange$distfromrange)) , FALSE, TRUE
)
#converts logical vector to numeric vector where true=0
species1['distanceflag3sd'] <- ifelse(species1$distanceflag3sd, 0, 1)

##LAND FLAG##

countries <- getData('countries')
o = overlay(species1_xy,countries)
waterpts <- which(is.na(o)) #index of pts that are not on land
species1['landflag'] <- 0
species1[waterpts,'landflag'] <- 1

##ZOO FLAG##

wordcheck <- as.character(species1$locality)
wordvector <- strsplit(wordcheck, split=" ")
zoomatches1 <- grep('Zoo', wordvector)
zoomatches <- grep('zoo', wordvector)
ALLzoomatches <- c(zoomatches1,zoomatches)
species1['zooflag'] <- 0
species1[ALLzoomatches,'zooflag'] <- 1

##Climate flag##

#Get world clim data
WorldClimStack <- getData('worldclim', var='bio', res=10)
#extracts all the clim data from inside the range
rangeextract <- extract(WorldClimStack, Rangesp1, na.rm=TRUE)
rextract <- as.data.frame(rangeextract)

#test for normality
clip <- rextract[sample(nrow(rextract), 5000), ]
shapiro <- matrix(nrow=19, ncol=2)
shapiro[,2] <- 1:19
for (i in 1:19) {
  shapiro[i,1] <- shapiro.test(clip[,i])$p.value
}
shapiro <- as.data.frame(shapiro)
shapiro1 <- order(shapiro[,1], shapiro[,2])
indexnorm <- matrix(ncol=1, nrow=10)
indexnorm <- shapiro1[10:19] #picks only to 10 'most normal' distributions to compare. 
#extracts from only 'normal' variables
normextract <- matrix(nrow=nrow(rextract), ncol=10)
normextract <- rextract[,indexnorm]
#means and sd of normal variables
rangemeans <- apply(normextract, 2, mean, na.rm=TRUE)
rangesd <- apply(normextract, 2, sd, na.rm=TRUE)
#extract values for all pts
sp1ALL_climateextracts<-extract(WorldClimStack, species1[, c('lon','lat')])
sp1_norm <- sp1ALL_climateextracts[,indexnorm]

sp1_scores <- matrix(ncol=10, nrow=(nrow(sp1_norm)))

for (i in 1:10){
  sp1_scores[,i] <- (sp1_norm[,i] - rangemeans[i]) / rangesd[i]
}

abs_sp1_scores <- abs(sp1_scores)
abs_sp1_scores <- as.data.frame(abs_sp1_scores)
abs_sp1_scores['scoresum'] <- apply(abs_sp1_scores, 1, sum, na.rm=TRUE)
abs_sp1_scores['climflag'] <- ifelse(abs_sp1_scores$scoresum > mean(abs_sp1_scores$scoresum) + (3*sd(abs_sp1_scores$scoresum)), 1,0)

species1allflags <- cbind(species1, abs_sp1_scores['climflag'])

#unflag points that are within a certain distance of good points

climflagpts <- subset(species1allflags, !climflag==0)
climnoflagpts <- subset(species1allflags, !climflag==1)

x<-matrix(nrow=nrow(climflagpts), ncol=2)
flagfix <- vector(length=nrow(x))
x[,1] <- climflagpts$lon
x[,2] <- climflagpts$lat
##search for unflagged pts close to flagged pts to signal flag may be false
for (h in 1:nrow(x)) {
  point <- matrix(ncol=2, nrow=1)
  point[1,1] <- x[h,1]
  point[1,2] <- x[h,2]
  window <- disc(radius=1, centre=point)
  pointasppp <- ppp(climnoflagpts$lon, climnoflagpts$lat, window=window, marks=climnoflagpts$index)
  as.ppp(pointasppp)
  closepts <- as.data.frame(pointasppp)
  flagfix[h] <- ifelse(nrow(closepts)>0, TRUE, FALSE)
}

##TRUE = fix the clim flag
climflagpts <- cbind(climflagpts, flagfix)
subindex <- ifelse(climflagpts$flagfix, climflagpts$index, 0)
subindex <- subset(subindex, subindex>0)
species1allflags['flagfix'] <- 0
species1allflags[subindex,'flagfix'] <- 1

species1allflags['climflag2'] <- species1allflags$climflag + species1allflags$flagfix

species1allflags['climflag2'] <- ifelse(species1allflags$climflag2 == 1, 1, 0)


#subset points to only include those without flags
species1_noflagpts <- subset(species1allflags, distanceflag3sd!=1 & zooflag!=1 & landflag!=1 & climflag2!=1)

#Writes 2 files: one with all the points and flags and one without flagged points
write.csv(species1_noflagpts, file=paste(Genus, species,"_noflagpts",".csv", sep=""), row.names=FALSE)
write.csv(species1allflags, file=paste(Genus, species,"_allpts",".csv", sep=""), row.names=FALSE )


