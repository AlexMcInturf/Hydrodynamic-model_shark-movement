##### CURRENT AND SHARK SIMULATION for ANALYSIS #######
##### master script for running the code on all sharks ####

library(ncdf4)
library(maptools)
library(OpenStreetMap)
library(rgdal)

### here we define which shark we are working with
#sharkID <- 8446
#sharkID <- 8444
sharkID <- 8448

### define variables we want to work with


### now we can incorporate the shark into the code, first by defining the working directory
#setwd(paste0("~/Desktop/7 Gill Files/Shark",sharkID))
setwd(paste0("/Volumes/McInturf/7G_Final Analysis/Shark",sharkID))
getwd()


####################### PLOTTING CURRENT: ########################
# Note - do not need any Shark ID for this part of the code, because all of these files are present in the working directory 

#1. Need to incorporate netCDF files in correct format (files generated in Delft3d)
ncin <- nc_open("trim-ggate2.nc") #pulling in data from the netCDF file of the body of water we want to look at 
udat <- ncvar_get(ncin, varid="U1") #coordinates of current vector (east/west)
vdat <- ncvar_get(ncin, varid="V1") #coordinates of current vector velocities (north/south)
xcoor <- ncvar_get(ncin, varid="XCOR") #intersection of node in grid system
ycoor <- ncvar_get(ncin, varid="YCOR") #intersection of node in grid system
initialTimeVec <- ncvar_get(ncin, varid="time") #time within the simulation
alpha <- ncvar_get(ncin, varid="ALFAS") #corrects for the angle of every grid cell - alters the bias of the grid system
#east to west 
originTime <- as.POSIXct(strsplit(ncin[[14]][[28]]$dim[[4]]$units, split="seconds since ")[[1]][2], tz="GMT") #start time of simulation
attr(originTime, "tzone") <- "US/Pacific"
timevec <- originTime+ initialTimeVec #time vector = start time of simulation (IRL date/time)+ time within the simulation 

#some checks
timevec[1430:1441]
plot(xcoor[xcoor!=0],ycoor[xcoor!=0], cex=.2)
dim(udat)


#2. Need to convert current coordinates so can be plotted onto map

#create function to convert  between coordinate reference systems (crs)
conv.coord <- function(dd, crs1, crs2)
{
  temp <- SpatialPoints(dd, proj4string=CRS(crs1)) #create object containing information on CRS1
  data.frame(spTransform(temp, CRS(crs2))) #create dataframe from converting points from CRS1 to CRS2
}

#some checks
udat[1:10,1:10,1000]
udat[100:110,100:110,1000]
tmp <- udat[,,1000]
hist(tmp[tmp>0], nclass=100)
#3. Have readjusted time stamp so it reflects the correct time in the correct order for
#the simulation
#this is particularly important because we are sampling on every round number for the current data, 
#but our shark data takes place every second (does not line up without conversion)
# set up the timestamp to run a sample simulation (in order to view whether our simulation is working )

#timestamp8446 <-  "2008-09-09 13:15:00"
#timestamp8444 <- "2008-07-26 17:10:00"
timestamp8448 <- "2008-10-13 17:10:00"

#timestamp <- timestamp8446
#timestamp <- timestamp8444
timestamp <- timestamp8448

#4. Create a dataframe with all of the necessary coordinate info (both IRL and grid)
# NOTE: each time stamp will provide information on x,y,u,v, and alpha for each grid cell 
index <- 1000
cdatall <- data.frame(x=xcoor[xcoor!=0], y=ycoor[xcoor!=0], u=udat[,,index][xcoor!=0], v=vdat[,,index][xcoor!=0], alpha=alpha[xcoor!=0])
#in this case, we are randomly sampling data in order to run the simulation faster
subsampling <- 1
cdatallsub <- cdatall[sort(sample(1:nrow(cdatall), nrow(cdatall)/subsampling)),] 

#5. here we convert coordinates to correct time zone
#we will convert the coordinates of the first two columns of data (x,y) and then combine them with the rest of the data file (u,v, alpha) into one df
cdat <- data.frame(conv.coord(cdatallsub[,1:2], crs1="+proj=utm +zone=10 +north +datum=WGS84", 
                              crs2="+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"), 
                   cdatallsub[,3:5])

#in this case, we use alpha to convert the u and v velocity vector points to fit within the x,y grid system using alpha and trig
#we  define the function (v); v converts the data cdat 
arrowsize <- 500
arrowlength=0.02
arrowcol=2
arrowwidth=.5


u2 <- apply(cdat, 1, function(v) {v[3]*cos(-v[5]*pi/180)+v[4]*sin(-v[5]*pi/180)}) #v[3] and v[4] are u and v, while v[5] is the angle used in the trig
v2 <- apply(cdat, 1, function(v) {-v[3]*sin(-v[5]*pi/180)+v[4]*cos(-v[5]*pi/180)})
cdat <- cbind(cdat[,1:2], u2, v2, cdat[,3:5]) #then we combine u2 and v2 into the original dataframe with x,y,u,v, and alpha
head(cdat)




##6. here we create our map of the SF bay

#set the proportion of the map (note: these vary for each shark)

## shark 8446
#sizex <- 0.20
#sizey <- 0.1

#sharks 8444, 8448
sizex <- 0.1
sizey <- 0.07


#identify the center based on the mean of the columns of the x,y coordinates
center <- apply(cdatall[,1:2],2,mean)
#need to concatenate the objects in order to make them into a matrix with two rows (center with x and y), and then unlist so we can operate on them separately
#as two different objects
centerdeg <- unlist(c(conv.coord(matrix(center, ncol=2), crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=longlat +datum=WGS84"))) 
map <-  openmap(c(lat= centerdeg[2]+sizey,   lon= centerdeg[1]-sizex),c(lat= centerdeg[2]-sizey,   lon= centerdeg[1]+sizex), minNumTiles=9,type="bing")
plot(map, main="title")
arrows(cdat[,1], cdat[,2], cdat[,1]+ cdat[,3]*arrowsize, cdat[,2]+ cdat[,4]*arrowsize, angle=30, length=.02, col=2, lwd=.5)###pretty!
# argument 'angle' and 'length' describes the arrow head (not the angle or length of the arrow; these are done with the coordinates in argumen 3 and 4)  
#arrows(coordinates from which to draw (cdat 1 and 2), coordinates of points to which to draw (cdat 3 and 4), define the arrow features)

### now to map the grid system (for figure 2 in paper) #####
# first, need to get rid of alcatraz (grid points are actually overlaid on top)
grid <- data.frame(x=c(xcoor), y=c(ycoor)) #make a dataframe containing each vector of coordinates
head(grid)
gridMerc <- conv.coord(grid, crs1="+proj=utm +zone=10 +north +datum=WGS84", 
                              crs2="+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs") #convert all values to mercator
xMerc <- matrix(gridMerc[,1], nrow=dim(xcoor)[1], ncol=dim(xcoor)[2]) #now convert both x and y to a matrix
yMerc <- matrix(gridMerc[,2], nrow=dim(xcoor)[1], ncol=dim(xcoor)[2])
plot(xMerc, yMerc, cex=.1) # plot values

xMerc[xMerc<(-14000000)] <- NA #remove values where there is probably land (all listed as 0/NA in the vectors)
yMerc[yMerc<(1000)] <- NA

xMerc[sqrt(udat[,,1000]^2+vdat[,,1000]^2)<.001] <- NA
yMerc[sqrt(udat[,,1000]^2+vdat[,,1000]^2)<.001] <- NA

plot(xMerc, yMerc, cex=.1, xlim=c(-13650000, -13640000), ylim=c(4530000, 4540000)) #look at subset of map

#now to make the map (zoomed out)
center <- apply(grid[,1:2],2,mean)
sizex <- 0.22
sizey <- 0.22
map <-  openmap(c(lat= centerdeg[2]+sizey,   lon= centerdeg[1]-sizex),c(lat= centerdeg[2]-sizey,   lon= centerdeg[1]+sizex), minNumTiles=9,type="bing")
plot(map, main="title")
segments(xMerc[-nrow(xMerc),], yMerc[-nrow(xMerc),], xMerc[-1,], yMerc[-1,], col=2, lwd=.3)
segments(xMerc[,-ncol(xMerc)], yMerc[,-ncol(xMerc)], xMerc[,-1], yMerc[,-1], col=2, lwd=.3)

## map zoomed in
center <- apply(grid[,1:2],2,mean)
sizex <- 0.08
sizey <- 0.07
map <-  openmap(c(lat= centerdeg[2]+sizey,   lon= centerdeg[1]-sizex),c(lat= centerdeg[2]-sizey,   lon= centerdeg[1]+sizex), minNumTiles=9,type="bing")
plot(map, main="title")
segments(xMerc[-nrow(xMerc),], yMerc[-nrow(xMerc),], xMerc[-1,], yMerc[-1,], col=2, lwd=.7)
segments(xMerc[,-ncol(xMerc)], yMerc[,-ncol(xMerc)], xMerc[,-1], yMerc[,-1], col=2, lwd=.7)



#### now extract this image ###

pdf("Figure_X_map.pdf", width=8, height=16)# guidelines for j an ecol state the figures should be eps or pdf, 600dpi, 80mm wide or 1800px

map <-  openmap(c(lat= centerdeg[2]+sizey,   lon= centerdeg[1]-sizex),c(lat= centerdeg[2]-sizey,   lon= centerdeg[1]+sizex), minNumTiles=9,type="bing")
plot(map, main="title")
segments(xMerc[-nrow(xMerc),], yMerc[-nrow(xMerc),], xMerc[-1,], yMerc[-1,], col=2, lwd=.7)
segments(xMerc[,-ncol(xMerc)], yMerc[,-ncol(xMerc)], xMerc[,-1], yMerc[,-1], col=2, lwd=.7)

dev.off()


###7. Map the currents, by creating a function...again, have to convert time so that our times
#for plotting match up with the times for the current data
mapCurrent <- function(map, udat, vdat, alpha, xcoor, ycoor, timestamp, timevec, arrowsize, arrowcol, arrowlength, arrowwidth){
  timestamp2 <- as.POSIXct(timestamp, "%Y-%m-%d %H:%M:%S", tz="US/Pacific")
  minval <- floor(as.numeric(format(timestamp2+150, "%M"))/5)*5
  timestamp3 <- as.POSIXct(paste0(format(timestamp2+150, "%Y-%m-%d %H:"), minval, ":00"),tz="US/Pacific") #making sure that our values are all on the :00
  
  index <- which(format(timevec, "%Y-%m-%d %H:%M:%S")== format(timestamp3, "%Y-%m-%d %H:%M:%S"))
  cdatall <- data.frame(x=xcoor[xcoor!=0], y=ycoor[xcoor!=0], u=udat[,,index][xcoor!=0], v=vdat[,,index][xcoor!=0], alpha=alpha[xcoor!=0])
  cdatallsub <- cdatall[sort(sample(1:nrow(cdatall), nrow(cdatall))),]
  cdat <- data.frame(conv.coord(cdatallsub[,1:2], crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"), cdatallsub[,3:5])
  #hist(cdat$v[cdat$v>0], nclass=100)

  u2 <- apply(cdat, 1, function(v) {v[3]*cos(-v[5]*pi/180)+v[4]*sin(-v[5]*pi/180)})
  v2 <- apply(cdat, 1, function(v) {-v[3]*sin(-v[5]*pi/180)+v[4]*cos(-v[5]*pi/180)})
  cdat <- cbind(cdat[,1:2], u2, v2, cdat[,3:5])
  par(mar=c(0,0,0,0), oma=c(0,0,3,0))
  plot(map)
  mtext(timestamp, side=3, line=1, cex=1.5)
  arrows<- arrows(cdat[,1], cdat[,2], cdat[,1]+ cdat[,3]*arrowsize, cdat[,2]+ cdat[,4]*arrowsize, angle=30, length= arrowlength, col= arrowcol, lwd= arrowwidth)###pretty!
}

#double check to make sure that everything is working in the function (because everything is defined by timestamp (x,y,u,v,etc.) then we can leave most arguments as default except
#timestamp)

#8446
mapCurrent(map=map, udat=udat, vdat=vdat, alpha=alpha, xcoor=xcoor, ycoor=ycoor, timestamp="2008-09-11 05:20:23", timevec=timevec, arrowsize=500, arrowcol="red", arrowlength=0.02, arrowwidth=.5)
#8444
mapCurrent(map=map, udat=udat, vdat=vdat, alpha=alpha, xcoor=xcoor, ycoor=ycoor, timestamp="2008-07-30 12:02:35", timevec=timevec, arrowsize=500, arrowcol="red", arrowlength=0.02, arrowwidth=.5)
#8448
mapCurrent(map=map, udat=udat, vdat=vdat, alpha=alpha, xcoor=xcoor, ycoor=ycoor, timestamp="2008-10-13 17:10:00", timevec=timevec, arrowsize=500, arrowcol="red", arrowlength=0.02, arrowwidth=.5)


######## Interpolate Shark Data #############
#######if you already have data from interpolated points, skip this step and go further down to loading that data#############
########################################


#if you need to interpolate points:
dat <- read.csv(paste0("7G_", sharkID,".csv"))
head(dat, 100)


dat$timestamp <- as.POSIXlt(paste(dat$DATE, dat$TIME), "%m/%d/%y %H:%M:%S", tz="US/Pacific")

#conv.coord(dat[1,5:4], crs1="+proj=longlat +datum=WGS84" , crs2="+proj=utm +zone=10 +north +datum=WGS84") #convert the coordinates into UTM for easier measurements on map
dat[!is.na(dat[,4]),5:4] <- conv.coord(dat[!is.na(dat[,4]),5:4], crs1="+proj=longlat +datum=WGS84" , crs2="+proj=utm +zone=10 +north +datum=WGS84") #get rid of NAs

head(dat)
dat <- dat[order(dat$timestamp),]
plot(dat[,5:4], type="o", col=2, cex=.5, asp=1)

####sort data in chronological order

sum(is.na(dat$timestamp))##cool, no missing timestamp
dat1 <- dat[!is.na(dat$LAT),] #get rid of NAs in latitude
##calculate time intervals
head(dat1$timestamp)
interv <- dat1$timestamp[-1]-dat1$timestamp[-length(dat1$timestamp)] #calculate the difference between timestamps
interv2 <- as.numeric(interv) #make timestamp intervals numeric
dev.off() #do this to remove any limitations on plot size (allow us to see the full histogram)
hist(interv2[interv2<50], nclass=100)#
table(interv2)
#taking a look at the distribution of intervals in the data
dat1$interv <- c(NA, interv2) #shows the intervals between every timestamp
head(dat1[dat1$interv==0,], 50) #look at the data
dat1[41:99,]

#####remove timestamp replicates
dat2 <- dat1[!duplicated(dat1$timestamp),]
##now recalculate my time intervals
interv3 <- dat2$timestamp[-1]-dat2$timestamp[-length(dat2$timestamp)]
table(interv3)
plot(dat2[,5:4], type="o", col=2, cex=.5, asp=1)
dim(dat2)

##interpolate locations for every second
dat3 <- dat2[order(dat2$timestamp),] #order the data by timestamp
npts <- as.numeric(diff(range(dat3[,7])) , units="secs")+20000 #trying to determine the range of time in seconds for which we will interpolate data (plus a cushion)
res1 <- data.frame(id=character(npts), interptime=character(npts), #in this case, creating an empty dataframe in which we will fill interpolated points
                   interpx=numeric(npts), interpy=numeric(npts), depth=numeric(npts), stringsAsFactors=F)
iterator <- 1 ##how to keep track of time as the analysis is going 
for(j in 1:(nrow(dat3)-1)){ #going through all the rows in the dataframe, from 1 to n-1
  #for(j in 130:140){
  print(j/nrow(dat3)) #prints the percentage of the iteration that is finished
  temp <- dat3[j:(j+1),] #creates a temporary object that allows you to look at row x and the next row every time the function loops over a row
  diffTime <- trunc(as.numeric(temp[2,7]-temp[1,7], units="secs"))-1 #looks at the difference in time between those two rows (as above) - should be 0 if the values are 1 sec apart
  if(diffTime>0){ #if the difference in time between rows is greater than 0 (so interval greater than 1 sec), apply the following function to interpolate a point
    interpx <- seq(temp[1,5], temp[2,5], length.out= diffTime+2) #calculate a latitude in between that of the two rows (add 2 to diff time because including the true points on either end)
    interpy <- seq(temp[1,4], temp[2,4], length.out= diffTime+2) #calculate a longitude in between that of the two rows
    interptime <- seq(temp[1,7], temp[2,7], length.out= diffTime+2) #calculate a new time point htat is between that of the two rows 
    depth <- seq(temp[1,6], temp[2,6], length.out= diffTime+2) #calculate a new time point htat is between that of the two rows 
    #res will contain all of the rows: starting with the actual data from the iterator and adding in the interpolated points at the end into the data frame we have created
    res1[iterator:(iterator+length(interpx)-1),] <- data.frame(id=sharkID,interptime =format(interptime), interpx, interpy, depth, stringsAsFactors=F) #go through and fill in blanks with next row
    #put results of interpolation in res1 dataframe
    iterator <- iterator+length(interpx)
  }
}
res1 <- res1[res1$interpx!=0,] #res datframe will now only include points where interpolated points are not zero (don't need duplicates of existing points with correct 1 s interval)
dat33 <- dat3[,c(1,7,5,4,6)] #here we need to make sure that dat3 only contains the the columns that we want to match up to our res interpolated data 
names(res1) <- names(dat33) #the names from the original data are applied to the original+interpolated columns
lapply(res1, class)
res22 <- data.frame(res1,interp=TRUE) #create a new column, in which interpolated points are all labelled as true (however, these also contain the points on either side the interpolated ones)
dat4 <- data.frame(dat33, interp=FALSE, stringsAsFactors=F) #creates column in which interpolations are recorded as false (this is our "true" raw data)
lapply(dat4, class)
dat4$timestamp <- format(dat4$timestamp) 
res3 <- rbind(dat4, res22) #combine dat4 and res22, with dat4 on top so when we remove duplicates they will be removed from res (where they are currently recorded as false interpolations)
res3 <- res3[!duplicated(res3[,-6]),] #get rid of duplicates from res, minus the last column 
res3 <- res3[order(res3$timestamp),] #then re-order both interpolated and raw data based on timestamp, so that they are merged and correctly labeled as interpolated or not
head(res3)
res3[50:100,]

#for 8448, need to get rid of mistakes in data
res3 <- res3[res3$TRANSMITTER==8448,] #had a few 8446 and 8444
head(res3)
table(res3$TRANSMITTER)


write.csv(res3, paste0(sharkID,"_with_interpolated_points.csv"), row.names=F)

################# If already have interpolated points, then can move on to this step ###############
###########################################################################################################################
################# superimpose interpolated shark movement and current vectors ######################
###########################################################################################################################


shark <- read.csv(paste0(sharkID,"_with_interpolated_points.csv"))

head(shark)
shark$timestamp <- as.POSIXlt(shark$timestamp, tz="US/Pacific")
dim(shark)

###some checks
shark <- shark[order(shark$timestamp),]
interv <- shark $timestamp[-1]-shark $timestamp[-nrow(shark)]
table(interv)
which(interv==0)

table(shark$TRANSMITTER)


###########################
#Now let's plot everything
###########################
dim(shark)
xrange <- range(shark[,3])
yrange <- range(shark[,4])

#creating a function for color transparency on the tail
transp <- function(col, alphas=.5){
  temp <- rbind(col2rgb(col), alphas)
  res <- apply(temp,2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, c[4]))
  return(res)
}
####################
times <- shark$timestamp
shark2 <- shark

#here we need to convert coordinates from UTM to mercator, applied to lon and lat
shark2[,3:4] <- conv.coord(shark2[,3:4], crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs")


####################
####################
####################
####################
# for 8446
#sizex <- 0.10
#sizey <- 0.06
#identify the center based on the mean of the columns of the x,y coordinates
#center <- apply(shark[,3:4],2,mean)

#for other sharks:
sizex <- 0.13
sizey <- 0.08
#identify the center based on the mean of the columns of the x,y coordinates
center <- apply(cdatall[,1:2],2,mean)
#need to concatenate the objects in order to make them into a matrix with two rows (center with x and y), and then unlist so we can operate on them separately
#as two different objects

#need to concatenate the objects in order to make them into a matrix with two rows (center with x and y), and then unlist so we can operate on them separately
#as two different objects
centerdeg <- unlist(c(conv.coord(matrix(center, ncol=2), crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=longlat +datum=WGS84"))) 
#centerdeg <- centerdeg+c(-.01, -.008) #8446 
map <-  openmap(c(lat= centerdeg[2]+sizey,   lon= centerdeg[1]-sizex),c(lat= centerdeg[2]-sizey,   lon= centerdeg[1]+sizex), minNumTiles=9,type="bing")

# make a sample image 
jpeg(filename=paste0(sharkID,"sample.jpg"), width= 1600, height=1600* sizey/sizex+300)

par(mar=c(0,0,0,0), oma=c(0,0,3,0), fig=c(0,1, 0, 1))
plot(map, main="title")
arrows(cdat[,1], cdat[,2], cdat[,1]+ cdat[,3]*arrowsize, cdat[,2]+ cdat[,4]*arrowsize, angle=30, length=.02, col=2, lwd=.5)###pretty!
lines(shark2[,3:4], type='l', col="yellow")
par(fig = c(0.0,0.12, 0.75, 1), new = T ) 
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey", border="grey")
par(fig = c(0.0,0.12, 0.75, 1), new = T ) 
plot(1,1, type="n", xlab='', ylab='', xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1))
arrows(.5,.1,.5,.9, angle=90, length=.5, code=3)
text(.5,.1,pos=1, labels=paste0("-",max(shark2$DEPTH)," m"), cex=1.5, offset=.9)
text(.5,.9,pos=3, labels=paste0("0"," m"), cex=1.5, offset=.9)
dev.off()

####################
####################
####################
####################
####################
indices <- which(rep(c(rep(F, 49), T), length.out=length(times))) #rep() repeating a vector - in this case, only taking points every 50 
iter <- 1 #start at 1 in the dataframe and go row to row
for(i in indices){
#for(i in indices[1:5]){
  print(iter/length(indices)) #this gives information on where you are going 
    jpeg(filename=paste0(sharkID, '_movement1min_',formatC(iter,width = 6, format = "d", flag = "0"),'.png'), width= 1600, height=1600* sizey/sizex+300) #this provides information on where file will be saved
  #define the area of the plot htat will be filled by the current map
  plot(1,1, xlim=xrange, ylim=yrange, type="n", main=shark2$timestamp[i])
  # define the printed curent map
  mapCurrent(map=map, udat=udat, vdat=vdat, alpha=alpha, xcoor=xcoor, ycoor=ycoor, timestamp=shark2$timestamp[i], timevec=timevec, arrowsize=500, arrowcol="red", arrowlength=0.02, arrowwidth=.5)
  mtext(shark2$timestamp[i], side=3, line=1, cex=1.5)
  #going to have a tail for every 30 minutes, so here we define that
  minTime <- shark2$timestamp[i]-as.difftime("0:30:0")
  xy <- shark2[shark2$timestamp>minTime & shark2$timestamp<=shark2$timestamp[i],3:6] #this takes only the lat, long, interpolated data and depth
  if(nrow(xy)>1){
    coltail <- ifelse(xy[-nrow(xy),4],colors()[53] , "yellow")
    cols <- transp(coltail, alphas=seq(0.1,1, length.out=nrow(xy)-1))
    segments(xy[-nrow(xy),1], xy[-nrow(xy),2], xy[-1,1], xy[-1,2], col=cols, lwd=8)
  }
  if(nrow(xy)>0){
    points(xy[nrow(xy),1:2], col= "yellow", pch=16, cex=2) 
    par(fig = c(0.0,0.12, 0.75, 1), new = T ) 
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey", border="grey")
    par(fig = c(0.0,0.12, 0.75, 1), new = T ) 
    plot(1,1, type="n", xlab='', ylab='', xaxt="n", yaxt="n", xlim=c(0.2,1), ylim=c(0,1))
    arrows(.5,.1,.5,.9, angle=90, length=.1, code=3)
    text(.5,.1,pos=1, labels=paste0("-",max(shark2$DEPTH)," m"), cex=1.5, offset=1)
    text(.5,.9,pos=3, labels=paste0("0"," m"), cex=1.5, offset=1)
    depthy <- .9-(0.8*xy[nrow(xy),3]/max(shark2$DEPTH))
    points(0.5, depthy, col= colors()[132], pch=16, cex=3)
    text(.5,depthy,pos=4, labels=paste0("-",round(xy[nrow(xy),3],digits=1)," m"), cex=1.5, offset=1.5)
  }
  iter <- iter+1
  dev.off()
}

################################################
#################smoothing algo#################
################################################
head(shark)
plot(shark$timestamp, 1:length(shark$timestamp))##nice, looking for a point every second (no gaps in time)
library(devtools)
#install_github("edzer/trajectories")
library(trajectories)
library(spacetime)
shark1 <- shark[!is.na(shark$LAT),] #get rid of NAs in latitude
##take one point every 10 seconds
subsampling <- 10
timeseq <- seq(min(shark1[,2]), max(shark1[,2]), subsampling)
shark11 <- shark1[format(shark1[,2])%in%format(timeseq),] #taking out the points we want in the dataframe for ever 10 s

smth = function(x, y, xout,...) {
  predict(smooth.spline(as.numeric(x), y, ...), as.numeric(xout))
} #creates a function with our x,y values, then smooths them - makes them all numeric to do so
#now smooth actual data - using approxTrack to represent the track of trajectory data, with the smoothing function
#applied to dat2, defining size of the new track dataset (n) and the degree of smoothing
plot(shark[,3:4], asp=1, type="l")
stidf <- STIDF(SpatialPoints(shark11[,3:4],proj4string=CRS("+proj=utm +zone=10 +north +datum=WGS84")), time=shark11[,2], shark11[,c(1,5,6)])
shark2 <- Track(stidf)
shark.smooth <-  approxTrack(shark2, FUN = smth, n = length(shark2), spar=.5)##spar modifies the degree of smoothing
##plots
plot(shark2@sp@coords, type="l", asp=1) #tracks are s4 objects, which can be subset by @ instead of $. Here, we plot the coordintes of the track
lines(shark.smooth@sp@coords, col=2)
###
plot(shark2@sp[1500:1550]@coords, type="o", asp=1) # looking at smoothing at a smaller scale (hard to see at full-track scales)
lines(shark.smooth@sp[1500:1550]@coords, col=2)
length(shark2) #making sure the length of both the original and new dataframes are the same 
length(shark.smooth)
shark3 <- data.frame(shark11, shark.smooth@sp@coords)
names(shark3)[7:8] <- c("LON.smoothed", "LAT.smoothed")
head(shark3)

#install.packages("CircStats")
library(CircStats)

circAnglesCalc <- function(xy, interp){
	tmp <- cbind(xy[-nrow(xy),], xy[-1,]) #here, we are combining the rows we want for vectors
	                                      #we need to have one row of normal coordinates, and then shift the other row up to
	                                      # mark the first point of the vector and the point to which it is going 
	tmp <- tmp[!interp[-1] & !interp[-length(interp)],] #then we are getting rid of any interpolated data that may bias our results
	#cosAngle <- apply(tmp, 1, function(v) (v[1]*v[3]+v[2]*v[4])/(abs(v[3]-v[1])*(sqrt((v[3]-v[1])^2+(v[4]-v[2])^2))))
	angle <- apply(tmp, 1, function(v) { #we are going to apply function v to calculate the arc tan angle 
		v0 <- c(1,0) #v0 are the coordinates of the starting point of the vector
		v1 <- c(v[3]-v[1],v[4]-v[2]) #v1 is the difference between the first point of the vector and the second
		x <- atan2(v0[1]*v1[2]-v0[2]*v1[1], v0[1]*v1[1]+v0[2]*v1[2]) #here we calculate arc tan between the two points
		if(x <0) x <- x +2*pi #now we need to normalize, make sure that we get all of the negative numbers
		return(x)
		}
		)
	return(angle)
}

angles <- circAnglesCalc(xy=shark3[,7:8], interp=shark3$interp) #here is our function to calculate angles for the rose diagram
rose.diag(angles, 100, prop=5) # here's our rose diagram; angles, number of bins, prop= radius
hist(angles, nclass=100) #histogram to look at distribution of angles, to make sure nothing is off 
table(angles)

res <- hist(angles, nclass=100) #need to save this as a histogram because it has the frequency of angles (essentially want to wrap the histogram around in a circle)

# library(ggplot2)
# 
# ggplot(data.frame(x=1:length(res$counts), y=res$counts), aes(x,y)) +
#   geom_bar(stat="identity", col=2, fill=4) + #change color by altering the color of lines and fill of bins in plot
#   coord_polar(start = -pi/2, direction=-1) + #need to adjust the axis to reflect what we are looking for (going counter-clockwise, and starting at 0 deg)
#   scale_x_continuous("",breaks = seq(1, length(res$counts)))+
# theme(axis.text = element_blank(), # gets rid of numbers of axes
#       axis.ticks = element_blank())
      #panel.grid  = element_blank())

#########################
#### now to determine step length and break up into longer lengths and shorter lengths #####
#### Barnett et al., 2010: average of actively tracked cow sharks = .48 m/s #####
#########################
head(shark3) #using data from smoothing algorithm

# looking to set up a data frame with the initial point of the vector and the next point to determine the distance between them
# ANNA CODE
#sharkvec1 <- data.frame(timestamp=shark[-nrow(shark),2], x=shark[-1,3]-shark[-nrow(shark),3], y=shark[-1,4]-shark[-nrow(shark),4], interv=as.numeric(shark[-1,2]-shark[-nrow(shark),2]))
#head(sharkvec1)
# this shows the difference between each x and y value (or lat and long values) for each timestamp, with an interval of 10 seconds because the data is interpolated
#and smoothed

# Alex Breakdown
sharkvec2 <- data.frame(timestamp=shark3[-nrow(shark3),2], # timestamp is the timestamp of the shark
                        x2 = shark3[-1,7], #shark minus the first row, all of lon column (moving them up a row)
                        x1 = shark3[-nrow(shark3),7], #shark minus the last row, all of lon column
                        y2 = shark3[-1,8], # minus the first row, all of lat column
                        y1 = shark3[-nrow(shark3),8], 
                       # depth = shark3[-nrow(shark3),5], #depth at the beginning of the vector (x1 and y1 location)
                        interp1=shark3[-1,6], #need to also include interp data, so these correspond to the x2 and y2 data
                        interp2=shark3[-nrow(shark3),6],
                        angles = circAnglesCalc(xy=shark3[,7:8], interp=rep(FALSE,nrow(shark3))) #need to list these as false to "trick" R into thinking there isn't interpolated data
                                                                                                # circAngles takes out interpolated data, but we want to wait until after we calculate distance
                        ) #these correspond to the x1 and y1 data
# this provides the x and y coordinates of each arrow, so then we will need to find the distance between the two using sqrt((x2-x1)^2 + (y2-y1)^2)
# in this case, because data is interpolated, each second is the interval over which this distance occurs
 #here is our function to calculate angles for the rose diagram
head(sharkvec2)

dist <- apply(sharkvec2[,2:5], 1, function(v) {sqrt((v[1]-v[2])^2+(v[3]-v[4])^2)}) #we were getting a "non-numeric argument to binary operator" argument
#so, in the above dataframe, needed to specify the columns we were operating on, which were all numeric, so they didn't convert to character due to presence of timestamp
#need to make each component of the data frame numeric in order to conduct mathematical operations
shark4 <- cbind(sharkvec2, dist)
head(shark4)
sharkvec3 <- shark4[!shark4$interp1==TRUE & !shark4$interp2==TRUE,] #then we are getting rid of any interpolated data that may bias our results
histDist <- hist(sharkvec3$dist, nclass=100, xlim =c(0,50))
rose.diag(sharkvec3$angles, 100, prop=5) # here's our rose diagram; angles, number of bins, prop= radius

#### trying to separate out angles to plot based on the distance of the step length (shorter localized movements, longer distance movemements)
par(mfrow=c(1,3))
rose.diag(sharkvec3$angles[sharkvec3$dist<quantile(sharkvec3$dist, .15)], 100, prop=5) #plotting rose diagrams from 0 - 15th quantile
rose.diag(sharkvec3$angles[sharkvec3$dist>=quantile(sharkvec3$dist, .15) & #plotting rose diagrams from 15-median
                             sharkvec3$dist<quantile(sharkvec3$dist, .50)], 100, prop=5) # here's our rose diagram; angles, number of bins, prop= radius
rose.diag(sharkvec3$angles[sharkvec3$dist>=median(sharkvec3$dist)], 100, prop=5) # here's our rose diagram; angles, number of bins, prop= radius

########################################################################################################################
######### now need to set up current code in order to get direction of current for each timestamp and location #####
########################################################################################################################

getCurrent <- function(udat, vdat, alpha, xcoor, ycoor, timestamp, timevec, shark.x, shark.y){ 
  timestamp2 <- as.POSIXct(timestamp, "%Y-%m-%d %H:%M:%S", tz="US/Pacific")
  minval <- floor(as.numeric(format(timestamp2+150, "%M"))/5)*5
  timestamp3 <- as.POSIXct(paste0(format(timestamp2+150, "%Y-%m-%d %H:"), minval, ":00"),tz="US/Pacific") #making sure that our values are all on the :00, reformatting as expected
  
  index <- which(format(timevec, "%Y-%m-%d %H:%M:%S")== format(timestamp3, "%Y-%m-%d %H:%M:%S"))
  cdatall <- data.frame(x=xcoor[xcoor!=0], y=ycoor[xcoor!=0], u=udat[,,index][xcoor!=0], v=vdat[,,index][xcoor!=0], alpha=alpha[xcoor!=0]) #again, reformatting as above
  
  u2 <- apply(cdatall, 1, function(v) {v[3]*cos(-v[5]*pi/180)+v[4]*sin(-v[5]*pi/180)}) #calculating the angles of the current
  v2 <- apply(cdatall, 1, function(v) {-v[3]*sin(-v[5]*pi/180)+v[4]*cos(-v[5]*pi/180)})
  cdat <- cbind(cdatall[,1:2], u2, v2, cdatall[,3:5]) #creating a new df with x,y, u2, v2, and vector information
  
  distVec <- apply(cdat[,1:2], 1, function(v) sqrt((v[1]-shark.x)^2+(v[2]-shark.y)^2)) #now need to calculate the distance between shark vector and current vector  
  return(cdat[order(distVec)[1],3:4]) #order by least distance between shark vector and grid point so as to find where to compare vectors within the grid
}

# a test
getCurrent(udat=udat, vdat=vdat, alpha=alpha, xcoor=xcoor, ycoor=ycoor, timestamp="2008-10-14 17:10:00", timevec=timevec, shark.x=543828.3, shark.y=4162973)
#this seems to work

sharkvec4 <- cbind(sharkvec3, data.frame(current.x=NA, current.y=NA)) #now need to add  current vectors to sharkvec3 data frame
for(i in 1:nrow(sharkvec4)){ #fill in data frame with current vector points that correspond to shark location
#for(i in 1:100){
  print(i/nrow(sharkvec4))
  sharkvec4[i,10:11] <- getCurrent(udat=udat, vdat=vdat, alpha=alpha, xcoor=xcoor, ycoor=ycoor, timestamp=sharkvec4[i,1], timevec=timevec, shark.x=sharkvec4[i,3], shark.y=sharkvec4[i,5])
}
  

##### Now trying to find out if sharks are ACTIVELY moving with or against the current #############
# need to create dataframe now with new columns that sutract current vector from shark vector
sharkvec5 <- cbind(sharkvec4, shark.x.minus.curr=NA, shark.y.minus.curr=NA)
#create shark vectors alone (no current) - note: this an important step in getting a scalar product
# note: current.x is the entire vector, not a point? 
sharkvec5$shark.x.minus.curr <- (sharkvec5$x2- sharkvec5$x1)/10-sharkvec5$current.x
sharkvec5$shark.y.minus.curr <- (sharkvec5$y2- sharkvec5$y1)/10-sharkvec5$current.y
head(sharkvec5)
plot((sharkvec5$x2-sharkvec5$x1)/10, sharkvec5$current.x)
hist(sharkvec5$shark.x.minus.curr) #note: there are lots of 0, so likely that the shark went with the current in many places
hist(sharkvec5$shark.y.minus.curr)


# we will then calculate the angles of the shark WITHOUT the current present and plot these on a rose diagram
angle.minus.current <- apply(sharkvec5[,12:13], 1, function(v) { #we are going to apply function v to calculate the arc tan angle 
  v0 <- c(1,0) #v0 are the coordinates of the starting point of the vector, where x is the length of the vector (1) and y is 0
  v1 <- v #v1 is the difference between the first point of the vector and the second
  x <- atan2(v0[1]*v1[2]-v0[2]*v1[1], v0[1]*v1[1]+v0[2]*v1[2]) #here we calculate arc tan between the two points
  if(x <0) x <- x +2*pi #now we need to normalize, make sure that we get all of the negative numbers
  return(x)
}
)

########################################

head(angle.minus.current)
hist(angle.minus.current)
rose.diag(angle.minus.current, 100, prop=3) 
rose.diag(angles, 100, prop=3) 
plot(angle.minus.current, angles)

# Now plotting just the current vectors, no sharks
current.angle <- apply(sharkvec5[,10:11], 1, function(v) { #we are going to apply function v to calculate the arc tan angle 
  v0 <- c(1,0) #v0 are the coordinates of the starting point of the vector, where x is the length of the vector (1) and y is 0
  v1 <- v #v1 is the difference between the first point of the vector and the second
  x <- atan2(v0[1]*v1[2]-v0[2]*v1[1], v0[1]*v1[1]+v0[2]*v1[2]) #here we calculate arc tan between the two points
  if(x <0) x <- x +2*pi #now we need to normalize, make sure that we get all of the negative numbers
  return(x)
}
)

head(current.angle)
summary(current.angle)
rose.diag(current.angle, 100, prop=2)
sharkvec6 <- cbind(sharkvec5, current.angle) #add current angle to the sharkvec dataframe
head(sharkvec6)

######## Now let's do a scatterplot of shark angles with the current and the current itself ##########
plot(sharkvec6$angles,sharkvec6$current.angle, pch=16) 
sharkvec6$current.speed <- sqrt((sharkvec6$current.x)^2 + (sharkvec6$current.y)^2) #now calculate the distance of each current vector and add as "current speed" into the dataframe (meters/second)
head(sharkvec6)

# plotting according to current speed, so plotting current angle and shark angle when current speed is greater than the median speed or less than the median speed
plot(current.angle~angles,data=sharkvec6[sharkvec6$current.speed>median(sharkvec6$current.speed),]) 
plot(current.angle~angles,data=sharkvec6[sharkvec6$current.speed<median(sharkvec6$current.speed),])

# now plotting the cosine, which provides the correlation between angles
# if the angles are entirely correlated, there should be no difference between them and cosine will thus be + 1. If they are in opposite directions, cosine = -1 
plot(cos(current.angle)~cos(angles),data=sharkvec6[sharkvec6$current.speed>median(sharkvec6$current.speed),]) 
#plot(sin(current.angle)~sin(angles),data=sharkvec6[sharkvec6$current.speed>median(sharkvec6$current.speed),]) 

plot(cos(current.angle)~cos(angles),data=sharkvec6[sharkvec6$current.speed<median(sharkvec6$current.speed),]) 
#plot(sin(current.angle)~sin(angles),data=sharkvec6[sharkvec6$current.speed<median(sharkvec6$current.speed),]) 

#scalar product of current and shark vectors
# going to use scalar product to multiple each vector by the cosine of the angle between them
# again, this will tell us how much the vector angles agree 
#(if cosine is 1, and both vectors are going quickly in the same direction, would expect to see highest numbers)
#(if cosine is -1 and both vectors are going quickly in the opposite direction, would expect to see highly negative numbers)
sharkvec6$shark.vec.x <- sharkvec6$x2-sharkvec6$x1
sharkvec6$shark.vec.y <- sharkvec6$y2-sharkvec6$y1
sharkvec6$dot.product <- sharkvec6$shark.vec.x*sharkvec6$current.x + sharkvec6$shark.vec.y*sharkvec6$current.y

#create a histogram of the dot.product with a line at zero to look at the distribution of angles around 0 (hope to see a lot of angles that are 0)
hist(sharkvec6$dot.product, nclass=100)
segments(0,0,0,1700, lwd=2, col=2, lty=2)


#now we decide to look at the cosine once again, to find a correlation between angles
# to do this, we subtract the current angle from the shark angle and find the cosine of the difference
hist(cos(sharkvec6$angles-sharkvec6$current.angle), nclass=100)
## ideally, almost all angles have cosine of 1, positively correlated

#now comparing the cosines of the angles for current speeds that are higher vs. lower current speeds
hist(cos((sharkvec6$angles-sharkvec6$current.angle)[sharkvec6$current.speed>median(sharkvec6$current.speed)]), nclass=100)
hist(cos((sharkvec6$angles-sharkvec6$current.angle)[sharkvec6$current.speed<median(sharkvec6$current.speed)]), nclass=100)
#using a rose diagram to show that the difference b/w the angles is largely 0 
rose.diag((sharkvec6$angles-sharkvec6$current.angle)[sharkvec6$current.speed>median(sharkvec6$current.speed)], 100, prop=4)
rose.diag((sharkvec6$angles-sharkvec6$current.angle)[sharkvec6$current.speed<median(sharkvec6$current.speed)], 100, prop=4)


# finding median of all sharks
#dat1 <- read.csv("currentspeed_8448.csv")
#dat2 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8444/currentspeed_8444.csv")
#dat3 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8446/currentspeed_8446.csv")
#head(dat1)
#head(dat2)
#head(dat3)

#removing first column
#dat1 <- dat1[,-1]
#dat2 <- dat2[,-1]
#dat3 <- dat3[,-1]

#full.dat <- c(dat1,dat2,dat3)
#head(full.dat)
#tail(full.dat)
#tail(dat3)
#head(dat1)

#double check
#length(dat1)
#length(dat2)
#length(dat3)
#sum(length(dat1), length(dat2), length(dat3))
#length(full.dat)

#summary(full.dat)

median.curr <- 0.452903 #calculated from the distribution of all current speeds from all sharks


################### FINAL ROSE PLOTS AND DIAGRAMS! ################
###################################################################
#write.csv(sharkvec6$current.speed, "currentspeed_8446.csv")
#write.csv(sharkvec6$current.speed, "currentspeed_8444.csv")
#write.csv(sharkvec6$current.speed, "currentspeed_8448.csv")



####################################
## first look at faster current speed
rose.diag((angle.minus.current-sharkvec6$current.angle)[sharkvec6$current.speed>median.curr], 100, prop=4) # lots of -180, likely by chance...fighting the current? trying to forage in face of current, or selecting areas in which the current is not as strong as the average so it seems like the shark is going against the current 
rose.diag((angle.minus.current)[sharkvec6$current.speed>median.curr], 100, prop=4)
rose.diag((sharkvec6$current.angle)[sharkvec6$current.speed>median.curr], 100, prop=2)
rose.diag((angles)[sharkvec6$current.speed>median.curr], 100, prop=4)
plot((angle.minus.current)[sharkvec6$current.speed>median.curr], (sharkvec6$current.angle)[sharkvec6$current.speed>median.curr],pch=16)
# no real pattern, looks like shark may be drifting


## now look at lower speeds (below median shark speed)
rose.diag((angle.minus.current-sharkvec6$current.angle)[sharkvec6$current.speed<median.curr], 100, prop=4) # for smaller current speeds, looks random
rose.diag((angle.minus.current)[sharkvec6$current.speed<median.curr], 100, prop=4) # shark does not swim randomly at lower speeds
rose.diag((sharkvec6$current.angle)[sharkvec6$current.speed<median.curr], 100, prop=3)
rose.diag((angles)[sharkvec6$current.speed<median.curr], 100, prop=4)
 plot((angle.minus.current)[sharkvec6$current.speed<median.curr], (sharkvec6$current.angle)[sharkvec6$current.speed<median.curr],pch=16)

 

### plotting cosine versus sine for the angle minus the current compared to the current angle
par(mfrow=c(2,1))
# comparing sin/cosin of slow current plots
plot(cos(angle.minus.current)[sharkvec6$current.speed<median(sharkvec6$current.speed)], (cos(sharkvec6$current.angle)[sharkvec6$current.speed<median(sharkvec6$current.speed)]),pch=16) #slow current
plot(sin(angle.minus.current)[sharkvec6$current.speed<median(sharkvec6$current.speed)], (sin(sharkvec6$current.angle)[sharkvec6$current.speed<median(sharkvec6$current.speed)]),pch=16)

# of fast current plots
plot(cos(angle.minus.current)[sharkvec6$current.speed>median(sharkvec6$current.speed)], (cos(sharkvec6$current.angle)[sharkvec6$current.speed>median(sharkvec6$current.speed)]),pch=16) #fast current
plot(sin(angle.minus.current)[sharkvec6$current.speed>median(sharkvec6$current.speed)], (sin(sharkvec6$current.angle)[sharkvec6$current.speed>median(sharkvec6$current.speed)]),pch=16)

# of slow current rose diagrams
rose.diag(cos(angle.minus.current)[sharkvec6$current.speed<median(sharkvec6$current.speed)], 100,prop=2) #slow current
rose.diag(sin(angle.minus.current)[sharkvec6$current.speed<median(sharkvec6$current.speed)], 100,prop=2)

# of fast current rose diagrams
rose.diag(cos(angle.minus.current)[sharkvec6$current.speed>median(sharkvec6$current.speed)], 100, prop=2) #fast current
rose.diag(sin(angle.minus.current)[sharkvec6$current.speed>median(sharkvec6$current.speed)], 100, prop=2)

##### Final movement numbers #####
## first, to determine mean/range of shark speed
# need to find distance/10 seconds for meters/sec
head(sharkvec6)
sharkvec6$sharkspeed1 <- sharkvec6$dist/10 #not corrected for current 
summary(sharkvec6$sharkspeed1)
summary(sharkvec6$sharkspeed1[sharkvec6$current.speed>median.curr])
summary(sharkvec6$sharkspeed1[sharkvec6$current.speed<median.curr])
sharkvec6$sharkspeed2 <- (sharkvec6$dist-sharkvec6$current.speed)/10 #corrected for current
summary(sharkvec6$sharkspeed2)
summary(sharkvec6$sharkspeed2[sharkvec6$current.speed>median.curr])
summary(sharkvec6$sharkspeed2[sharkvec6$current.speed<median.curr])

## now try to do angles (of shark)


summary(sharkvec6$sharkspeed)
summary(sharkvec6$sharkspeed[sharkvec6$current.speed>median(sharkvec6$current.speed)]) #look at metrics when currents are fast
summary(sharkvec6$sharkspeed[sharkvec6$current.speed<median(sharkvec6$current.speed)]) # look at metrics when currents are slow 


############### SAVE ALL DATA AS .csv ####################
write.csv(sharkvec6, paste0("/Volumes/McInturf/7G_Final Analysis/Shark", sharkID, "_fulldata.csv"))

############## FINAL FIGURES #############
# first, a 3x2 panel of slow vs. fast
# a = shark bearing (uncorrected) - current
# b = shark bearing (corrected) - current
# c = current vectors

install.packages("circular")
library(circular)

## Fig 2

### alex = fix median and repeat for other sharks
tiff("Figure_2_rosediags_8446.tiff", width=8, height=16, units="cm", res=600, compression="lzw")# guidelines for j an ecol state the figures should be eps or pdf, 600dpi, 80mm wide or 1800px
pdf("Figure_2_rosediags_8446.pdf", width=8, height=16, units="cm")# guidelines for j an ecol state the figures should be eps or pdf, 600dpi, 80mm wide or 1800px

par(mfcol=c(4,2), mar=c(1,1,1,1))
# faster than median (shark w/ current, current alone, shark swimming movement in current, difference b/w shark and current)
sharktrue <- (sharkvec6$angles)[sharkvec6$current.speed>median.curr]
str(sharktrue)
sharktrue = circular(sharktrue, units="rad")
rose.diag(sharktrue, bins=45, shrink=1, prop=1.8, axes=F)
  
currentonly <- (sharkvec6$current.angle)[sharkvec6$current.speed>median.curr]
str(currentonly)             
currentonly = circular(currentonly, units="rad")
rose.diag(currentonly, bins=45, shrink=1, prop=1.8, axes=F)  
  
sharkmbear1 <- (sharkvec6$angles-sharkvec6$current.angle)[sharkvec6$current.speed>median.curr]
str(sharkmbear1)             
sharkmbear1 = circular(sharkmbear1, units="rad")
rose.diag(sharkmbear1, bins=45, shrink=1, prop=1.8, axes=F)  
  
sharkmbear2 <- deg((angle.minus.current-sharkvec6$current.angle)[sharkvec6$current.speed>median.curr])
str(sharkmbear2)             
sharkmbear2 = circular(sharkmbear2, units="deg")
rose.diag(sharkmbear2, bins=45, shrink=1, prop=1.8, tcl.text=0.2)  
  
  # slower than median 
sharktrue <- (sharkvec6$angles)[sharkvec6$current.speed<median.curr]
str(sharktrue)
sharktrue = circular(sharktrue, units="rad")
rose.diag(sharktrue, bins=45, shrink=1, prop=1.8, axes=F)
  
currentonly <- (sharkvec6$current.angle)[sharkvec6$current.speed<median.curr]
str(currentonly)             
currentonly = circular(currentonly, units="rad")
rose.diag(currentonly, bins=45, shrink=1, prop=1.8, axes=F)  
  
sharkmbear1 <- (sharkvec6$angles-sharkvec6$current.angle)[sharkvec6$current.speed<median.curr]
str(sharkmbear1)             
sharkmbear1 = circular(sharkmbear1, units="rad")
rose.diag(sharkmbear1, bins=45, shrink=1, prop=1.8, axes=F)  
  
sharkmbear2 <- deg((angle.minus.current-sharkvec6$current.angle)[sharkvec6$current.speed<median.curr])
str(sharkmbear2)             
sharkmbear2 = circular(sharkmbear2, units="deg")
rose.diag(sharkmbear2, bins=45, shrink=1, prop=1.8, tcl.text=0.2)  
  
dev.off()

install.packages("ggplot2")
library(ggplot2)

### Fig 3 ###
# to calculate amount of time interpolated vs. non-interpolated
table(shark$interp)
total.seconds <- length(shark$interp)
total.seconds
interp.seconds <- sum(shark$interp[shark$interp=="TRUE"])
real.seconds <- total.seconds - interp.seconds
# convert to minutes
percent.interp <- interp.seconds/total.seconds
percent.interp
percent.real <- real.seconds/total.seconds
percent.real

# make for each shark, rbind the three trashdf into a single df; if you want to change the way you refernto sharks, put that in here instead of the current shark ID column
pdf("Figure_3_interppts_8444.pdf", width=8, height=16)# guidelines for j an ecol state

#shark8448pts = data.frame(sharkid = rep(3,nrow(sharkvec6)), truetimestamp = sharkvec6$timestamp, expttimestamp.10min = as.numeric(difftime(sharkvec6$timestamp, sharkvec6$timestamp[1]))/(10*60))
#write.csv(shark8448pts, "shark8448pts.csv")

#shark8444pts = data.frame(sharkid = rep(1,nrow(sharkvec6)), truetimestamp = sharkvec6$timestamp, expttimestamp.10min = as.numeric(difftime(sharkvec6$timestamp, sharkvec6$timestamp[1]))/(10*60))
#write.csv(shark8444pts, "shark8444pts.csv")

#shark8446pts = data.frame(sharkid = rep(2,nrow(sharkvec6)), truetimestamp = sharkvec6$timestamp, expttimestamp.10min = as.numeric(difftime(sharkvec6$timestamp, sharkvec6$timestamp[1]))/(10*60))
#write.csv(shark8446pts, "shark8446pts.csv")

pts8444 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8444/shark8444pts.csv")
pts8446 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8446/shark8446pts.csv")
pts8448 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8448/shark8448pts.csv")

pts8444 <- pts8444[,-1]
pts8446 <- pts8446[,-1]
pts8448 <- pts8448[,-1]

pts.full <- rbind(pts8444,pts8446, pts8448)
head(pts.full)
tail(pts.full)

ggplot(data=pts.full, aes(x=factor(sharkid), y=expttimestamp.10min)) + geom_point(size=.2) + coord_flip() + ylab("Duration of Tracking (10 min)") + xlab("Shark Identification") + theme_bw() +theme(panel.grid.minor = element_blank()) + theme(text = element_text(size=12)) 
ggsave("Figure_3_TrackConsistency.pdf", scale=1.5, width=8, height=5)

#### Fig 4 ####

# of fast current plots
#fast.true = ggplot(data=sharkvec6, aes(x=angle.minus.current, y=current.angle, color=current.speed)) + geom_point(pch=16) + scale_color_viridis() + ylim(0,2*pi) + theme_bw()#fast current
install.packages("viridis")
library(viridis)

cosdata.shark = cos(angle.minus.current)[sharkvec6$current.speed>median.curr]
cosdata.cur = cos(sharkvec6$current.angle)[sharkvec6$current.speed>median.curr]
cosdata = data.frame(cosdata.shark, cosdata.cur, Speed = sharkvec6$current.speed[sharkvec6$current.speed>median.curr])
fast.cos = ggplot(data=cosdata, aes(x=cosdata.shark, y=cosdata.cur, color=Speed)) + geom_point(pch=16) + scale_color_gradient() + ylim(-1,1) + theme_bw() + xlab("Cos of Shark Movement") + ylab("Current Movement (Fast)") + labs(colour = "") #fast current

sindata.shark = sin(angle.minus.current)[sharkvec6$current.speed>median.curr]
sindata.cur = sin(sharkvec6$current.angle)[sharkvec6$current.speed>median.curr]
sindata = data.frame(sindata.shark, sindata.cur, Speed = sharkvec6$current.speed[sharkvec6$current.speed>median.curr])
fast.sin = ggplot(data=sindata, aes(x=sindata.shark, y=sindata.cur, color=Speed)) + geom_point(pch=16) + scale_color_gradient() + ylim(-1,1) + theme_bw() + xlab("Sin of Shark Movement") + ylab("") + labs(colour = "")  # to make grey scale, check scale_color_continous() or scale_color_bnrewer()


cosdatas.shark = cos(angle.minus.current)[sharkvec6$current.speed<median.curr]
cosdatas.cur = cos(sharkvec6$current.angle)[sharkvec6$current.speed<median.curr]
cosdatas = data.frame(cosdatas.shark, cosdatas.cur, Speed = sharkvec6$current.speed[sharkvec6$current.speed<median.curr])
slow.cos = ggplot(data=cosdatas, aes(x=cosdatas.shark, y=cosdatas.cur, color=Speed)) + geom_point(pch=16) + scale_color_gradient() + ylim(-1,1) + theme_bw() + xlab("Cos of Shark Movement") + ylab("Current Movement (Slow)") + labs(colour = "")

sindatas.shark = sin(angle.minus.current)[sharkvec6$current.speed<median.curr]
sindatas.cur = sin(sharkvec6$current.angle)[sharkvec6$current.speed<median.curr]
sindatas = data.frame(sindatas.shark, sindatas.cur, Speed = sharkvec6$current.speed[sharkvec6$current.speed<median.curr])
slow.sin = ggplot(data=sindatas, aes(x=sindatas.shark, y=sindatas.cur, color=Speed)) + geom_point(pch=16) + scale_color_gradient() + ylim(-1,1) + theme_bw() + xlab("Sin of Shark Movement") + ylab("") + labs(colour = "")

# sindata.shark = sin(angle.minus.current)
# sindata.cur = sin(sharkvec6$current.angle)
# sindata = data.frame(sindata.shark, sindata.cur, curspd = sharkvec6$current.speed)
# all.sin = ggplot(data=sindata, aes(x=sindata.shark, y=sindata.cur, color=curspd)) + geom_point(pch=16) + scale_color_viridis() + ylim(-1,1) + theme_bw()#fast current

#install.packages("cowplot")
library(cowplot)  

tiff("Figure_4_cossin_8446.tiff", width=16, height=16, units="cm", res=600, compression="lzw")# guidelines for j an ecol state
plot_grid(fast.cos, fast.sin, slow.cos, slow.sin)
dev.off()

############### Fig. 1 Panels ###################
head(sharkvec6)
## note: these do NOT have interpolated points
track8444 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8444_fulldata.csv")
track8446 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8446_fulldata.csv")
track8448 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8448_fulldata.csv")


# need to remove interpolated data
# coordinates, sharkID, current speed
sharkvec7_8444 <- track8444[,c("timestamp","x1","y1", "current.speed", "sharkspeed1")]
head(sharkvec7_8444)
sharkvec7_8446 <- track8446[,c("timestamp","x1","y1", "current.speed", "sharkspeed1")]
head(sharkvec7_8448)
sharkvec7_8448 <- track8448[,c("timestamp","x1","y1", "current.speed", "sharkspeed1")]
head(sharkvec7_8448)

# create a new column that will break the track into groups between missing segments
sharkvec7_8444$timedif <- as.numeric(difftime(sharkvec7_8444$timestamp,lag(sharkvec7_8444$timestamp)))
sharkvec7_8444$trackseg <- c(NA,(cumsum(sharkvec7_8444$timedif[2:nrow(sharkvec7_8444)]>600)+1) )

sharkvec7_8446$timedif <- as.numeric(difftime(sharkvec7_8446$timestamp,lag(sharkvec7_8446$timestamp)))
sharkvec7_8446$trackseg <- c(NA,(cumsum(sharkvec7_8446$timedif[2:nrow(sharkvec7_8446)]>600)+1) )

sharkvec7_8448$timedif <- as.numeric(difftime(sharkvec7_8448$timestamp,lag(sharkvec7_8448$timestamp)))
sharkvec7_8448$trackseg <- c(NA,(cumsum(sharkvec7_8448$timedif[2:nrow(sharkvec7_8448)]>600)+1) )


## map 
library(ggmap)
library(viridis)
library(cowplot)

#sizex <- 0.0002
#sizey <- 0.0001
center <- apply(sharkvec7_8446[,2:3],2,mean)
centerdeg <- unlist(c(conv.coord(matrix(center, ncol=2), crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=longlat +datum=WGS84"))) 
centerdeg


center = c(-122.500, 37.820)

sharkvec7_8444[,c("lon","lat")] = conv.coord(sharkvec7_8444[,2:3], crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=longlat +datum=WGS84")
sharkvec7_8446[,c("lon","lat")] = conv.coord(sharkvec7_8446[,2:3], crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=longlat +datum=WGS84")
sharkvec7_8448[,c("lon","lat")] = conv.coord(sharkvec7_8448[,2:3], crs1="+proj=utm +zone=10 +north +datum=WGS84", crs2="+proj=longlat +datum=WGS84")



map = get_map(location=c(lon=center[1], lat=center[2]), maptype=c("satellite"), zoom=12)
maplayer = ggmap(map)

shark8448map <- maplayer + geom_point(data=sharkvec7_8448, aes(x=lon, y=lat, group=trackseg, color=trackseg), size=.5) + 
  #scale_color_viridis(direction=-1, option="A", begin=0.1, end=0.95) + 
  xlab("Longitude") +  scale_color_gradient(high="yellow", low="yellow2") + 
  ylab("") + #ylab("Latitude") +
  theme(legend.position="none")
shark8446map <- maplayer + geom_point(data=sharkvec7_8446, aes(x=lon, y=lat, group=trackseg, color=trackseg), size=.5) + 
  scale_color_gradient(high="cyan", low="cyan4") + 
  xlab("Longitude") + 
  ylab("") + #ylab("Latitude") + 
  theme(legend.position="none")
shark8444map <- maplayer + geom_point(data=sharkvec7_8444, aes(x=lon, y=lat, group=trackseg, color=trackseg), size=.5) + 
  scale_color_gradient(high="palegreen", low="palegreen2") +  
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(legend.position="none") 

pdf("Figure_1_panel.pdf", width=16, height=8)
plot_grid(shark8444map, shark8446map, shark8448map, align="hv", ncol=3, labels=c("Shark 1: July 29-31, 2008", "Shark 2: September 09-11, 2008", "Shark 3: October 14-17, 2008"), vjust=7, hjust=-.4)
dev.off()


###### simply trying to extrapolate where interpolated points are located
s8444 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8444/8444_with_interpolated_points.csv") 
head(s8444)
interp <- s8444[s8444$interp=="TRUE",]

########################################################################################
# Next, trying to determine how a water molecule would move (total vector length) over the course of each tide
####################################################################################
# first, read in full data  
track8444 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8444_fulldata.csv")
track8446 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8446_fulldata.csv")
track8448 <- read.csv("/Volumes/McInturf/7G_Final Analysis/Shark8448_fulldata.csv")

#here, we are looking at current.x and current.y to find the current direction vectors
# we want to summarize the x and the y and add them together to find the total current vector length for that tide
# second, need to divide each full data by tide
# timestamps for 8444 tides (found by hand in the data to identify rows):
#7/29 - 14:08, 20:24 (row 96-756)
#7/30 - 03:36, 10:56, 15:09, 21:21
#7/31 - 04:23, 11:36
track8444_t1 <- track8444[1:96,] #first tidal segment
track8444_t2 <- track8444[96:756,] #second tidal segment
track8444_t3 <- track8444[757:1795,] #third tidal segment
track8444_t4 <- track8444[1796:2457,] 
track8444_t4 <- track8444[2458:2843,] # note - we are missing data here during 21:21 on the 30th and 04:23, so chose the place where the dates split (which is at approximately 11:41, so skip one full tidal cycle - in, out, in)


#timestamps for 8446 tides:
#9/9 - 13:04, 18:37
#9/10 - 01:54, 09:22, 14:39, 20:28
#9/11 - 02:36, 9:49, 13:56, 19:36
track8446_t1 <- track8446[1:69,] #first tidal segment
track8446_t2 <- track8446[70:445,] #second tidal segment (also missing chunk of data here)
track8446_t3 <- track8446[446:998,] #third tidal segment
track8446_t4 <- track8446[999:2184,] 
track8446_t5 <- track8446[2185:3487,]
track8446_t6 <- track8446[2185:3487,]
track8446_t7 <- track8446[3488:3896,]
track8446_t8 <- track8446[3897:5294,]
track8446_t9 <- track8446[5295:5980,]
track8446_t10 <- track8446[5981:6680,]
track8446_t11 <- track8446[5982:6680,]#missing chunk of data
track8446_t12 <- track8446[6680:7112,] #don't quite make it to the end of this last segment

#timestamps for 8448 tides:
#10/14 - 16:50, 23:36
#10/15 - 04:44, 11:04, 17:30
#10/16 - 00:32, 05:23, 11:40, 18:21
#10/17- 01:32, 06:06
track8448_t1 <- track8448[1:1112,] #first tidal segment
track8448_t2 <- track8448[1113:2353,] #second tidal segment
track8448_t3 <- track8448[2354:3630,]
track8448_t4 <- track8448[3631:4772,]
track8448_t5 <- track8448[4773:5703,] #missing chunk of data from 17:30-07:30
track8448_t6 <- track8448[5704:6684,] #07:30 until 11:40
track8448_t7 <- track8448[6685:8244,] 
track8448_t8 <- track8448[8245:9866,] 
track8448_t9 <- track8448[9867:10474,] #not quite at end of segment 



