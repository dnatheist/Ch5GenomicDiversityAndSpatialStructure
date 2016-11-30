## copied script



#maps in R

#set your working directory to a writable folder (ALA4R stores its maps there)
setwd("c:/temp")


#ALA4R package is the Atlas for living Australia package for R, here we get our sample data from
require(ALA4R)

#if this procudes an error then the package is not installed yet to install the pacakge uncomment the following code and run it

#install.packages(c("httr","stringr","plyr","digest","RCurl","jsonlite","assertthat","sp"))
#install.packages(c("data.table"))

#


#this downloads data from the Atlas server (you can change the species name to your favourite one)
x <- occurrences(taxon="gehyra variegata",download_reason_id=10) 

#creates a plot, check your working directory for a pdf file afterwards
occurrences_plot(x)


#look at the x object....
#create some plots...
plot( x$data$longitude,x$data$latitude)


#########now our own first map
#we use the image provided from the package
image(aus, col = "grey")
#and put our points manually on the map
points(x$data$longitude,x$data$latitude)


#dismo package
require(dismo)
#or you not installed you need to install the package:
#install.packages("dismo")


#### this line, gets a map from the google server. Check ?gmap to see how you can download different versions (e.g. satellite map, street map etc.)
g = gmap('Australia' )
plot(g)

#the next lines of codes show you how to reproject our data from the Atlas. 
#because Google needs them in the Mercator projection and Atlas gives them as latlong
# reprojections need to be done quite often when you have different data sources on a map

# transform your coordinates....
#get the coordinates from the Atlas
long <- x$data$longitude
lat <- x$data$latitude
#and store them in a new data.frame
latlong <- cbind(long, lat)

#here we do the reprojection into the Google Mercator format
mercxy <- Mercator(latlong)
#this plots the reprojected points on the map
points(mercxy)


#to customise our plot we can interactively select the extent with the mouse
#find the right extent....
#plot the Australia map
plot(g)
#once you send this command you neet to click on the map and select the upper-left and lower right corner of the area you want the map to crop to
e<- drawExtent()

#plot the croped map
ge <- gmap(e)
plot(ge)
points(mercxy, pch=16, col=rgb(1,0,0,0.2))


####excursion: do some spatial statistics
#this part is only to demonstrate that it is quite easy to do some spatial statistics with R
#e.g. we put a grid on the map and count the number of points within a grid cell
#propbably you need to install the package via

#install.packages("spatstat")
require(spatstat)

#remove missing values (as Atlas Data are sometimes messy)
xy <- mercxy[!is.na(mercxy[,1]),]
#define the area
pp <-ppp(xy [,1], xy[,2], window=owin(c(e@xmin, e@xmax), c(e@ymin, e@ymax)), unitname=c("metre","metres"))

#here we count the number of points in a 10 x 10 grid
qc <-quadratcount(pp,nx=10, ny=10)

plot(ge)
points(mercxy, pch=16, col=rgb(1,0,0,0.2))
plot(qc, add=TRUE)



#ggmap
#another package to do maps 
#install.packages("ggmap") # might be necessary
require(ggmap)

map <- get_map("australia", zoom=4)
ggmap(map)
latlongd <- data.frame(latlong)


ggmap(map, extent = 'device') +
  geom_point(aes(x = long, y = lat), data = latlongd, alpha = .5)

# contourlines
ggmap(map) +
  geom_point(aes(x = long, y = lat ), data = latlongd, alpha = .5, color="orange") + 
  stat_density2d(aes(x = long, y = lat), data = latlongd, bins = 5, lwd=2) 



#and some other maps to play with....
map <- get_map("university of canberra, australia", zoom=18, maptype="satellite")
ggmap(map) 

#to annotate click on the map at the intended point
(clicks <- gglocator(1) )

ggmap(map) +
  annotate('text', x=clicks[1,1], y=clicks[1,2], label = 'Building 2', colour = I('green'), size = 6)
