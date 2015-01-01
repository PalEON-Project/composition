#!/usr/bin/Rscript

# code to intersect Charlie's townships with Simon's 8 km grid
# result is 180x296x(1372+471) array of proportion intersection
# where 180 is y dimension and 296 is x dimension
# open issue: intersection only uses discrete approximation with 100 points


#require(rgeos)
#require(maptools)
require(rgdal)
require(raster)

source("config")

easternDataDir <- "eastern"
ohioDataDir <- "ohio"

####################################################################
# read in shape file info and create raster for PalEON Albers grid
####################################################################

eastern_townships <- readOGR(file.path(dataDir, easternDataDir), paste0(easternVersionID, 'polygonsver', easternVersion))
ohio_townships <- readOGR(file.path(dataDir, ohioDataDir), paste0('OH', ohioVersionID, 'polygonsver', ohioVersion))

#proj4string(ohio_townships) <- CRS('+init=epsg:4326')  # seems to have a lat/lon proj now, so don't need this
ohio_townships <- spTransform(ohio_townships, CRSobj=CRS('+init=epsg:3175'))  # transform to Albers

nTowns <- length(eastern_townships) + length(ohio_townships)

source(file.path(codeDir, "set_domain.R"))

rast <- raster(crs = CRS('+init=epsg:3175'),
               xmn = xRange[1], xmx = xRange[2],
               ymn = yRange[1], ymx = yRange[2],
               ncols = xRes, nrows = yRes)

# this raster has rows as y and columns as x and starts with 1,1 in the NW corner, so 180,1 is the SW corner

# NOTE: do not subset as townships@polygons[[index]] as the sorting of this is different than townships[index, ]


####################################################################
# intersect grid with townships ------------------------------------
####################################################################

mini <- min(unique(eastern_townships$ID))

# intersect with eastern townships
for(i in sort(unique((eastern_townships$ID)))){
  aa <- rasterize(x = eastern_townships[eastern_townships$ID == i, ], y = rast, getCover = TRUE)  
  if(i == mini){
    poly.stack <- stack((aa)) 
  }
  else{
    poly.stack <- addLayer(poly.stack, aa)
  }
}

# intersect with Ohio townships
for(i in sort(unique((ohio_townships$ID)))){
  aa <- rasterize(x = ohio_townships[ohio_townships$ID == i, ], y = rast, getCover = TRUE)
  poly.stack <- addLayer(poly.stack, aa)
}

interTmp <- as.array(poly.stack) * (64/100)  # 100 converts to proportion; 64 has to do with 8x8?

# check
if(FALSE) {
  area <- unlist(sapply(eastern_townships@polygons, function(x) x@area/1000000))
  plot(area, apply(interTmp, 3, sum))
}

inter <- array(0, c(xRes, yRes, nTowns))
for(i in 1:nTowns)
  inter[ , , i] <- t(interTmp[ , , i]/sum(interTmp[ , , i]))


# inter goes NW to NE and proceeds in lat bands southward. Last cell is the SE corner
# for plotting, I'll need to reverse the columns

nCells <- xRes*yRes
ids <- 1:nCells
usedIds <- c(matrix(ids, xRes, yRes)[easternDomainX, easternDomainY])

inter <- inter[easternDomainX, easternDomainY, ]

save(inter, usedIds, nTowns, file = file.path(dataDir, paste0('intersection_eastern_', runID, '.Rda')))


