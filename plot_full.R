#!/usr/bin/Rscript

source("config")

source(file.path(codeDir, "plot.R"))
source(file.path(codeDir, "set_domain.R"))

require(ncdf4)

require(reshape) # reshape2?
require(ggplot2)
require(gridExtra)
require(maptools)
require(gpclib)
require(sp)
require(rgdal)
require(fields)
require(raster)
require(maps)
require(ColorBrewer)

usShp <- readShapeLines(file.path(dataDir, 'us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')




########################################################################################
## get raw data in form to use
########################################################################################

# don't I need rev(easternDomainY)?
region = t(as.matrix(raster(file.path(dataDir, 'paleonDomain.tif'))))
water = t(as.matrix(raster(file.path(dataDir, 'water.tif'))))
# t()  manipulates matrix so plots correctly W-E and N-S in R

# region[region %in% c(2,3,5,6,11,12)] <- NA
water[water == 100] <- NA
mask = is.na(region)
maskWater = is.na(water)

# western data/results

load(file.path(dataDir, 'westernData.Rda'))
coordsWest <- coord

rawWest <- matrix(0, nrow = nCells, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$cell[data$taxon == p])
  rawWest[as.numeric(names(tbl)) , p] <- tbl
}
total <- rowSums(rawWest)
rawWest[total == 0] = NA
total[total == 0] <- 1
rawWest <- rawWest / total
dimnames(rawWest)[[2]] <- gsub("/", "ZZZ", taxa$taxonName)  # otherwise both / and " " become "." so can't distinguish when I substitute back in for "."

finalNcdfName <- paste0('PLScomposition_western_', productVersion, '-release.nc')

ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))
test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nSamples <- dim(test)[3]

if(nCells != prod(dim(test)[1:2]))
  stop("nCells does not match first dimension of netCDF file.")

preds <- array(0, c(nCells, nTaxa, nSamples))
#dimnames(preds)[[2]] <- taxa$taxonName
for(p in 1:nTaxa) {
  preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1))
}

attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

pmWest <- apply(preds, c(1, 2), 'mean')
psdWest <- apply(preds, c(1, 2), 'sd')



# eastern data/results

load(file.path(dataDir, 'easternData.Rda'))
load(file.path(dataDir, 'intersection.Rda'))

  

finalNcdfName <- paste0('PLScomposition_eastern_', productVersion, '-release.nc')

ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))
test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nSamples <- dim(test)[3]

nCells <- m1*m2

if(nCells != prod(dim(test)[1:2]))
  stop("nCells does not match first dimension of netCDF file.")


preds <- array(0, c(nCells, nTaxa, nSamples))
for(p in 1:nTaxa) {
  preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1))
}

attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

pmEast <- apply(preds, c(1, 2), 'mean')
psdEast <- apply(preds, c(1, 2), 'sd')

## combine across grids


pmFull <- matrix(NA, c(xRes*yRes), ncol(pmEast))
dimnames(pmFull)[[2]] <- dimnames(pmEast)[[2]]
fullTmp <- psdFull <- rawFull <- pmFull

ids <- matrix(1:(xRes*yRes), nrow = xRes, ncol = yRes)
eastSubset <- c(ids[easternDomainX, easternDomainY])
westSubset <- c(ids[westernDomainX, westernDomainY])

eastOnly <- dimnames(pmEast)[[2]]
eastOnly <- eastOnly[!(eastOnly %in% dimnames(pmWest)[[2]])]

pmFull[eastSubset, ] <- pmEast
fullTmp[westSubset, dimnames(pmWest)[[2]]] <- pmWest
fullTmp[!(region %in% c(6,12,11, 5, 2, 3)), ] <- NA
pmFull[!is.na(fullTmp)] <- fullTmp[!is.na(fullTmp)]
pmFull[region %in% c(5, 12), eastOnly] <- NA

psdFull[eastSubset, ] <- psdEast
fullTmp[westSubset, dimnames(psdWest)[[2]]] <- psdWest
fullTmp[!(region %in% c(6,12,11, 5, 2, 3)), ] <- NA
psdFull[!is.na(fullTmp)] <- fullTmp[!is.na(fullTmp)]
psdFull[region %in% c(5, 12), eastOnly] <- NA

rawFull[westSubset, dimnames(rawWest)[[2]]] <- rawWest

pmFull[mask, ] <- NA
psdFull[mask, ] <- NA
rawFull[mask, ] <- NA

#tmp <- pmFull
# tmp<tmp[1:296, 180:1]
#image.plot(1:296,1:180, matrix(tmp,296,180))

coordFull <- expand.grid(X = xGrid, Y = rev(yGrid))

#########################################################################################
## plot the data, plots in large grid
#########################################################################################

propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

figHgt = 16
figWth = 22
figHgtIndiv = 16*.5
figWthIndiv = 22*.5

pdf(file.path(outputDir, paste0('PLScomposition_full_', productVersion, '_fits.pdf')), height = figHgt, width = figWth)
make_veg_map(data = pmFull, breaks = propBreaks, coords = coordFull, legendName = 'fitted proportions', map_data = usFortified, facet = TRUE)
dev.off()

pdf(file.path(outputDir, paste0('PLScomposition_full_', productVersion, '_fits_indiv.pdf')), height = figHgtIndiv, width = figWthIndiv)
make_veg_map(data = pmFull, breaks = propBreaks, coords = coordFull, legendName = 'fitted proportions', map_data = usFortified, facet = FALSE)
dev.off()

psdFull[psdFull > .25] = 0.25
psdBreaks = c(0, 0.01, 0.03, 0.05, 0.075, 0.10, 0.15, 0.2, 0.25)
#psd[psd > .3] = 0.3
#psdBreaks = c(0, 0.01, 0.03, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3)
  
pdf(file.path(outputDir, paste0('PLScomposition_full_', productVersion, '_uncertainty.pdf')), height = figHgt, width = figWth)
make_veg_map(data = psdFull, breaks = psdBreaks, coords = coordFull, legendName = 'std. error', map_data = usFortified, col = heat.colors, facet = TRUE)
dev.off()

pdf(file.path(outputDir, paste0('PLScomposition_full_', productVersion, '_uncertainty_indiv.pdf')), height = figHgtIndiv, width = figWthIndiv)
make_veg_map(data = psdFull, breaks = psdBreaks, coords = coordFull, legendName = 'std. error', map_data = usFortified, col = heat.colors, facet = FALSE)
dev.off()


# this is not working and not sure exactly how to combine the raster info in West with polygon info in east - probably convert rasters to polygons and then combine the polygons
if(FALSE) { 
# plot raw data as colored polygons
# it would be nice to have this as a function, but there is a lot of pre-processing...

# west

  data_binned <- matrix(0, nrow = nrow(rawWest), ncol = ncol(rawWest))
  colnames(data_binned) <- colnames(rawWest)
  
  for (p in seq_len(ncol(data))) {
    data_binned[ , p] <- cut(rawWest[ , p], breaks, include.lowest = TRUE, labels = FALSE)
  }
    taxon_dat <- data.frame(X = coordsWest$X, Y = coordsWest$Y, data_binned)
    taxon_dat_long_west <- melt(taxon_dat, c('X', 'Y'))

# east


nTowns <- dim(townCellOverlap)[1]
easternDataDir <- paste0(easternVersionID, '.', easternVersion)
ohioDataDir <- paste0(ohioVersionID, '.', ohioVersion)
eastern_townships <- readOGR(file.path(dataDir, easternDataDir), paste0(easternVersionID, 'polygonsver', easternVersion))
ohio_townships <- readOGR(file.path(dataDir, ohioDataDir), paste0('OH', ohioVersionID, 'polygonsver', ohioVersion))
ohio_townships <- spTransform(ohio_townships, CRSobj=CRS('+init=epsg:3175'))  # transform to Albers

raw <- matrix(0, nrow = nTowns, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$town[data$taxon == p])
  raw[as.numeric(names(tbl)) , p] <- tbl
}
raw <- raw / rowSums(raw)
attributes(raw)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

for(p in seq_len(nTaxa))
  raw[ , p] <- cut(raw[ , p], propBreaks, include.lowest = TRUE, labels = FALSE)

east_fort <- fortify(eastern_townships)
ohio_fort <- fortify(ohio_townships)

idMap = c(eastern_townships@data$ID, ohio_townships@data$ID + 1 + length(eastern_townships))
# this is dangerous as it relies on data$town corresponding to $id+1 values where $id is from the @polygons 
east_fort$id = as.numeric(east_fort$id) + 1
ohio_fort$id = as.numeric(ohio_fort$id) + 1 + length(eastern_townships)
fort <- rbind(east_fort, ohio_fort)
fort$id = idMap[fort$id]
fort <- cbind(fort, raw[fort$id,])
names(fort)[1:2] <- c('X', 'Y')

#ggplot(fort, aes(long, lat, group = id)) + geom_polygon(aes(fill = Oak))

names(fort) <- gsub("ZZZ", "/", names(fort))
taxon_dat_long <- melt(fort, c('X', 'Y', 'order', 'hole', 'piece', 'group', 'id' ))

breaks = propBreaks
breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                       function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

legendName = "raw proportions"

# make_veg_map(data = rawFull, breaks = propBreaks, coords = coord, map_data = usFortified, legendName = 'raw proportions', facet = TRUE)

# ggplot(taxon_dat_long, aes(X, Y, group = id)) + geom_polygon(aes(fill = as.factor(value)))

d <- ggplot(taxon_dat_long, aes(X, Y, group = id)) +
  geom_polygon(aes(fill = as.factor(value))) + 
  scale_fill_manual(values = tim.colors(length(breaks)), labels = breaklabels, name = legendName) +
  theme(strip.text.x = element_text(size = 16), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  coord_fixed() +
  geom_raster(data = taxon_dat_long_west, aes(x = X, y = Y, fill = factor(value))) +
  facet_wrap( ~ variable, ncol=5)
d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = taxon_dat_long)
d <- theme_clean(d)

pdf(file.path(outputDir, paste0('PLScomposition_eastern_', productVersion, '_rawData.pdf')), height = figHgt, width = figWth)
print(d)
dev.off()

}