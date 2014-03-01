#!/usr/bin/Rscript

source("config")

source(file.path(codeDir, "plot.R"))
source(file.path(codeDir, "set_domain.R"))

library(ncdf4)

library(reshape) # reshape2?
library(ggplot2)
library(gridExtra)
library(maptools)
library(gpclib)
library(sp)
library(rgdal)
library(fields)
library(raster)
 
 

usShp <- readShapeLines(file.path(dataDir, 'us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')


########################################################################################
## get raw data and model output in form to use
########################################################################################

# don't I need rev(westernDomainY)?
region = t(as.matrix(raster(file.path(dataDir, 'paleonDomain.tif'))))[westernDomainX, westernDomainY]
region[region == 9] <- NA # Ohio
water = t(as.matrix(raster(file.path(dataDir, 'water.tif'))))[westernDomainX, westernDomainY]
# t() manipulates matrix so plots correctly W-E and N-S in R

mask = is.na(region)


load(file.path(dataDir, 'westernData.Rda'))

raw <- matrix(0, nrow = nCells, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$cell[data$taxon == p])
  raw[as.numeric(names(tbl)) , p] <- tbl
}
total <- rowSums(raw)
raw[total == 0] = NA
total[total == 0] <- 1
raw <- raw / total
dimnames(raw)[[2]] <- gsub("/", "ZZZ", taxa$taxonName)  # otherwise both / and " " become "." so can't distinguish when I substitute back in for "."

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

pm <- apply(preds, c(1, 2), 'mean')
psd <- apply(preds, c(1, 2), 'sd')

pm[mask, ] <- NA
psd[mask, ] <- NA
raw[mask, ] <- NA

#########################################################################################
## make plots
#########################################################################################

propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

figHgt = 20
figWth = 20

pdf(file.path(outputDir, paste0('PLScomposition_western_', productVersion, '_fits.pdf')), height = figHgt, width = figWth)
make_veg_map(data = pm, breaks = propBreaks, coords = coord, legendName = 'fitted proportions', map_data = usFortified, facet = TRUE)
dev.off()

psd[psd > .25] = 0.25
psdBreaks = c(0, 0.01, 0.03, 0.05, 0.075, 0.10, 0.15, 0.2, 0.25)
  
pdf(file.path(outputDir, paste0('PLScomposition_western_', productVersion, '_uncertainty.pdf')), height = figHgt, width = figWth)
make_veg_map(data = psd, breaks = psdBreaks, coords = coord, legendName = 'std. error', map_data = usFortified, facet = TRUE)
dev.off()
  

pdf(file.path(outputDir, paste0('PLScomposition_western_', productVersion, '_rawData.pdf')), height = figHgt, width = figWth)
make_veg_map(data = raw, breaks = propBreaks, coords = coord, map_data = usFortified, legendName = 'raw proportions', facet = TRUE)
dev.off()


######################################################################
### mixing
######################################################################


if(FALSE) {
  if(!exists('uniqueRunID'))
    uniqueRunID <- ""
  if(uniqueRunID == "")
    fnAdd <- "" else fnAdd <- paste0("-run", uniqueRunID)

  load(file.path(outputDir, paste0('sigma2-western', fnAdd, '.Rda')))
  par(mfrow = c(5, 5), mai = c(.3, .3, .1, .1))
  for(p in 1:ncol(sigma2store))
    tsplot(sigma2store[ , p])
}


