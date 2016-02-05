#!/usr/bin/Rscript

source("config")

# grab from input args after sourcing config as we want to overwrite burnin
args <- commandArgs(TRUE)
# now args is a character vector containing the arguments
burnin <- as.numeric(args[1])


require(ncdf4)

source(file.path(codeDir, "netCDF.R"))
source(file.path(codeDir, "set_domain.R"))

stitch <- function(east, west, m1, m2, numS) {
    full <- as.numeric(NA); length(full) <- m1*m2*numS; dim(full) <- c(m1, m2, numS)
    for(s in seq_len(numS)) {
        tmpe <- as.numeric(NA); length(tmpe) <- m1*m2; dim(tmpe) <- c(m1, m2)
        tmpw <- as.numeric(NA); length(tmpw) <- m1*m2; dim(tmpw) <- c(m1, m2)
        
        # the code in the y-dimension reverses the fact that
        # in the stat modeling the y-dimensions started with 1 in the N to 180 in S
        # remove_burnin.R flipped the y-axis but easternDomainY is not flipped
        # so need to rev() and then put in correct position in full grid
        # 1-17 empty, 18-157 with data, 158-180 empty
        tmpe[easternDomainX, m2+1-rev(easternDomainY)] <- east[ , , s]
        tmpw[1:m1_west, ] <- west[ , , s]

        tmpe[!regions %in% east_regions] <- 0
        tmpw[!regions %in% west_regions] <- 0

        tmp <- (tmpe + tmpw) 

        tmp[regions == 0] <- NA
        full[ , , s] <-  tmp
    }
    return(full)
}



outputNcdfName <- paste0("composition_v", productVersion, ".nc")

load(file.path(dataDir, paste0('data_western_', runID, '.Rda')))
taxa_west <- taxa
m1_west <- m1; m2_west <- m2
load(file.path(dataDir, paste0('data_eastern_', runID, '.Rda')))
taxa_east <- taxa
m1_east <- m1; m2_east <- m2

taxaNames <- unique(taxa_east$taxonName, taxa_west$taxonName)

setwd(outputDir)

numSamples <- floor((S-burnin)/(thin*secondThin))

makeAlbersNetCDF(name = 'proportion', longname = 'relative composition, relative to all tree taxa,', fn = outputNcdfName, dir = outputDir, x = xGrid, y = yGrid, taxa = taxaNames, numSamples = numSamples)

outputNcdfPtr <- nc_open(file.path(outputDir, outputNcdfName), write = TRUE)

inputNcdfPtr_west <- nc_open(file.path(outputDir, paste0('PLScomposition_western_', runID, '.nc')))
inputNcdfPtr_east <- nc_open(file.path(outputDir, paste0('PLScomposition_eastern_', runID, '.nc')))

mask <- nc_open(file.path(dataDir, 'paleonMask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1),c(-1,-1))
wat <- ncvar_get(mask, 'water', c(1,1),c(-1,-1))
dom <- ncvar_get(mask, 'domain', c(1,1),c(-1,-1))

west_regions <- c(2,3,5,6,11,12 )
east_regions <- c(1,4,7,8,9,10)

for(p in seq_along(taxaNames)) {
    east <- try(ncvar_get(inputNcdfPtr_east, taxaNames[p], c(1,1,1), c(-1,-1,-1)))
    west <- try(ncvar_get(inputNcdfPtr_west, taxaNames[p], c(1,1,1), c(-1,-1,-1)))
    if(class(east) == 'try-error') {
        east <- array(as.numeric(NA), c(m1_east, m2_east, numSamples))
        cat("Filling NA for ", taxaNames[p], " in east.\n")
    }
    if(class(west) == 'try-error') {
        west <- array(as.numeric(NA), c(m1_west, m2_west, numSamples))
        cat("Filling NA for ", taxaNames[p], " in west\n")
    }
    tmp <- stitch(east, west, xRes, yRes, numSamples)
    ncvar_put(outputNcdfPtr, taxaNames[p], tmp, start = c(1, 1, 1), count = c(xRes, yRes, numSamples))
}

nc_close(outputNcdfPtr)
nc_close(inputNcdfPtr_west)
nc_close(inputNcdfPtr_east)

    
