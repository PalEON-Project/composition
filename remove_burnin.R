#!/usr/bin/Rscript
# subsets to post-burnin samples

# also flips N-S so that values are increasing from S to N
# interim netCDF files go from N to S becaus that is how the tif files
# are set up

source("config")

# grab from input args after sourcing config as we want to overwrite burnin
args <- commandArgs(TRUE)
# now args is a character vector containing the arguments
burnin <- as.numeric(args[1])
domain <- args[2]

runID <- paste0(domain, "_", runID)

#source("tmp.config")
source(file.path(codeDir, 'netCDF.R'))
source(file.path(codeDir, 'set_domain.R'))
require(ncdf4)

load(file.path(dataDir, paste0('data_', runID, '.Rda')))

outputNcdfName <- paste0('PLScomposition_', runID, '_full.nc')

finalNcdfName <- paste0('PLScomposition_', runID, '.nc')

outputNcdfPtr <- nc_open(file.path(outputDir, outputNcdfName))

if(domain == "western") 
    makeAlbersNetCDF(name = 'proportion', units = 'unitless (proportion from 0 to 1)', longname = 'relative composition, relative to all tree taxa,', fn = finalNcdfName, dir = outputDir, x = sort(unique(coord$X)), y = sort(unique(coord$Y)), taxa = taxa$taxonName, numSamples = floor((S-burnin)/(thin*secondThin)))
if(domain == "eastern")
    makeAlbersNetCDF(name = 'proportion', units = 'unitless (proportion from 0 to 1)', longname = 'relative composition, relative to all tree taxa,', fn = finalNcdfName, dir = outputDir, x = xGrid[easternDomainX], y = yGrid[easternDomainY], taxa = taxa$taxonName, numSamples = floor((S-burnin)/(thin*secondThin)))

finalNcdfPtr <- nc_open(file.path(outputDir, finalNcdfName), write = TRUE)

currentSamples <- floor(S/(thin*secondThin))
finalSamples <- floor((S-burnin)/(thin*secondThin))
for(p in seq_len(nTaxa)) {
    tmp <- ncvar_get(outputNcdfPtr, varid = taxa$taxonName[p], start = c(1, 1, (currentSamples - finalSamples + 1)), count = c(-1, -1, finalSamples))
    # m2:1 flips the y-axis so that y values go from S to N
    ncvar_put(finalNcdfPtr, taxa$taxonName[p], tmp[ , m2:1, ], start = c(1, 1, 1), count = c(-1, -1, finalSamples))
}

nc_close(outputNcdfPtr)
nc_close(finalNcdfPtr)

