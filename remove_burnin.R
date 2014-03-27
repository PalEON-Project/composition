#!/usr/bin/Rscript
# subsets to post-burnin samples 
source("config")
source("tmp.config")
source(file.path(codeDir, 'netCDF.R'))
source(file.path(codeDir, 'set_domain.R'))
require(ncdf4)

if(!exists('uniqueRunID'))
  uniqueRunID <- ""
if(uniqueRunID == "")
  fnAdd <- "" else fnAdd <- paste0("-run", uniqueRunID)


  load(file.path(dataDir, paste0(domain, 'Data.Rda')))


  outputNcdfName <- paste0('PLScomposition_', domain, '_', productVersion, fnAdd, '.nc')

  finalNcdfName <- paste0('PLScomposition_', domain, '_', productVersion, '-release.nc')

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
    ncvar_put(finalNcdfPtr, taxa$taxonName[p], tmp, start = c(1, 1, 1), count = c(-1, -1, finalSamples))
  }
  
  nc_close(outputNcdfPtr)
  nc_close(finalNcdfPtr)

