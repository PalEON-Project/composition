#!/usr/bin/Rscript

source("config")

# Jody - I'd suggest you comment out the line above (which is set up for my full workflow) and uncomment these next lines:
# productVersion <- "0.3"
# outputDir <- insert_Jody_directory_here
# 

source(file.path(codeDir, 'netCDF.R'))

require(ncdf4)

inputNcdfName <- paste0('composition_v', productVersion, '.nc')
inputNcdfPtr <- nc_open(file.path(outputDir, inputNcdfName))

pmNcdfName <- paste0('composition_v', productVersion, '_fits.nc')
psdNcdfName <- paste0('composition_v', productVersion, '_uncertainty.nc')


taxa <- names(inputNcdfPtr$var)

makeAlbersNetCDFsummary(name = 'proportion posterior mean', units = 'unitless (proportion from 0 to 1)', longname = 'posterior mean of relative composition, relative to all tree taxa,', fn = pmNcdfName, dir = outputDir, dims = inputNcdfPtr$dim[1:2], taxa = taxa)
makeAlbersNetCDFsummary(name = 'proportion posterior standard deviation', units = 'unitless (standard deviation of proportions', longname = 'posterior standard deviation of relative composition, relative to all tree taxa,', fn = psdNcdfName, dir = outputDir, dims = inputNcdfPtr$dim[1:2], taxa = taxa)

pmNcdfPtr <- nc_open(file.path(outputDir, pmNcdfName), write = TRUE)
psdNcdfPtr <- nc_open(file.path(outputDir, psdNcdfName), write = TRUE)


for(v in taxa) {
    dat <- ncvar_get(inputNcdfPtr, v, c(1, 1, 1), c(-1, -1, -1))
    pm <- apply(dat, c(1, 2), mean)
    psd <- apply(dat, c(1, 2), sd)
    # could add other summaries here, e.g. 2.5 and 97.5%iles for uncertainty intervals, using quantile()
    ncvar_put(pmNcdfPtr, v, pm, start = c(1, 1), count = c(-1, -1))
    ncvar_put(psdNcdfPtr, v, psd, start = c(1, 1), count = c(-1, -1))
}

nc_close(pmNcdfPtr)
nc_close(psdNcdfPtr)

