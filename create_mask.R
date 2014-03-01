#!/usr/bin/Rscript

# code to create netCDF file with Paleon mask for Paleon Albers grid

source('config')

require(raster)
require(ncdf4)

region = t(as.matrix(raster(file.path(dataDir, 'paleonDomain.tif'))))
water = t(as.matrix(raster(file.path(dataDir, 'water.tif'))))

region[is.na(region)] <- 0

domain <- region > 0

source(file.path(codeDir, 'set_domain.R'))

# image.plot(1:296, 1:180, region[,180:1])

x_dim <-  ncdim_def("x", "meters_east", xGrid)
y_dim <-  ncdim_def("y", "meters_north", yGrid)
vars <- list()
vars[[1]] <- ncvar_def(name = 'water', dim = list(x_dim, y_dim),
                       units = 'unitless (proportion from 0 to 1)',
                       longname = 'percentage of water in grid cell',
                       prec = 'integer')
vars[[2]] <- ncvar_def(name = 'subregion', dim = list(x_dim, y_dim),
                       units = 'unitless (region index)',
                       longname = 'PalEON subregions',
                       prec = 'integer')
vars[[3]] <- ncvar_def(name = 'domain', dim = list(x_dim, y_dim),
                       units = 'unitless (0/1 indicator)',
                       longname = 'indicator of whether grid cell is in PalEON domain',
                       prec = 'integer')

ptr <- nc_create(file.path(dataDir, 'paleonMask.nc'), vars)

ncvar_put(ptr, 'water', water[ , dim(water)[2]:1], c(1, 1), c(-1, -1))
ncvar_put(ptr, 'subregion', region[ , dim(region)[2]:1], c(1, 1), c(-1, -1))
ncvar_put(ptr, 'domain', domain[ , dim(domain)[2]:1], c(1, 1), c(-1, -1))

nc_close(ptr)
