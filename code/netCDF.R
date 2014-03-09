# auxiliary code for creating a netCDF file and returning a pointer to it
require(ncdf4)

makeAlbersNetCDF <- function(name = NULL, units = '', longname = '', fn = NULL, dir = '', x, y, taxa, numSamples) {
  if(is.null(fn)) stop("makeAlbersNetCDF: requires filename as argument.")

  x_dim <-  ncdim_def("x", "meters_east", x)
  y_dim <-  ncdim_def("y", "meters_north", y)
  draw_dim <- ncdim_def("sample", "number", 1:numSamples, longname = "MCMC sample")
  vars <- list()
  length(vars) == length(taxa)
  for(k in seq_along(taxa)) {
    vars[[k]] <- ncvar_def(name = taxa[k], dim=list(x_dim, y_dim, draw_dim), units = units, missval=1e20, longname = paste0(longname, " for taxon ", taxa[k]), prec="double")
  }
  ncdfPtr <- nc_create(file.path(dir, fn), vars)
  nc_close(ncdfPtr)
  invisible(NULL)
}