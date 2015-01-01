#!/usr/bin/Rscript

# assumes inputs in easternData.Rda

source("config")

runID <- paste0("eastern_", runID)

source(file.path(codeDir, "mcmc.R"))  
source(file.path(codeDir, "netCDF.R"))
source(file.path(codeDir, "set_domain.R"))

library(RhpcBLASctl)

omp_set_num_threads(1)

tmp <- chol(as.spam(diag(rep(1,3)))) # chol.spam failing in runMCMC for some reason unless do this as a sort of initialization step

# fit model --------------------------------------------------

load(file.path(dataDir, paste0('data_', runID, '.Rda')))


latentNcdfName <- paste0('PLScomposition_raw_', runID, '.nc')

if(!resumeRun) {
  set.seed(seed)
  makeAlbersNetCDF(name = 'latent', units = 'unitless', longname = 'latent multivariate logit values', fn = latentNcdfName, dir = dataDir, x = xGrid[easternDomainX], y = yGrid[easternDomainY], taxa = taxa$taxonName, numSamples = floor(S/thin))
}

# this creates netCDF with draws of the latent variables
out = runMCMC(y = data$taxon, cell = NULL, C = nbhd, town = data$town,
  townCellOverlap = townCellOverlap, townCellIds = townCellIds,
  S = S, thin = thin, resumeRun = resumeRun, hyperpar = c(-0.5, 0),
  nbhdStructure = nbhdStructure,
  areallyAggregated = TRUE, outputNcdfName = latentNcdfName, taxa = taxa,
  numCores = numCoresToUse, runID = runID, dataDir = dataDir,
  outputDir = outputDir)


# post process to get draws of proportions
outputNcdfName <- paste0('PLScomposition_', runID, '_full.nc')
makeAlbersNetCDF(name = 'proportion', units = 'unitless (proportion from 0 to 1)',
                 longname = 'relative composition, relative to all tree taxa,',
                 fn = outputNcdfName, dir = outputDir, x = xGrid[easternDomainX],
                 y = yGrid[easternDomainY], taxa = taxa$taxonName,
                 numSamples = floor(S/(thin*secondThin)))

latentNcdfPtr <- nc_open(file.path(dataDir, latentNcdfName))
outputNcdfPtr <- nc_open(file.path(outputDir, outputNcdfName), write = TRUE)

# this draws the proportions based on the draws of the latent variables
set.seed(seed)
drawProportions(latentNcdfPtr, outputNcdfPtr, numMCsamples = numSamplesForProps,
                numInputSamples = floor(S/thin), secondThin = secondThin,
                I = m1*m2, taxa = taxa$taxonName)

warning("remember to augment eastern domain to include the areas with no info, or just warn Jackie")

nc_close(latentNcdfPtr)
nc_close(outputNcdfPtr)

