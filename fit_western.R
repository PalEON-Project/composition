#!/usr/bin/Rscript

# assumes inputs in westernData.Rda

source("config")

source(file.path(codeDir, "mcmc.R"))
source(file.path(codeDir, "netCDF.R"))
source(file.path(codeDir, "set_domain.R"))


if(!exists('uniqueRunID'))
  uniqueRunID <- ""

# fit model --------------------------------------------------

load(file.path(dataDir, 'westernData.Rda'))

if(uniqueRunID == "")
  fnAdd <- "" else fnAdd <- paste0("-run", uniqueRunID)

latentNcdfName <- paste0('PLScomposition_raw_western_', productVersion, fnAdd, '.nc')

if(!resumeRun) {
  set.seed(0)
  makeAlbersNetCDF(name = 'latent', units = 'unitless', longname = 'latent multivariate logit values', fn = latentNcdfName, dir = dataDir, x = xGrid[westernDomainX], y = yGrid[westernDomainY], taxa = taxa$taxonName, numSamples = floor(S/thin))
}


# this creates netCDF with draws of the latent variables
# change to signature of new runMCMC that allows township data
out <- runMCMC(y = data$taxon, cell = data$cell, C = nbhd, Cindices = nbhdIndices, S = S, thin = thin, resumeRun = resumeRun, hyperpar = c(-0.5, 0), areallyAggregated = FALSE, outputNcdfName = latentNcdfName, taxa = taxa, numCores = numCoresToUse, jointSample = jointSample, adaptStartDecay = 5000, runID = paste0("-western-run", uniqueRunID), dataDir = dataDir, outputDir = outputDir)


# post process to get draws of proportions
outputNcdfName <- paste0('PLScomposition_western_', productVersion, fnAdd, '.nc')

makeAlbersNetCDF(name = 'proportion', units = 'unitless (proportion from 0 to 1)', longname = 'relative composition, relative to all tree taxa,', fn = outputNcdfName, dir = outputDir, x = sort(unique(coord$X)), y = sort(unique(coord$Y)), taxa = taxa$taxonName, numSamples = floor(S/(thin*secondThin)))

latentNcdfPtr <- nc_open(file.path(dataDir, latentNcdfName))
outputNcdfPtr <- nc_open(file.path(outputDir, outputNcdfName), write = TRUE)

# this draws the proportions based on the draws of the latent variables
set.seed(0)
drawProportions(latentNcdfPtr, outputNcdfPtr, numMCsamples = numSamplesForProps, numInputSamples = floor(S/thin), secondThin = secondThin, I = m1*m2, taxa = taxa$taxonName)

nc_close(latentNcdfPtr)
nc_close(outputNcdfPtr)
