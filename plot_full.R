#!/usr/bin/Rscript
source("config")

source(file.path(codeDir, "plot.R"))
source(file.path(codeDir, "prep_plot_data.R"))

require(ggplot2)

#########################################################################################
## plot the data, plots in large grid
#########################################################################################

propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

figHgt = 16
figWth = 23
figHgtIndiv = figHgt*.5
figWthIndiv = figWth*.5

pdf(file.path(outputDir, paste0('composition_fits_v', productVersion, '.pdf')), height = figHgt, width = figWth)
make_veg_map(data = pm, breaks = propBreaks, coords = coord, legendName = 'fitted proportions', map_data = usFortified, facet = TRUE)
dev.off()

pdf(file.path(outputDir, paste0('composition_fits_indiv_v', productVersion, '.pdf')), height = figHgtIndiv, width = figWthIndiv)
make_veg_map(data = pm, breaks = propBreaks, coords = coord, legendName = 'fitted proportions', map_data = usFortified, facet = FALSE)
dev.off()

psd[psd > .25] = 0.25
psdBreaks = c(0, 0.01, 0.03, 0.05, 0.075, 0.10, 0.15, 0.2, 0.25)
#psd[psd > .3] = 0.3
#psdBreaks = c(0, 0.01, 0.03, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3)

pdf(file.path(outputDir, paste0('composition_uncertainty_v', productVersion, '.pdf')), height = figHgt, width = figWth)
make_veg_map(data = psd, breaks = psdBreaks, coords = coord, legendName = 'std. error', map_data = usFortified, col = heat.colors, facet = TRUE)
dev.off()

pdf(file.path(outputDir, paste0('composition_uncertainty_indiv_v', productVersion, '.pdf')), height = figHgtIndiv, width = figWthIndiv)
make_veg_map(data = psd, breaks = psdBreaks, coords = coord, legendName = 'std. error', map_data = usFortified, col = heat.colors, facet = FALSE)
dev.off()

figHgt = 13
figWth = 20
figHgtIndiv = figHgt*.5
figWthIndiv = figWth*.59

pdf(file.path(outputDir, paste0('composition_rawData_v', productVersion, '.pdf')), height = figHgt, width = figWth)
make_areal_map(data = taxon_dat_long, variables = taxaNames, breaks = propBreaks, legendName = 'raw proportions', map_data = usFortified, facet = TRUE, legend = FALSE)
dev.off()

pdf(file.path(outputDir, paste0('composition_rawData_indiv_v', productVersion, '.pdf')), height = figHgtIndiv, width = figWthIndiv)
make_areal_map(data = taxon_dat_long, variables = taxaNames, breaks = propBreaks, legendName = 'raw proportions', map_data = usFortified, facet = FALSE, legend = TRUE)
dev.off()

