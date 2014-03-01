#!/usr/bin/Rscript
# creates input datasets for MCMC fitting for western domain

# assumes data in western.csv
source("config")


source(file.path(codeDir, "graph.R"))

########################################################################
# read data --------------------------------------------------
########################################################################

fn <- file.path(dataDir, 'western.csv')
data <- read.csv(fn, header = TRUE, stringsAsFactors=FALSE, na.strings = "NA")
nm <- scan(pipe(paste0("head -n 1 ", fn)), sep = ',', what = 'character')
substring(nm, 1, 1) <- toupper(substring(nm, 1, 1))
names(data) <- nm

cat(paste0("Read ", nrow(data), " rows from western.csv, with field names: "))
cat(names(data), sep = ',')
cat("\n")

coordNames <- c('X', 'Y')
coord <- data[ , names(data) %in% coordNames]
coord$ID <- seq_len(nrow(coord))
coord <- coord[ , c('ID', 'X', 'Y')]

data <- data[ , !(names(data) %in% coordNames)]
data[is.na(data)] <- 0


########################################################################
# subset and manipulate taxa ------------------------------------
########################################################################

taxaUse <- scan("westernTaxaUse", what = "character", sep = ',')
taxaOther <- scan("westernTaxaOtherHardwood", what = "character", sep = ',')
taxaNonTree <- scan("westernTaxaNonTree", what = "character", sep = ',')

tmp <- c(taxaUse, taxaOther, taxaNonTree)
if(length(tmp) != length(unique(tmp)))
  warning("One or more taxa in the 'use', 'other', 'non-tree' categories overlap.")

nameConversions <- read.csv("westernTaxaNameConversions", header = TRUE, stringsAsFactors = FALSE)

otherNames <- names(data)[!(names(data) %in% c(taxaUse, taxaOther, taxaNonTree))]

if(length(otherNames)) {
  cat("Found these additional taxa categories, adding them to 'other':")
  cat(otherNames)
}

taxaOther <- c(taxaOther, otherNames)

other <- rowSums(data[ , taxaOther])
data <- data[ , taxaUse]

cat("Grouping the following taxa into an 'Other hardwood' category:")
cat(taxaOther, sep = ',')
cat("\n")

numConv <- nrow(nameConversions)
if(numConv) {
  for(i in seq_len(numConv)) {
    if(!(nameConversions$to[i] %in% taxaUse) || !(nameConversions$from[i] %in% taxaUse)) {
      warning(paste0("one or more names in name conversions file not in file of taxa to use:", nameConversions$to[i], " ", nameConversions$from[i]))
    } else {
      data[nameConversions$to[i]] <- data[nameConversions$to[i]] + data[nameConversions$from[i]]
    }
  }
  data <- data[ , !(names(data) %in% nameConversions$from)]
  taxaUse <- taxaUse[!(taxaUse%in% nameConversions$from)]
}

ord <- order(names(data))

data <- data[ , ord]
data$"Other hardwood" <- other

nTaxa <- ncol(data)
taxaUse <- names(data)

taxa <- data.frame(taxonID = 1:nTaxa, taxonName = taxaUse, stringsAsFactors = FALSE)

cat("Using the following", nTaxa, "taxa, with", sum(data), "trees, of which", round(100*sum(data$Other)/sum(data), 2), "% are in the 'other' category:")
cat(taxaUse, sep = ",")
cat("\n")

########################################################################
# subset to portion of PalEON Albers grid and create graph -------------
########################################################################

source(file.path(codeDir, 'set_domain.R'))

data <- data[ coord$X <= easternLimitOfWesternDomain, ]
coord <- coord[ coord$X <= easternLimitOfWesternDomain, ]
  
data <- data[order(coord$ID), ]
coord <- coord[order(coord$ID), ]
dimnames(coord)[[1]] <- coord$ID

m1 <- length(westernDomainX)
m2 <- length(westernDomainY)
nCells <- m1 * m2

tmp <- nbhdStructure
substring(tmp, 1 ,1) <- toupper(substring(tmp, 1, 1))
graphFileName <- paste0('graph', tmp, m1, 'x', m2, '.csv')
if(!file.exists(file.path(dataDir, graphFileName)))
  graphCreate(m1, m2, type = nbhdStructure, dir = dataDir, fn = graphFileName)
nbhd <- graphRead(file.path(dataDir, graphFileName), m1, m2)


########################################################################
# create data objects for MCMC fitting ---------------------------------
########################################################################


data <- as.matrix(data)

taxon <- rep(rep(1:nTaxa, nCells), times = c(t(data)))
total <- rowSums(data)
cell <- rep(which(total > 0), times = total[total > 0])

data <- data.frame(taxon = taxon, cell = cell)

save(nbhd, m1, m2, nTaxa, nCells, data, coord, taxa, file = paste0(file.path(dataDir, 'westernData.Rda')))
