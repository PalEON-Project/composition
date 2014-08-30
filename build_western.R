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

taxaInfo <- read.csv(paste0('level3s_v', productVersion, '.csv'),
                     stringsAsFactors = FALSE)
taxaInfo[ , ncol(taxaInfo)] <- gsub("\\s", "", taxaInfo[ , ncol(taxaInfo)])  # strip any (trailing) whitespace

taxaOtherHardwood <- taxaInfo[["Level.3a"]][taxaInfo[["Level.3s"]] == "Other hardwood"]
taxaUse <- unique(taxaInfo[["Level.3s"]][!(taxaInfo[["Level.3s"]] %in% c("Other hardwood", "None")) & taxaInfo[["omit.western"]] == "no"])
taxaNonTree <- taxaInfo[["Level.3a"]][taxaInfo[["Level.3s"]] == "None"]

otherNames <- unique(names(data))[!(names(data) %in% taxaInfo[["Level.3a"]])]
if(length(otherNames)) {
  cat("Warning: Found these additional taxa categories, ignoring them:")
  cat(otherNames, sep = ',')
}

cat("Grouping the following taxa into an 'Other hardwood' category:")
cat(taxaOtherHardwood, sep = ',')
cat("\n")

nameConversions <- taxaInfo[c("Level.3a", "Level.3s")][taxaInfo[["Level.3a"]] != taxaInfo[["Level.3s"]] & !(taxaInfo[["Level.3s"]] %in% c("Other hardwood", "None")), ]
names(nameConversions) <- c('from', 'to')

numConv <- nrow(nameConversions)
if(numConv) {
  for(i in seq_len(numConv)) {
    if(nameConversions$from[i] %in% names(data)) {
      cat("Grouping ", nameConversions$from[i], " into ", nameConversions$to[i], ".\n")
      if(nameConversions$to[i] %in% names(data)) {
        data[nameConversions$to[i]] <- data[nameConversions$to[i]] + data[nameConversions$from[i]]
      } else data[nameConversions$to[i]] <- data[nameConversions$from[i]]
    }
  }
  data <- data[ , !(names(data) %in% nameConversions$from)]
}


other <- rowSums(data[ , taxaOtherHardwood])

data <- data[ , taxaUse]
ord <- order(names(data))
data <- data[ , ord]

data$"Other hardwood" <- other

nTaxa <- ncol(data)
taxaUse <- names(data)

taxa <- data.frame(taxonID = 1:nTaxa, taxonName = taxaUse, stringsAsFactors = FALSE)

cat("Using the following", nTaxa, "taxa, with", sum(data), "trees, of which", round(100*sum(data$"Other hardwood")/sum(data), 2), "% are in the 'other' category:")
cat(taxaUse, sep = ",")
cat("\n")

########################################################################
# subset to portion of PalEON Albers grid and create graph -------------
########################################################################
source(file.path(codeDir, 'set_domain.R'))


if(buffer > 0) {
  coordFull <- expand.grid(X = xGrid, Y=rev(yGrid))
  coordFull$fullID <- seq_len(nrow(coordFull))
  
  coord$origID <- seq_len(nrow(coord))
  
  dataFull <- as.data.frame(matrix(0, nrow = nrow(coordFull), ncol = nTaxa))
  names(dataFull) <- names(data)
  
  tmp <- merge(coordFull, coord,
               all.x = TRUE, all.y = FALSE)
  
  tmp <- tmp[order(tmp$fullID), ]
  origRows <- !is.na(tmp$origID)
  dataFull[origRows, ] <- data[tmp$origID[origRows], ]

  data <- dataFull[ coordFull$X <= easternLimitOfWesternDomain, ]
  coord <- coordFull[ coordFull$X <= easternLimitOfWesternDomain, ]
  coord$ID <- coord$fullID
} else {
  data <- data[ coord$X <= easternLimitOfWesternDomain, ]
  coord <- coord[ coord$X <= easternLimitOfWesternDomain, ]
}

data <- data[order(coord$ID), ]
coord <- coord[order(coord$ID), ]
dimnames(coord)[[1]] <- coord$ID

m1 <- length(westernDomainX)
m2 <- length(westernDomainY)
nCells <- m1 * m2

type <- nbhdStructure
substring(type, 1 ,1) = toupper(substring(type, 1, 1))
fns <- rep("", 2)
fns[1] <- paste('graph', type, '-',  m1, 'x', m2, '.csv', sep='')
fns[2] <- paste('graphCats', type, '-', m1, 'x', m2, '.csv', sep='')

if(!file.exists(file.path(dataDir, fns[1])) || (nbhdStructure != 'bin' && !file.exists(file.path(dataDir, fns[2])))) {
  fns <- graphCreate(m1, m2, type = nbhdStructure, dir = dataDir)
} 

nbhd <- graphRead(fns[1], fns[2], m1, m2, type = nbhdStructure, dir = dataDir)

if(nbhdStructure == "lindgren_nu1") {
  # remove boundary stuff for now while wait to hear from Finn about boundary correction
  nbhd@entries[nbhd@entries %in% c(-4, -6)] <- -8
  nbhd@entries[nbhd@entries %in% c(4, 10, 11, 18, 19)] <- 20
  nbhdIndices <- list()
  # determine which elements correspond to what types of neighbors for fast filling in MCMC
  nbhdIndices$self <- which(nbhd@entries == 20)
  nbhdIndices$cardNbr <- which(nbhd@entries == -8)
  nbhdIndices$otherNbr <- which(nbhd@entries %in% c(1,2))
} else nbhdIndices <- NULL


########################################################################
# create data objects for MCMC fitting ---------------------------------
########################################################################

data <- as.matrix(data)

taxon <- rep(rep(1:nTaxa, nCells), times = c(t(data)))
total <- rowSums(data)
cell <- rep(which(total > 0), times = total[total > 0])

data <- data.frame(taxon = taxon, cell = cell)

save(nbhdStructure, nbhd, nbhdIndices, m1, m2, nTaxa, nCells, data, coord, taxa, file = paste0(file.path(dataDir, 'westernData.Rda')))
