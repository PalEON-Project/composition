#!/usr/bin/Rscript
source("config")


source(file.path(codeDir, "graph.R"))


########################################################################
# deal with town-cell intersections-------------------------------------
########################################################################


load(file.path(dataDir, 'intersection.Rda'))
# reads in 'inter'

maxPossibleCells <- 70
townCellOverlap <- matrix(0, nTowns, maxPossibleCells)
townCellIds <- matrix(1, nTowns, maxPossibleCells)

for(k in 1:nTowns) {
  tmp <- inter[ , , k]
  ids <- which(tmp > 0)
  nMatch <- length(ids)
  townCellIds[k, 1:nMatch] <- ids
  townCellOverlap[k, 1:nMatch] <- tmp[ids]
}

tmp <- colSums(townCellOverlap)

townCellOverlap <- townCellOverlap[ , which(tmp > 0)]
townCellIds <- townCellIds[ , which(tmp > 0)]

########################################################################
# read graph for grid --------------------------------------------------
########################################################################


source(file.path(codeDir, 'set_domain.R'))

m1 <- length(easternDomainX)
m2 <- length(easternDomainY)

tmp <- nbhdStructure
substring(tmp, 1 ,1) <- toupper(substring(tmp, 1, 1))
graphFileName <- paste0('graph', tmp, m1, 'x', m2, '.csv')
if(!file.exists(file.path(dataDir, graphFileName)))
  graphCreate(m1, m2, type = nbhdStructure, dir = dataDir, fn = graphFileName)
nbhd <- graphRead(file.path(dataDir, graphFileName), m1, m2)



########################################################################
# read data ------------------------------------------------------------
########################################################################

easternDataDir <- paste0(easternVersionID, '.', easternVersion)
ohioDataDir <- paste0(ohioVersionID, '.', ohioVersion)

fn <- file.path(dataDir, easternDataDir, paste0(easternVersionID, 'polygonsver', easternVersion, '.csv'))
data1 <- read.csv(fn)
names(data1) <- scan(pipe(paste0("head -n 1 ", fn)), sep = ',', what = 'character')
data1 <- data1[order(data1$ID), 8:31]

cat(paste0("Read ", nrow(data1), " rows from ", fn, ", with field names: "))
cat(names(data1), sep = ',')
cat("\n")

fn <- file.path(dataDir, ohioDataDir, paste0("OH", ohioVersionID, "polygonsver", ohioVersion, ".csv"))
data2 <- read.csv(fn)
names(data2) <- scan(pipe(paste0("head -n 1 ", fn)), sep = ',', what = 'character')
data2 <- data2[order(data2$ID), 10:33]

cat(paste0("Read ", nrow(data2), " rows from ", fn, ", with field names: "))
cat(names(data2), sep = ',')
cat("\n")

data1[is.na(data1)] <- 0
data2[is.na(data2)] <- 0



########################################################################
# subset and manipulate taxa ------------------------------------
########################################################################


taxaOther <- scan("easternTaxaOtherHardwood", what = "character", sep = ',')
taxaUse <- scan("easternTaxaUse", what = "character", sep = ',')
taxaNonTree <- ""

otherNames <- unique(c(names(data1), names(data2))[!(c(names(data1), names(data2)) %in% c(taxaUse, taxaOther, taxaNonTree))])
if(length(otherNames)) {
  cat("Found these additional taxa categories, adding them to 'other':")
  cat(otherNames, sep = ',')
}
taxaOther <- c(taxaOther, otherNames)
cat("Grouping the following taxa into an 'Other hardwood' category:")
cat(taxaOther, sep = ',')
cat("\n")

tmp <- c(taxaUse, taxaOther, taxaNonTree)
if(length(tmp) != length(unique(tmp)))
  warning("One or more taxa in the 'use', 'other', 'non-tree' categories overlap.")

nameConversions <- read.csv("easternTaxaNameConversions", header = TRUE, stringsAsFactors = FALSE)

numConv <- nrow(nameConversions)
if(numConv) {
  for(i in seq_len(numConv)) {
    if(!(nameConversions$to[i] %in% taxaUse) || !(nameConversions$from[i] %in% taxaUse)) {
      warning(paste0("one or more names in name conversions file not in file of taxa to use:", nameConversions$to[i], " ", nameConversions$from[i]))
    } else {
      data1[nameConversions$to[i]] <- data1[nameConversions$to[i]] + data1[nameConversions$from[i]]
      data2[nameConversions$to[i]] <- data2[nameConversions$to[i]] + data2[nameConversions$from[i]]
    }
  }
  data1 <- data1[ , !(names(data1) %in% nameConversions$from)]
  data2 <- data2[ , !(names(data2) %in% nameConversions$from)]
  taxaUse <- taxaUse[!(taxaUse%in% nameConversions$from)]
}


allTaxa <- c(taxaOther, taxaUse)
zeros <- allTaxa[!(allTaxa %in% names(data1))]
data1[ , zeros] <- 0
zeros <- allTaxa[!(allTaxa %in% names(data2))]
data2[ , zeros] <- 0

other1 <- rowSums(data1[ , taxaOther])
other2 <- rowSums(data2[ , taxaOther])

data1 <- data1[ , taxaUse]
data2 <- data2[ , taxaUse]
  
data1 <- data1[ , order(names(data1))]
data2 <- data2[ , order(names(data2))]

data <- rbind(data1, data2)

data$"Other hardwood" <- c(other1, other2)

taxaUse <- names(data)
nTaxa <- length(taxaUse)

taxa <- data.frame(taxonID = 1:nTaxa, taxonName = taxaUse, stringsAsFactors = FALSE)


cat("Using the following", nTaxa, "taxa, with", sum(data), "trees, of which", round(100*sum(data$Other)/sum(data), 2), "% are in the 'other' category:")
cat(taxaUse, sep = ",")
cat("\n")

########################################################################
# create data objects for MCMC fitting ---------------------------------
########################################################################


nTrees <- rowSums(data)

town <- rep(seq_len(nTowns), times = nTrees)

tmp <- c(t(as.matrix(data)))

taxon <- rep(rep(1:nTaxa, nTowns), times = tmp)

data <- data.frame(taxon = taxon, town = town)

save(data, townCellOverlap, townCellIds, taxa, nbhd, m1, m2, nTaxa, nTowns, file = paste0(file.path(dataDir, 'easternData.Rda')))
