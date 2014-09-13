#!/usr/bin/Rscript

# assumes inputs in westernData.Rda

library(ncdf4)
source("config")

if(!exists('uniqueRunID'))
  uniqueRunID <- ""

load(file.path(dataDir, 'westernData.Rda'))


finalNcdfName <- paste0('PLScomposition_western_', productVersion, '-', uniqueRunID, '.nc')
ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))

test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nSamples <- dim(test)[3]


preds <- array(0, c(nCells, nTaxa, nSamples))
#dimnames(preds)[[2]] <- taxa$taxonName
for(p in 1:nTaxa) {
  preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1))
}

attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

pm <- apply(preds, c(1, 2), 'mean')
psd <- apply(preds, c(1, 2), 'sd')


dataCellTest <- as.matrix(dataFull[dataFull$cell %in% holdOutCells, ])[ , c(2,1)]
dataTreeTest <- as.matrix(treeHoldOut)[ , c(2, 1)]
# need columns reversed for picking out relevant probs from preds array

logDensCellPoint <- sum(log(pm[dataCellTest]))
logDensTreePoint <- sum(log(pm[dataTreeTest]))


#tmpP <- preds[ , , 30]
preds[preds == 0] <- 0.0001  # to deal with loss of digits because of sampling of probs given alphas

nSamples <- dim(preds)[3]

logDensCell <- rep(0, nSamples)
logDensTree <- rep(0, nSamples)

for(s in seq_len(nSamples)) {
  logDensCell[s] <- sum(log(preds[ , , s][dataCellTest]))
  logDensTree[s] <- sum(log(preds[ , , s][dataTreeTest]))
}

nullLogDens <- nrow(dataCellTest)*log(1/nTaxa)

heldOutEmpProbs <- matrix(0, nCells, nTaxa)
for(p in seq_len(nTaxa)) {
  tbl <- table(dataCellTest[dataCellTest[ , 2] == p, 1])
  heldOutEmpProbs[as.numeric(names(tbl)) , p] <- tbl
}


treesPerCell <- rowSums(heldOutEmpProbs)
heldOutEmpProbs <- heldOutEmpProbs / treesPerCell

treeCountCut <- 19
tmp1 <- heldOutEmpProbs[treesPerCell > treeCountCut, ]
tmp2 <- apply(preds[treesPerCell > treeCountCut, , ], c(1, 2), mean)

mspe = mean((tmp1 - tmp2)^2)
mae = mean(abs(tmp1 - tmp2))

mspeSamples <- maeSamples <- rep(0, nSamples)

tmp2 <- preds[treesPerCell > treeCountCut, , ]

for(s in seq_len(nSamples)) {
  mspeSamples[s] <- mean((tmp1 - tmp2[ ,, s])^2)
  maeSamples[s] <- mean(abs(tmp1 - tmp2[ ,, s]))
#  logDensTree[s] <- sum(log(preds[ , , s][dataTreeTest]))
}

nullMspe = mean((tmp1 - 1/nTaxa)^2)
nullMae = mean(abs(tmp1 - 1/nTaxa))


# weighted mae/mspe

treeCountCut <- 0
tmp1 <- heldOutEmpProbs[treesPerCell > treeCountCut, ]
tmp2 <- apply(preds[treesPerCell > treeCountCut, , ], c(1, 2), mean)
numTrees <- treesPerCell[treesPerCell > treeCountCut]

mspeWgt = sum(numTrees*(tmp1 - tmp2)^2)/sum(numTrees*nTaxa)
maeWgt = sum(numTrees*abs(tmp1 - tmp2))/sum(numTrees*nTaxa)




