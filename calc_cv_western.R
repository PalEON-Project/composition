#!/usr/bin/Rscript

library(ncdf4)
source("config")

if(!exists('uniqueRunID'))
  uniqueRunID <- ""

load(file.path(dataDir, paste0('westernData_', productVersion, '-', uniqueRunID, '.Rda')))


finalNcdfName <- paste0('PLScomposition_western_', productVersion, '-', uniqueRunID, '_release.nc')
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

nC <- nrow(dataCellTest)
nT <- nrow(dataTreeTest)

### checks

if(FALSE) {
  tmp = matrix(pm[,15],m1,m2)[,m2:1]
  image.plot(1:m1,1:m2,tmp)
  
  heldOutEmpProbs <- matrix(0, nCells, nTaxa)
  for(p in seq_len(nTaxa)) {
    tbl <- table(dataCellTest[dataCellTest[ , 2] == p, 1])
    heldOutEmpProbs[as.numeric(names(tbl)) , p] <- tbl
  }
  
  treesPerCell <- rowSums(heldOutEmpProbs)
  heldOutEmpProbs <- heldOutEmpProbs / treesPerCell
  
  tmp = matrix(heldOutEmpProbs[,15],m1,m2)[,m2:1]
  image.plot(1:m1,1:m2,tmp)

  heldOutEmpProbs <- matrix(0, nCells, nTaxa)
  for(p in seq_len(nTaxa)) {
    tbl <- table(dataTreeTest[dataTreeTest[ , 2] == p, 1])
    heldOutEmpProbs[as.numeric(names(tbl)) , p] <- tbl
  }

  treesPerCell <- rowSums(heldOutEmpProbs)
  heldOutEmpProbs <- heldOutEmpProbs / treesPerCell

  tmp = matrix(heldOutEmpProbs[,15],m1,m2)[,m2:1]
  image.plot(1:m1,1:m2,tmp)


}




######################################################
### log posterior density
######################################################

### using posterior mean

logDensCellPoint <- sum(log(pm[dataCellTest]))
logDensTreePoint <- sum(log(pm[dataTreeTest]))

if(is.infinite(logDensCellPoint)) {
  tmp <- pm[dataCellTest]
  tmp[tmp == 0] <- .00001
  logDensCellPoint <- sum(log(tmp))
}
if(is.infinite(logDensTreePoint)) {
  tmp <- pm[dataTreeTest]
  tmp[tmp == 0] <- .00001
  logDensTreePoint <- sum(log(tmp))
}

### using posterior samples
predsTrunc <- preds
predsTrunc[predsTrunc == 0] <- 0.00001  # to deal with loss of digits because of sampling of probs given alphas


logDensCell <- rep(0, nSamples)
logDensTree <- rep(0, nSamples)

for(s in seq_len(nSamples)) {
  logDensCell[s] <- sum(log(predsTrunc[ , , s][dataCellTest]))
  logDensTree[s] <- sum(log(predsTrunc[ , , s][dataTreeTest]))
}

nullLogDensCell <- nC*log(1/nTaxa)
nullLogDensTree <- nT*log(1/nTaxa)

mean(logDensCell)
mean(logDensTree)

######################################################
### Brier score
######################################################

### using posterior mean
phat <- pm[dataCellTest[ , 1], ]
yvals <- matrix(0, nr = nC, nc = nTaxa)
yvals[cbind(seq_len(nC), dataCellTest[, 2])] <- 1
brierCellPoint <- sum((phat - yvals)^2)

phat <- pm[dataTreeTest[ , 1], ]
yvals <- matrix(0, nr = nT, nc = nTaxa)
yvals[cbind(seq_len(nT), dataTreeTest[, 2])] <- 1
brierTreePoint <- sum((phat - yvals)^2)

### using posterior samples

brierCell <- rep(0, nSamples)
brierTree <- rep(0, nSamples)

yvals <- matrix(0, nr = nC, nc = nTaxa)
yvals[cbind(seq_len(nC), dataCellTest[,2])] <- 1
for(s in seq_len(nSamples)) {
  phat <- preds[ , , s][dataCellTest[ , 1], ]
  brierCell[s] <- sum((phat - yvals)^2)
}

yvals <- matrix(0, nr = nT, nc = nTaxa)
yvals[cbind(seq_len(nT), dataTreeTest[,2])] <- 1
for(s in seq_len(nSamples)) {
  phat <- preds[ , , s][dataTreeTest[ , 1], ]
  brierTree[s] <- sum((phat - yvals)^2)
}

mean(brierCell)
mean(brierTree)

########################################################
### MSPE/MAE on cell-level
#########################################################

heldOutEmpProbs <- matrix(0, nCells, nTaxa)
for(p in seq_len(nTaxa)) {
  tbl <- table(dataCellTest[dataCellTest[ , 2] == p, 1])
  heldOutEmpProbs[as.numeric(names(tbl)) , p] <- tbl
}

treesPerCell <- rowSums(heldOutEmpProbs)
heldOutEmpProbs <- heldOutEmpProbs / treesPerCell

treeCountCut <- 0
numTrees <- treesPerCell[treesPerCell > treeCountCut]

### weighted mae/mspe using posterior mean

tmp1 <- heldOutEmpProbs[treesPerCell > treeCountCut, ]
tmp2 <- apply(preds[treesPerCell > treeCountCut, , ], c(1, 2), mean)

mspeWgt = sum(numTrees*(tmp1 - tmp2)^2)/sum(numTrees*nTaxa)
maeWgt = sum(numTrees*abs(tmp1 - tmp2))/sum(numTrees*nTaxa)

nullMspe = mean((tmp1 - 1/nTaxa)^2)
nullMae = mean(abs(tmp1 - 1/nTaxa))

#treeCountCut <- 19
#tmp1 <- heldOutEmpProbs[treesPerCell > treeCountCut, ]
#tmp2 <- apply(preds[treesPerCell > treeCountCut, , ], c(1, 2), mean)

#mspe = mean((tmp1 - tmp2)^2)
#mae = mean(abs(tmp1 - tmp2))

### using posterior samples

mspeSamples <- maeSamples <- mspeWgtSamples <- maeWgtSamples <- rep(0, nSamples)

tmp2 <- preds[treesPerCell > treeCountCut, , ]

for(s in seq_len(nSamples)) {
  mspeSamples[s] <- mean((tmp1 - tmp2[ ,, s])^2)
  maeSamples[s] <- mean(abs(tmp1 - tmp2[ ,, s]))
  mspeWgtSamples[s] <- sum(numTrees*(tmp1 - tmp2[ ,, s])^2)/sum(numTrees*nTaxa)
  maeWgtSamples[s] <- sum(numTrees*abs(tmp1 - tmp2[ ,, s]))/sum(numTrees*nTaxa)
#  logDensTree[s] <- sum(log(preds[ , , s][dataTreeTest]))
}

mean(mspeWgtSamples)
mean(maeWgtSamples)


########################################################
### coverage
#########################################################

cutoff = 50
wh <- which(treesPerCell > cutoff)

yhat = heldOutEmpProbs[wh, ]
pr = preds[wh, , ]

qu1 = apply(pr, c(1,2), quantile, .025)
qu2 = apply(pr, c(1,2), quantile, .975)

c(mean(qu2-qu1), median(qu2-qu1))

mean(yhat <= qu2 & yhat >= qu1)

trCnt = treesPerCell[wh]
prY = pr
for(s in 1:nSamples) 
  for(i in 1:nrow(prY)) 
    prY[i, ,s] = rmultinom(1, trCnt[i], pr[i, ,s])/trCnt[i]

qu1 = apply(prY, c(1,2), quantile, .025)
qu2 = apply(prY, c(1,2), quantile, .975)

c(mean(qu2-qu1), median(qu2-qu1))

mean(yhat <= qu2 & yhat >= qu1)
