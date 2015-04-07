#!/usr/bin/Rscript

source("config")

runID <- paste0('western_', runID)

require(ncdf4)

if(!exists('runID'))
  stop("should have 'runID' set")

load(file.path(dataDir, paste0('data_', runID, '.Rda')))

replaceZeroWithThis <- 1e-5 # 1/10 of the 10000 samples

finalNcdfName <- paste0('PLScomposition_', runID, '.nc')
ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))

test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nSamples <- dim(test)[3]

preds <- array(0, c(nCells, nTaxa, nSamples))
#dimnames(preds)[[2]] <- taxa$taxonName
# the m2:1 reverses the S->N orientation in the final netCDF file to match the cell numbering used in the modeling (which was N->S). Note that the netCDF files before burn-in are N->S
for(p in 1:nTaxa) {
  preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1))[ , m2:1, ]
}

attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

pm <- apply(preds, c(1, 2), 'mean')
psd <- apply(preds, c(1, 2), 'sd')


dataCellTest <- as.matrix(dataFull[dataFull$cell %in% holdOutCells, ])[ , c(2,1)]
dataTreeTest <- as.matrix(treeHoldOut)[ , c(2, 1)]
# need columns reversed for picking out relevant probs from preds array

nC <- nrow(dataCellTest)
nT <- nrow(dataTreeTest)

heldOutEmpProbs <- matrix(0, nCells, nTaxa)
for(p in seq_len(nTaxa)) {
  tbl <- table(dataCellTest[dataCellTest[ , 2] == p, 1])
  heldOutEmpProbs[as.numeric(names(tbl)) , p] <- tbl
}

treesPerCell <- rowSums(heldOutEmpProbs)
heldOutEmpProbs <- heldOutEmpProbs / treesPerCell

treeCountCut <- 0
numTrees <- treesPerCell[treesPerCell > treeCountCut]


### checks

if(FALSE) {
  tmp = matrix(pm[,15],m1,m2)[,m2:1]
  image.plot(1:m1,1:m2,tmp)
  
  tmp = matrix(heldOutEmpProbs[,15],m1,m2)[,m2:1]
  image.plot(1:m1,1:m2,tmp)

}


results <- array(NA, c(7, 2, 2))
attributes(results)$dimnames[[1]] <- c('brier','logdensity','mspe','mae','coverage','mean_int_len','median_int_len')
attributes(results)$dimnames[[2]] <- c('cell','tree')
attributes(results)$dimnames[[3]] <- c('post_mean_of_score','score_of_post_mean')

samples <- array(NA, c(7, 2, nSamples))
attributes(samples)$dimnames[[1]] <- c('brier','logdensity','mspe','mae','coverage','mean_int_len','median_int_len')
attributes(samples)$dimnames[[2]] <- c('cell','tree')


######################################################
### Brier score
######################################################

### using posterior mean
phat <- pm[dataCellTest[ , 1], ]
yvals <- matrix(0, nr = nC, nc = nTaxa)
yvals[cbind(seq_len(nC), dataCellTest[, 2])] <- 1
results['brier','cell','score_of_post_mean'] <- sum((phat - yvals)^2)

phat <- pm[dataTreeTest[ , 1], ]
yvals <- matrix(0, nr = nT, nc = nTaxa)
yvals[cbind(seq_len(nT), dataTreeTest[, 2])] <- 1
results['brier','tree','score_of_post_mean'] <- sum((phat - yvals)^2)

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

samples['brier','cell',] <- brierCell
samples['brier','tree',] <- brierTree

results['brier','cell','post_mean_of_score'] <- mean(brierCell)
results['brier','tree','post_mean_of_score'] <- mean(brierTree)


######################################################
### log posterior density
######################################################

### using posterior mean

logDensCellPoint <- sum(log(pm[dataCellTest]))
logDensTreePoint <- sum(log(pm[dataTreeTest]))

if(is.infinite(logDensCellPoint)) {
  tmp <- pm[dataCellTest]
  tmp[tmp == 0] <- replaceZeroWithThis
  logDensCellPoint <- sum(log(tmp))
}
if(is.infinite(logDensTreePoint)) {
  tmp <- pm[dataTreeTest]
  tmp[tmp == 0] <- replaceZeroWithThis
  logDensTreePoint <- sum(log(tmp))
}

results['logdensity','cell','score_of_post_mean'] <- logDensCellPoint
results['logdensity','tree','score_of_post_mean'] <- logDensTreePoint


### using posterior samples
predsTrunc <- preds
predsTrunc[predsTrunc == 0] <- replaceZeroWithThis  # to deal with loss of digits because of sampling of probs given alphas


logDensCell <- rep(0, nSamples)
logDensTree <- rep(0, nSamples)

for(s in seq_len(nSamples)) {
  logDensCell[s] <- sum(log(predsTrunc[ , , s][dataCellTest]))
  logDensTree[s] <- sum(log(predsTrunc[ , , s][dataTreeTest]))
}

samples['logdensity','cell',] <- logDensCell
samples['logdensity','tree',] <- logDensTree

nullLogDensCell <- nC*log(1/nTaxa)
nullLogDensTree <- nT*log(1/nTaxa)

results['logdensity','cell','post_mean_of_score'] <- mean(logDensCell)
results['logdensity','tree','post_mean_of_score'] <- mean(logDensTree)



########################################################
### MSPE/MAE on cell-level
#########################################################


### weighted mae/mspe using posterior mean

tmp1 <- heldOutEmpProbs[treesPerCell > treeCountCut, ]
tmp2 <- apply(preds[treesPerCell > treeCountCut, , ], c(1, 2), mean)

results['mspe','cell','score_of_post_mean'] <- sum(numTrees*(tmp1 - tmp2)^2)/sum(numTrees*nTaxa)
results['mae','cell','score_of_post_mean'] <- sum(numTrees*abs(tmp1 - tmp2))/sum(numTrees*nTaxa)

nullMspe = mean((tmp1 - 1/nTaxa)^2)
nullMae = mean(abs(tmp1 - 1/nTaxa))

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

samples['mspe','cell',] <- mspeWgtSamples
samples['mae','cell',] <- maeWgtSamples


results['mspe','cell','post_mean_of_score'] <- mean(mspeWgtSamples)
results['mae','cell','post_mean_of_score'] <- mean(maeWgtSamples)


########################################################
### coverage
#########################################################

cutoff = 50
wh <- which(treesPerCell >= cutoff)
yhat = heldOutEmpProbs[wh, ]

pr = preds[wh, , ]

qu1 = apply(pr, c(1,2), quantile, .025)
qu2 = apply(pr, c(1,2), quantile, .975)

intLength <- c(mean(qu2-qu1), median(qu2-qu1))
coverage <- mean(yhat <= qu2 & yhat >= qu1)

trCnt = treesPerCell[wh]
prY = pr
for(s in 1:nSamples) 
  for(i in 1:nrow(prY)) 
    prY[i, ,s] = rmultinom(1, trCnt[i], pr[i, ,s])/trCnt[i]

qu1 = apply(prY, c(1,2), quantile, .025)
qu2 = apply(prY, c(1,2), quantile, .975)

results['mean_int_len','cell','post_mean_of_score'] <- mean(qu2-qu1)
results['median_int_len','cell','post_mean_of_score'] <- median(qu2-qu1)
results['coverage','cell','post_mean_of_score'] <- mean(yhat <= qu2 & yhat >= qu1)

save(results, samples, file = file.path(outputDir, paste0('PLScomposition_', runID, '-cvResults.Rda')))
