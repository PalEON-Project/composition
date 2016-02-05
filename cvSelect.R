# code to explore selection of test/train split

dat <- read.csv('/var/tmp/paciorek/paleon/comp/data/western-0.3.csv')

dat$nTrees <- rowSums(dat[ , 3:ncol(dat)])

dat$xid <- 1 + (dat$x - min(dat$x)) / 8000
dat$yid <- 1 + (dat$y - min(dat$y)) / 8000



dat$nCardNeigh <- 0

for(i in which(dat$nTrees > 0)) {
  ind <- which(abs(dat$xid - dat$xid[i]) + abs(dat$yid - dat$yid[i]) == 1)
  dat$nCardNeigh[i] <- sum(dat[ind, 'nTrees'] > 0)
  if(i%%1000 == 0) print(c(date(), i))
}
