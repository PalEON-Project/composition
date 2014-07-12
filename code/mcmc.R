source('code/mcmcAuxil.R')

library(spam)

runMCMC <-function(y, cell = NULL, C, Cindices = NULL, town = NULL, townCellOverlap = NULL, townCellIds = NULL,
                   S, thin, resumeRun, hyperpar = c(-0.5, 0),
                   joint_sample = TRUE, adaptInterval = 100, adaptStartTaper = 3000,                   
                   areallyAggregated = FALSE, outputNcdfName, taxa, numCores = 1, runID = "",
                   dataDir, outputDir) {

  if(is.null(Cindices))
    # FIXME: deal with bin vs Lindgren better in terms of args to runMCMC() and conditionality in the code here
    stop("This code is set up to use the Lindgren model with nu=1")
  
  exclude <- is.na(y)
  if(areallyAggregated) {
    exclude <- exclude | is.na(town)
  } else exclude <- exclude | is.na(cell)
  
  y    <- y[!exclude]
  N<-length(y)
  Nvec <- 1:N

  if(sum(exclude))
    warning("NAs found in 'y', 'cell', or 'town'; excluding these observations.")

  if(!areallyAggregated) {
    cell <- cell[!exclude]
  } else {
    town <- town[!exclude]
    probs <- townCellOverlap / rowSums(townCellOverlap)
    cumProbs <- apply(probs, 1, cumsum)
    maxCells <- ncol(townCellOverlap)
    maxCellsByTown <- as.integer(apply(townCellOverlap, 1, function(x) sum(x>0)))
    possIds <- t(townCellIds[town,])
    prior <-  t(townCellOverlap[town, ]) 
    # could incorporate identification of cells within sampleCells
    cell <- possIds[cbind(sampleCells_cpp(prior), Nvec)]
  }
    
  P <- length(taxa$taxonName)
  I <- nrow(C)

  etaBounds <- c(0, 5)
  
  sigma_propSD <- rep(0.02, P)
  # Lindgren; eta = log(rho); rho = 1/kappa; kappa^2 = a-4
  eta_propSD <- rep(0.15, P)  # eta likely between 0 and 5 # .02
  # modify this to do adaptive
  numAcceptSigma2 <- rep(0, P)
  numAcceptEta <- rep(0, P)


  infs <- rep(Inf, N)
  negInfs <- rep(-Inf, N)
    
  S_keep <- floor(S/thin) 
  S      <- S - S%%S_keep
  
  if(resumeRun) {
    load(file.path(dataDir, paste0('lastState', runID, '.Rda')))
    .Random.seed <<- .Random.seed
    sampleIterates <- (s+1):S
    
  } else {
    sampleIterates <- 1:S
    storeIndex <- 1

    W <- matrix(0, N, P)
    
    sigma2store <- etaStore <- matrix(0, S_keep, P)
    #alpha_current <- alpha_next <- matrix(0, I, P)
    alpha_current <- alpha_next <- matrix(rnorm(I*P), I, P)

#    sigma2_current <- sigma2_next <- rep(1,P)
    sigma2_current <- sigma2_next <- runif(P, 0.1, 3)
    eta_current <- eta_next <- runif(P, 0, 5) # rep(log(3), P)
  }

  n <- rep(0, I)
  tbl <- table(cell)
  n[as.numeric(names(tbl))] <- tbl
  
  Wi.pbar<-matrix(0, I, P)
    
  Vinv <- C  # start with TPS
  a <- 4 + 1/exp(2*eta_current[1])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
  Vinv@entries[Cindices$self] <- 4 + a*a
  Vinv@entries[Cindices$cardNbr] <- -2 * a

  # do initial Cholesky based on sparseness pattern, with first taxon
  B <- Vinv
  diag(B) <- diag(B) + n
  if(!is.spam(B)) warning("B is not spam")
  U <- chol(B)  
    
  # Gathering some indices outside the loop

  treeInd <- treeNonInd <- cellInd <- cellNonInd <- list() 
  
  for(p in 1:P){
    treeInd[[p]]  <- which(y==p)
    treeNonInd[[p]] <- which(y!=p)
    cellInd[[p]]  <- cell[y==p]
    cellNonInd[[p]] <- cell[y!=p]
  }
  
  nTreeInd  <- lapply(treeInd, length)
  nTreeNonInd <- lapply(treeNonInd, length)
  
  for(s in sampleIterates){
    
    # sample the latent W's
    for(p in 1:P){
      if(nTreeInd[[p]] == 1){
          # need to use max, since don't have a matrix in this case; just use simple R version rtruncnorm
        W[treeInd[[p]], p] <- rtruncnorm(nTreeInd[[p]], alpha_current[cellInd[[p]], p],
                                         max(W[treeInd[[p]], -p]), Inf) 
      } else {
          # since this calc involves smaller objects, do simple rowmax_cpp as rowmax2_cpp doesn't help here
        W[treeInd[[p]], p] <- rtruncnorm_cpp(nTreeInd[[p]], alpha_current[cellInd[[p]], p],
                                             rowmax_cpp(W[treeInd[[p]],-p]), Inf)
      }
      if(numCores > 1) {
        # this saves another 0.5 s or so, but means that RNG is not same as with single core operation
 
        W[treeNonInd[[p]], p] <- rtruncnorm_cpp_mp(nTreeNonInd[[p]], alpha_current[cellNonInd[[p]], p],
                                                   -Inf, rowmax2_cpp_mp(W, treeNonInd[[p]], p))
        # rowmax2_cpp_mp does not seem to help in terms of speed - 3.4sec/it vs 3.7sec per it
      } else {
        # non-MP version
        W[treeNonInd[[p]], p] <- rtruncnorm_cpp(nTreeNonInd[[p]], alpha_current[cellNonInd[[p]], p],
                                                -Inf, rowmax2_cpp(W, treeNonInd[[p]], p))
      }
      
    }

    # summarize the W's
    Wi.pbar<-compute_cell_sums_cpp(W,cell,I,P)/c(n+1*(n==0))

    # update alpha's and sigma2's
   if(joint_sample) {
    if(!identical(hyperpar, c(-0.5, 0)))
      stop("Error (runMCMC): joint sampling not set up to use any prior other than flat on sd scale.")

    for(p in 1:P){
        
        sigma2_next[p] <- (rnorm(1, sqrt(sigma2_current[p]), sigma_propSD[p]))^2
        if(sigma2_next[p] < 0) {
          accept <- FALSE 
        } else {

          # terms for reverse proposal
          Vinv <- C  # start with TPS
          a <- 4 + 1/exp(2*eta_current[p])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
          Vinv@entries[Cindices$self] <- 4 + a*a
          Vinv@entries[Cindices$cardNbr] <- -2 * a
          Vinv <- C / sigma2_current[p]

          B <- Vinv
          diag(B) <- diag(B) + n
          
          U <- update.spam.chol.NgPeyton(U, B)
          
          denominator <- -(I)*log(sigma2_current[p])/2 - sum(log(diag(U)))
          UtWi <- forwardsolve(U, Wi.pbar[ , p] * n)
          denominator <- denominator + 0.5*sum(UtWi^2)
          
          # terms for forward proposal
          Vinv <- C  # start with TPS
          a <- 4 + 1/exp(2*eta_current[p])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
          Vinv@entries[Cindices$self] <- 4 + a*a
          Vinv@entries[Cindices$cardNbr] <- -2 * a
          Vinv <- C / sigma2_next[p]
          # simplify as Vinv from above *sigma2_current[p]/sigma2_next[p]
          
          B <- Vinv
          diag(B) <- diag(B) + n
          
          U <- update.spam.chol.NgPeyton(U, B)
          
          numerator <- -(I)*log(sigma2_next[p])/2 - sum(log(diag(U)))
          UtWi <- forwardsolve(U, Wi.pbar[ , p] * n)
          numerator <- numerator + 0.5*sum(UtWi^2)
          
          accept <- decide(numerator - denominator)
        }
        if(accept) {
          numAcceptSigma2[p] <- numAcceptSigma2[p] + 1
          # sample alphas
          means <- backsolve(U, UtWi)
          alpha_next[,p] <- means + backsolve(U, rnorm(I))
        } else {
          sigma2_next[p] <- sigma2_current[p]
          alpha_next[,p] <- alpha_current[,p]
        }
      }

    }
    sigma2_current <- sigma2_next
    alpha_current  <- alpha_next
    
                                        # for the moment, also include the non-joint sample as well so that alphas can move on their own; this adds about a second per iteration
    for(p in 1:P){
      
      Vinv <- C  # start with TPS
      a <- 4 + 1/exp(2*eta_current[p])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
      Vinv@entries[Cindices$self] <- 4 + a*a
      Vinv@entries[Cindices$cardNbr] <- -2 * a
      Vinv <- C / sigma2_current[p]
      
      B <- Vinv
      diag(B) <- diag(B) + n
      
      U <- update.spam.chol.NgPeyton(U, B)
      means <- backsolve(U, forwardsolve(U, Wi.pbar[ , p] * n))
      
      alpha_next[,p] <- means + backsolve(U, rnorm(I))
      
      ss <- sigma2_current[p] * t(alpha_next[,p]) %*% (Vinv %*% alpha_next[,p])
      sigma2_next[p] <- 1 / rgamma(1, shape = hyperpar[1] + (I)/2, scale = 1/(.5*ss + hyperpar[2])) 
      # need I-1 because of the zero eigenvalue in the precision matrix (but not for Lindgren)

      if(is.na(sigma2_next[p])) { # sigma2 too small
        sigma2_next[p] <- sigma2_current[p]
        alpha_next[,p] <- alpha_current[,p]
      }
    }
    sigma2_current <- sigma2_next
    alpha_current  <- alpha_next

    # eta sample
    for(p in 1:P){
      
      eta_next[p] <- rnorm(1, eta_current[p], eta_propSD[p])
      if(eta_next[p] < etaBounds[1] || eta_next[p] > etaBounds[2]) {
        accept <- FALSE 
      } else {
        
                                        # terms for reverse proposal
        Vinv <- C  # start with TPS
        a <- 4 + 1/exp(2*eta_current[p])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
        Vinv@entries[Cindices$self] <- 4 + a*a
        Vinv@entries[Cindices$cardNbr] <- -2 * a
        Vinv <- C / sigma2_current[p]
        
        B <- Vinv
        diag(B) <- diag(B) + n
        
        U <- update.spam.chol.NgPeyton(U, B)
        
        denominator <- -(I)*log(sigma2_current[p])/2 - sum(log(diag(U)))
        UtWi <- forwardsolve(U, Wi.pbar[ , p] * n)
        denominator <- denominator + 0.5*sum(UtWi^2)
        denominator <- denominator + eta_current[p]
        
                                        # terms for forward proposal
        Vinv <- C  # start with TPS
        a <- 4 + 1/exp(2*eta_next[p])  # a = 4 + kappa^2 = 4 + (1/rho)^2 = 4 + (1/exp(eta)^2)
        Vinv@entries[Cindices$self] <- 4 + a*a
        Vinv@entries[Cindices$cardNbr] <- -2 * a
        Vinv <- C / sigma2_current[p]
                                        # simplify as Vinv from above *eta_current[p]/eta_next[p]
        
        B <- Vinv
        diag(B) <- diag(B) + n
        
        U <- update.spam.chol.NgPeyton(U, B)
        
        numerator <- -(I)*log(sigma2_current[p])/2 - sum(log(diag(U)))
        UtWi <- forwardsolve(U, Wi.pbar[ , p] * n)
        numerator <- numerator + 0.5*sum(UtWi^2)
        numerator <- numerator + eta_next[p]
       
        accept <- decide(numerator - denominator)
      }
      if(accept) {
        numAcceptEta[p] <- numAcceptEta[p] + 1
                                        # sample alphas
        means <- backsolve(U, UtWi)
        alpha_next[,p] <- means + backsolve(U, rnorm(I))
      } else {
        eta_next[p] <- eta_current[p]
        alpha_next[,p] <- alpha_current[,p]
      }
    }

    eta_current <- eta_next
    alpha_current  <- alpha_next
#    print(c(s, exp(eta_current)))

    # sample cell indices
    if(areallyAggregated) {
      probs <- calcProbs(t(W), t(alpha_next), prior, possIds, maxCellsByTown, town)
      whichCell <- sampleCells_cpp(probs)
 
      cell <- possIds[cbind(whichCell, Nvec)];
      n <- table_cpp(as.integer(cell), I)

      for(p in 1:P){
        cellInd[[p]]  <- cell[treeInd[[p]]]
        cellNonInd[[p]] <- cell[treeNonInd[[p]]]
      }
    }

    if (s%%thin==0){
      outputNcdfPtr <- nc_open(file.path(dataDir, outputNcdfName), write = TRUE)
      sigma2store[storeIndex, ] <- sigma2_next
      etaStore[storeIndex, ] <- eta_next
      for(p in 1:P) 
        ncvar_put(outputNcdfPtr, taxa$taxonName[p], matrix(alpha_next[ , p], m1, m2) , start = c(1, 1, storeIndex), count = c(-1, -1, 1))
      nc_close(outputNcdfPtr)
      storeIndex <- storeIndex + 1
    }
    
 
    if(s%%100 == 0)
      cat("Finished MCMC iteration ", s, " at ", date(), ".\n", sep = "")

    if(s%%adaptInterval == 0) {
      cat("Acceptance rate for sigma/alpha joint proposals: ", round(numAcceptSigma2/adaptInterval, 2), ".\n", sep = " ")
      numAcceptSigma2 <- rep(0, nTaxa)
      cat("Acceptance rate for eta/alpha joint proposals: ", round(numAcceptEta/adaptInterval, 2), ".\n", sep = " ")
      numAcceptEta <- rep(0, nTaxa)
      print(c(s, exp(eta_current)))
    }
    
    if(s %% 250 == 0) {
      save(alpha_next, sigma2_next, eta_next, s, cell, .Random.seed, W, sigma2store, etaStore, storeIndex, file = file.path(dataDir, paste0('lastState', runID, '.Rda')))
    }
  } # end for s loop

  save(sigma2store, file = file.path(outputDir, paste0("sigma2", runID, ".Rda")))
  save(etaStore, file = file.path(outputDir, paste0("eta", runID, ".Rda")))
  invisible(NULL)
}

drawProportions <- function(latentNcdfPtr, outputNcdfPtr, numMCsamples = 1000, numInputSamples, secondThin = 1, I, taxa){

  samples <- seq(1, numInputSamples, by = secondThin)
  P <- length(taxa)

  tmp <- matrix(0, I, P)
  
  for(s in seq_along(samples)) {

    for(p in 1:P) 
      tmp[ , p] <- ncvar_get(latentNcdfPtr, varid = taxa[p], start = c(1, 1, samples[s]), count = c(-1, -1, 1))
    
    phat <- compute_cell_probabilities_cpp(tmp, numMCsamples, I, P)

    for(p in 1:P) 
      ncvar_put(outputNcdfPtr, taxa[p], phat[ , p], start = c(1, 1, s), count = c(-1, -1, 1))

    if(s%%10 == 0)
      cat("Finished sampling probabilities for (thinned) MCMC sample number ", samples[s], " at ", date(), ".\n", sep = "")
  }
  
  invisible(NULL)
}



