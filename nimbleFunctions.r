dCJS <- nimbleFunction(
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(1),
                 probCapture = double(1),
                 firstcap = integer(0),
                 nT = double(0, default = 0),
                 nrep = integer(0, default = 1), # Number of repeated identical capture histories with identical survival and capture probabilities. 
                 log = integer(0, default = 0) ## required log argument
  ) {
    if (nT != 0) {
      if (nT != length(x)) stop("Argument len must match length of data, or be 0.")
    }
    if (length(probSurvive) < length(x) - 1)
      stop("Length of probSurvive must be at least length of data minus 1.")
    if (length(x) != length(probCapture)) stop("Length of probCapture does not match length of data.")
    
    ind <- which(x > 0)
    if (firstcap != ind[1])
      stop('First encounter does not match firstcap')
    lastcap <- ind[length(ind)]
    lp <- 0
    if (lastcap > firstcap) {
      for (t in (firstcap + 1):lastcap) {
          lp <- lp + log(probSurvive[t-1]) + x[t] * log(probCapture[t]) + (1 - x[t]) * log(1-probCapture[t])
      }
    }
    
    if (lastcap < nT) {
      pNotSeenAgain <- 1
      for (t in 1:(nT - lastcap)) {
        pNotSeenAgain <- (1 - probSurvive[nT - t]) + probSurvive[nT - t] * (1-probCapture[nT - t + 1]) * pNotSeenAgain
      }                  
      lp <- lp + log(pNotSeenAgain)
    }

    lp <- nrep * lp
    if (log) {
      return(lp)
    }
    return(exp(lp))
    returnType(double())
  }
)
