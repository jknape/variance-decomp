library(nimble)
source('nimbleFunctions.r')

###############################
# IPM code (BUGS/NIMBLE)
###############################

cmrCode = nimbleCode({
  ###################################################################
  #Likelihood for capture-recapture data: CJS model (3 age classes)
  ###################################################################
  for (yr in 1:(nYears - 1)) {
    for (age in 1:3) {
      # Split linear predictor into different components. For computing contributions.
      lpc.phi[yr, age, 1:np.s] <-
        X.phi[yr, 1:np.s] * beta.phi[1:np.s, age]
      lpc.phi[yr, age, np.s + 1] <-
        mu.phi[age] + s.phi[age] * eps.phi[yr, age]
      lpc.phi[yr, age, np.s + 2] <-
        step(1.5 - age) * dd.phi * (N[yr, 1] + N[yr, 2] - 90) / 10
      phi[yr, age] <- ilogit(sum(lpc.phi[yr, age, 1:(np.s + 2)]))
    }
  }
  for (i in 1:nCap) {
    cmr[i, 1:nYears] ~ dCJS(
      survivalProb[i, 1:(nYears - 1)],
      captureProb[i, 1:nYears],
      firstcap = firstcap[i],
      nT = nYears,
      nrep = nCmrRep[i]
    )
    
    for (yr in firstcap[i]:(nYears - 1)) {
      survivalProb[i, yr] <- phi[yr, age.CMR[i, yr] + 1] *
        (step(age.CMR[i, yr] - .5) + step(.5 - age.CMR[i, yr]) * zeta * psi[yr])
    }
    
    for (yr in (firstcap[i] + 1):nYears) {
      captureProb[i, yr] <- p[age.CMR[i, yr]]
    }
  }
  for (yr in 1:(nYears - 1)) {
    eps.phi[yr, 1] ~ dnorm(0, 1)
    eps.phi[yr, 2] ~ dnorm(0, 1)
    eps.phi[yr, 3] ~ dnorm(0, 1)
  }
  mu.phi[1] ~ dnorm(0, sd = 1.5)
  mu.phi[2] ~ dnorm(0, sd = 1.5)
  mu.phi[3] ~ dnorm(0, sd = 1.5)
  dd.phi ~ dnorm(0, 1)
  
  for (j in 1:3) {
    beta.phi[1, j] ~ dnorm(0, 1) # Sahel rain
  }
  for (i in 2:np.s) {
    beta.phi[i, 1] <- 0 # No coefs for local weather for juveniles
    for (j in 2:3) {
      beta.phi[i, j] ~ dnorm(0, 1)
    }
  }
  s.phi[1] ~ dexp(20)
  s.phi[2] ~ dexp(20)
  s.phi[3] ~ dexp(20)
  
  p[1] ~ dunif(0, 1)
  p[2] ~ dunif(0, 1)
  
  ###################################################################
  #Likelihood for breeding success
  ###################################################################
  
  # Probability of individual juvenile surviving from age of ringing till fledging (psi),
  for (yr in 1:nYears) {
    lpc.psi[yr, 1:np.js] <- X.psi[yr, 1:np.js] * beta.psi[1:np.js]
    lpc.psi[yr, np.js + 1] <- mu.psi + s.psi * eps.psi[yr]
    logit(psi[yr]) <- sum(lpc.psi[yr, 1:(np.js + 1)])
  }
  
  # Probability of nest (=breeding attempt) surviving until age of ringing (nu).
  for (yr in 1:nYears) {
    for (age in 1:2) {
      lpc.nu[yr, age, 1:np.ns] <- X.nu[yr, 1:np.ns] * beta.nu[1:np.ns]
      lpc.nu[yr, age, np.ns + 1] <- mu.nu[age] + s.nu * eps.nu[yr]
      logit(nu[yr, age]) <- sum(lpc.nu[yr, age, 1:(np.ns + 1)])
    }
  }
  
  # Marginal and conditional (on nest survival until ringing) probability of breeding success (BS).
  for (yr in 1:nYears) {
    for (age in 1:2) {
      for (k in 1:8) {
        p.BS.conditional[yr, age, k] <- (1 - (1 - zeta) ^ k) * psi[yr]
      }
      p.BS[yr, age] <-
        nu[yr, age] * sum(p.C[1:8] * p.BS.conditional[yr, age, 1:8]) * (1 -
                                                                                                p.R[yr, age]) / (1 - nu[yr, age] * p.R[yr, age])
    }
  }
  zeta ~ dunif(0, 1)
  
  for (i in 1:np.js) {
    beta.psi[i] ~ dnorm(0, 1)
  }
  for (i in 1:np.ns) {
    beta.nu[i] ~ dnorm(0, 1)
  }
  for (yr in 1:nYears) {
    eps.psi[yr] ~ dnorm(0, 1)
    eps.nu[yr] ~ dnorm(0, 1)
  }
  mu.psi ~ dnorm(0, sd = 1.5)
  
  for (age in 1:2) {
    mu.nu[age] ~ dnorm(0, sd = 1.5)
  }
  
  s.psi ~ dexp(20) #dunif(0,5)
  s.nu ~ dexp(20) #dunif(0,5)
  
  # Categorical distribution for number of ringed chicks in nest
  ringedYoung[1:8] ~ dmulti(size = ringedYoung.phiS, prob = p.C[1:8])	
  p.C[1:8] ~ ddirch(alpha[1:8])
  alpha[1:8] <- rep(1, 8)
  # Compute mean and variance of categorical distribution	
  r <- sum(p.C[1:8] * 1:8)
  r.V <- sum(p.C[1:8] * (1:8 - r)^2)
  
  # Binomial model for breeding success
  for (i in 1:n.BS) {
    BS[i] ~ dbern((1 - BS.ringed[i]) * p.BS[yr.BS[i], age.BS[i]] + BS.ringed[i] * p.BS.conditional[yr.BS[i], age.BS[i], nrBS.ringed[i]])
    BS.ringed[i] ~ dbern(nu[yr.BS[i], age.BS[i]] * p.R[yr.BS[i], age.BS[i]])
  }
  
  for (yr in 1:nYears) {
    p.R[yr, 1] ~ dunif(0, 1)
    p.R[yr, 2] ~ dunif(0, 1)
  }
  
  ##########################################################################
  #Likelihood for counts of nr of breeding territories, pre-breeding census.
  ##########################################################################
  
  N[1, 1] <-
    round(N0[1]) 
  N[1, 2] <-
    round(N0[2]) 
  N0[1] ~ dnorm(50, sd = 40)
  N0[2] ~ dnorm(50, sd = 40)
  for (yr in 2:nYears) {
    for (age in 1:2) {
      # Number of nests surviving until the time of ringing
      N.nu[yr - 1, age] ~ dbin(prob = nu[yr - 1, age], size = N[yr - 1, age])
      # Normal approximation for total number of chicks	
      N.r[yr-1, age] ~ dnorm(N.nu[yr - 1, age] * r, sd = sqrt(N.nu[yr - 1, age] * r.V))

      # Number of recruits (one year old males)
      N.REC[yr, age] ~ dbin(prob = psi[yr - 1] * zeta * phi[yr - 1, 1], size = round(N.r[yr - 1, age]))
      # Number of surviving breeding males
      N.SURV[yr, age] ~ dbin(prob = phi[yr - 1, age + 1], size = N[yr - 1, age])
      eps.omega[yr - 1, age] ~ dnorm(0, 1)
      lpc.omega[yr - 1, age] <-
        log(r.omega[age]) + s.omega[age] * eps.omega[yr - 1, age]
      omega[yr - 1, age] <- exp(lpc.omega[yr - 1, age])
      # Number of immigrants
      N.IMM[yr, age] ~ dpois(omega[yr - 1, age] * (N[yr - 1, 1] + N[yr -
                                                                                  1 , 2]))
    }
    N[yr, 1] <- sum(N.REC[yr, 1:2]) + N.IMM[yr, 1]
    N[yr, 2] <- sum(N.SURV[yr, 1:2]) + N.IMM[yr, 2]
  }
  r.omega[1] ~ dexp(10)
  r.omega[2] ~ dexp(10)
  s.omega[1] ~ dexp(20)
  s.omega[2] ~ dexp(20)
  
  ## Population census model
  for (yr in 1:nYears) {
    for (age in 1:2) {
      count[yr, age] ~ dpois(N[yr, age])
    }
  }
  
  for (yr in 1:nYears)  {
    N.tot[yr] <- sum(N[yr, 1:2])
    N.struct[yr, 1:2] <- N[yr, 1:2] / N.tot[yr]
  }
  
  
  ##########################################################################
  # Compute realized growth rate, with and without demographic stochasticity
  ##########################################################################
  
  for (yr in 1:(nYears - 1))  {
    lam.ES[yr] <-
      (
        sum(nu[yr, 1:2] * r * zeta * psi[yr] * phi[yr, 1] * N[yr, 1:2]) +
          sum(phi[yr, 2:3] * N[yr, 1:2]) +
          sum(omega[yr, 1:2] * sum(N[yr, 1:2]))
      ) / N.tot[yr]
    lam[yr] <- N.tot[yr + 1] / N.tot[yr]
    lam.D[yr] <- lam[yr] - lam.ES[yr]
  }
})


###############################
# Fit IPM with NIMBLE
###############################

load('wheatear.Rdata')

bugsModel <- nimbleModel(cmrCode,
                         data = data,
                         constants = constants,
                         inits = inits)

bugsModel$initializeInfo()

cmodel = compileNimble(bugsModel)

mcmc = buildMCMC(configureMCMC(
  bugsModel,
  monitors = c('p', 'p.BS', 'p.BS.conditional', 
               'dd.phi', 'phi', 'mu.phi', 's.phi', 'eps.phi', 'beta.phi', 'lpc.phi',
               'psi', 'beta.psi', 'eps.psi', 'mu.psi', 's.psi', 'lpc.psi',
               'nu', 'beta.nu', 'eps.nu', 'mu.nu', 's.nu', 'lpc.nu',
               'r.omega', 'eps.omega', 's.omega', 'lpc.omega',
               'zeta', 'r', 'p.R',
               'N', 'N.REC', 'N.SURV', 'N.IMM', 'N.tot', 'N.struct',
               'lam.ES', 'lam.D','lam'
  )
))

cmcmc = compileNimble(mcmc, project = bugsModel)

# This will take a while...
# For final results, longer chains should be used.
sampleList = runMCMC(cmcmc,
                     2100,
                     nburnin = 100,
                     thin = 5,
                     nchains = 4)

rhat = coda::gelman.diag(lapply(sampleList, coda::as.mcmc) , multivariate = FALSE)
range(rhat[[1]][, 1], na.rm = TRUE)


###############################
# Compute contributions
###############################
source('contrib.r')

samples = do.call(rbind, sampleList)

cont  = apply(samples, 1, contrib)

# Posterior distribution of relative contributions from environmental variables
rel.contrib  = t(do.call(cbind, lapply(cont, `[[`, 'rel.contrib')))

# Posterior means and 50% and 95% credible intervals for relative contributions
colMeans(rel.contrib)
apply(rel.contrib, 2, quantile, c(.025, .25,.75, .975))

# Posterior distribution of the proportion of contribution from environmental stochasticity
env.contrib = sapply(cont, `[[`, 'rel.environmental.var')
hist(env.contrib)

# Posterior distribution of the proportion of demographic stochasticity
hist(1 - env.contrib)

# Posterior distribution of relative contributions to demographic stochasticity
rel.contrib.d  = t(do.call(cbind, lapply(cont, `[[`, 'rel.contrib.d')))

# Posterior means and 50% and 95% credible intervals for relative contributions
colMeans(rel.contrib.d)
apply(rel.contrib.d, 2, quantile, c(.025, .25,.75, .975))

# Posterior distribution of proportion of variance captured by Taylor approximation
taylor.diagnostic  = sapply(cont, `[[`, 'taylor.diagnostic')
hist(taylor.diagnostic)

# The relative error is the taylor diagnostic subtracted from 1
rel.error = 1 - taylor.diagnostic

