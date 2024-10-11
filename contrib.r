contrib = function(pars) {
  
  parNames = names(pars)
  
  ######################################################
  # Extract linear predictor terms
  ######################################################
  nyear = sum(grepl('lpc.phi', parNames))/5/3 + 1
  
  lpc.phi = array(pars[startsWith(parNames, 'lpc.phi[')], dim = c(nyear-1, 3, 5))
  lpc.psi = array(pars[startsWith(parNames, 'lpc.psi[')], dim = c(nyear, 3))
  lpc.nu = array(pars[startsWith(parNames, 'lpc.nu[')], dim = c(nyear, 2, 3))
  lpc.omega = array(pars[startsWith(parNames, 'lpc.omega[')], dim = c(nyear - 1, 2))
  N.struct = matrix(pars[startsWith(parNames, 'N.struct[')], ncol = 2)[1:(nyear - 1), ]; colnames(N.struct) = c('N1.struct', 'N2.struct')
  
  # Names of the linear predictor components, important that they match model specification.
  phinames = c('rain.winter' , 'rain.summer', 'temperature.summer', 'random', 'density')
  lpc.phi0 = lpc.phi[,1,]; colnames(lpc.phi0) = paste0('phi0.', phinames)
  lpc.phi1 = lpc.phi[,2,]; colnames(lpc.phi1) = paste0('phi1.', phinames)
  lpc.phi2 = lpc.phi[,3,]; colnames(lpc.phi2) = paste0('phi2.', phinames)
  
  # Sanity check of dimensions
  stopifnot(sum(grepl('lpc.phi', parNames)) ==  prod(dim(lpc.phi0)) +prod(dim(lpc.phi1)) + prod(dim(lpc.phi2)))
  
  psinames = c('rain.nestling', 'temperature.nestling', 'random')
  lpc.psi = lpc.psi[1:(nyear - 1),]; colnames(lpc.psi) = paste0('psi.', psinames)
  # Sanity check of dimensions
  stopifnot(sum(grepl('lpc.psi', parNames)) ==  prod(dim(lpc.psi)) + 3)
  
  nunames = c('rain.incubation', 'temperature.incubation', 'random')
  lpc.nu1 = lpc.nu[1:(nyear - 1),1,]; colnames(lpc.nu1) = paste0('nu1.', nunames)
  lpc.nu2 = lpc.nu[1:(nyear - 1),2,]; colnames(lpc.nu2) = paste0('nu2.', nunames)
  stopifnot(sum(grepl('lpc.nu', parNames)) ==  prod(dim(lpc.nu1)) +prod(dim(lpc.nu2)) + 2 * 3)
  
  omeganames = c('random')
  lpc.omega1 = lpc.omega[,1, drop = FALSE]; colnames(lpc.omega1) = paste0('omega1.', omeganames)
  lpc.omega2 = lpc.omega[,2, drop = FALSE]; colnames(lpc.omega2) = paste0('omega2.', omeganames)
  stopifnot(sum(grepl('lpc.omega', parNames)) ==  length(lpc.omega1) + length(lpc.omega2))
  
  ######################################################
  # Evaluate vital rates at means of linear predictors.
  ######################################################
  
  # Survival
  phi0 = plogis(sum(colMeans(lpc.phi0))) # First year survival
  phi1 = plogis(sum(colMeans(lpc.phi1))) # 1 year olds
  phi2 = plogis(sum(colMeans(lpc.phi2))) # adults
  
  # Survival of individual young in nest
  psi = plogis(sum(colMeans(lpc.psi))) # 1 year old parent
  
  # Survival of nest until time of ringing.
  nu1 = plogis(sum(colMeans(lpc.nu1))) # 1 year old parent
  nu2 = plogis(sum(colMeans(lpc.nu2))) # adult parent
  
  # Immigration rate
  omega1 = exp(mean(lpc.omega1)) # immigration rate into recruits
  omega2 = exp(mean(lpc.omega2)) # immigration rate into adults
  
  # Age structure
  N.struct1 = mean(N.struct[,1])
  N.struct2 = mean(N.struct[,2])
  
  
  ######################################################
  # Extract additional time constant parameters
  ######################################################
  
  # Expected number of chicks in nest
  r = getElement(pars,  'r')
  
  # Survival of nest from ringing of chicks untils fledging.
  zeta = getElement(pars, 'zeta')
  
  ######################################################
  # Compute gradient at means of linear predictors.
  ######################################################
  grad = c(
    phi0 = r * zeta * (psi * nu1 * N.struct1 + psi * nu2 * N.struct2) * phi0 * (1 - phi0), # First year survival
    phi1 = N.struct1 * phi1 * (1 - phi1),
    phi2 = N.struct2 * phi2 * (1 - phi2),
    psi = r * zeta * phi0 * (nu1 * N.struct1 + nu2 * N.struct2) * psi * (1 - psi), # Nestling survival for 1 year old parents
    nu1 = r * zeta * phi0 * psi * N.struct1 * nu1 * (1 - nu1), # Nest survival for 1 year old parents
    nu2 = r * zeta * phi0 * psi * N.struct2 * nu2 * (1 - nu2), # Nest survival for adult parents
    omega1 = 1 * omega1, # Immigration to recruits 
    omega2 = 1 * omega2,  # Immigration to adults
    N1 = r * zeta * phi0 * psi * nu1 + phi1 - # Age structure one year olds
    (r * zeta * phi0 * psi * nu2 + phi2) # Age structure adults
  )
  
  # Put all the annual linear predictor components in a matrix (#years x #components) for computing covariance matrix
  lpc.mat = cbind(lpc.phi0, lpc.phi1, lpc.phi2, 
                  lpc.psi, 
                  lpc.nu1, lpc.nu2, 
                  lpc.omega1, lpc.omega2, 
                  N1 = N.struct[, 1, drop = FALSE])
  
  Sigma = cov(lpc.mat)
  
  # We need the gradient for each linear predictor component, but the gradients are the same for all components of the same linear predictor.
  # We therefore map the gradient of the linear predictors above to the gradient of the linear predictor components using a matrix M:
  fac = sub('\\..*$', '', colnames(lpc.mat))
  fac = factor(fac, levels = unique(fac))
  M = model.matrix(~ f-1, data = data.frame(f = fac))
  rownames(M) = rownames(Sigma)
  
  # Now we can compute the contributions to environmental variance in growth.
  contributions = (M %*% grad) * (Sigma %*% (M %*% grad))
  
  # Relative contributions that sum to 1.
  rel.contrib = contributions / sum(contributions)
  
  # Realized environmental growth rates (without demographic stochasticity), from NIMBLE output
  lam.ES = pars[startsWith(parNames, 'lam.ES[')]
  
  ######################################################
  # Contributions from demographic stochasticity
  ######################################################
  
  # Demographic growth rate
  lam.D = pars[startsWith(parNames, 'lam.D[')]
  
  # Total growth rate including demographic stochasticity
  lam = pars[startsWith(parNames, 'lam[')]
  
  # Compute contributions to overall (both demographic and environmental) variance in growth 
  contributions.total = contributions + (M %*% grad) * apply(lpc.mat, 2, cov, y = lam.D) # A previous version was missing the multiplication with M. Thanks to Beth Ross for spotting the error!
  rel.contrib.total = contributions.total / var(lam)
  
  # Basic decomposition of lam.DD
  N.REC = array(pars[grepl('N.REC\\[', parNames)], dim = c(26,2))[-1,]
  N.SURV = array(pars[grepl('N.SURV\\[', parNames)], dim = c(26,2))[-1,]
  N.IMM = array(pars[grepl('N.IMM\\[', parNames)], dim = c(26,2))[-1,]
  N0 = array(pars[grepl('N\\[', parNames)], dim = c(26,2))[-26,]
  D.REC = N.REC- 
          cbind(plogis(rowSums(lpc.nu1)) * r * plogis(rowSums(lpc.psi)) * zeta * plogis(rowSums(lpc.phi0)) * N0[, 1], 
                plogis(rowSums(lpc.nu2)) * r * plogis(rowSums(lpc.psi)) * zeta * plogis(rowSums(lpc.phi0)) * N0[, 2])
  D.SURV = N.SURV - cbind(plogis(rowSums(lpc.phi1)) * N0[,1], plogis(rowSums(lpc.phi2)) * N0[,2])
  D.IMM = N.IMM - cbind(exp(lpc.omega1) * rowSums(N0), exp(lpc.omega2) * rowSums(N0))
  D = cbind(D.REC / rowSums(N0), D.SURV / rowSums(N0), D.IMM / rowSums(N0))
  colnames(D) = paste0(rep(c('REC', 'SURV', 'IMM'), each = 2), 1:2)
  
  contributions.d = as.matrix(rowSums(cov(D)))
  rel.contrib.d = contributions.d / var(lam.D) 
  
  list(rel.contrib = rel.contrib, 
       taylor.diagnostic = sum(contributions) / var(lam.ES),
       rel.environmental.var = (var(lam.ES)  + cov(lam.ES, lam.D)) / var(lam),
       rel.contrib.total = rel.contrib.total,
       rel.contrib.d = rel.contrib.d
       )
}
