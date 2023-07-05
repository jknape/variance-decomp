This repository contains code for running analyses described in (Knape et al. 2023).

Knape, J., Paquet, M., Arlt, D., Kačergytė, I., and Pärt, T. 2023 Partitioning variance in population growth for models with environmental and demographic stochasticity. Journal of Animal Ecology (preprint: https://doi.org/10.32942/X2FK5D).

# Tutorial

In this tutorial we aim to illustrate how to compute contributions of
environmental covariates to environmental variance in realized growth
rates, and how to compute contributions from demographic stochasticity,
for an integrated population model. We will not cover implementation of
integrated population models in detail, as there are several other
sources available for this. We will also show only the key components of
calculating contributions using code snippets. However, complete model
code, including data, for fitting the wheatear IPM, and code for
computing all contributions is available in Supplement 1.

All code below will be based on MCMC chains from NIMBLE, but the
tutorial should be useful for IPMs coded in other Bayesian modelling
software such as JAGS. In this case some minor modifications of the code
may be necessary because there are small difference in the BUGS code,
and the naming of parameters in the MCMC output differs slightly among
some of the software.

## Computing environmental contributions to variation in realized environmental growth rates

### Coding of linear predictors

We will compute contributions in R using MCMC output. To simplify this
we first set up the BUGS code in a specific way so that we have easy
access to things we need for calculating contributions. Apart from
simplifying the calculation of contributions, this has the additional
advantage of reducing the risk of mismatch between the BUGS code and the
R code used to compute contributions. E.g., if the model is changed and
the NIMBLE code is updated, we want to minimize the amount of changes we
need to make in the computation of contributions. The approach we take
to do this is to make all the components of the linear predictors of
each vital rate readily available in the MCMC output. As an example,
variation in nest survival is modeled via a logistic link function as
logit(*ν*<sub>*a*, *t*</sub>) = *α*<sub>*a*</sub> + *β*<sub>1</sub>rain<sub>*t*</sub> + *β*<sub>2</sub>temp<sub>*t*</sub> + *ϵ*<sub>*t*</sub>.
where *a* is 1 if the male is young and 2 if the male is old, and
*ϵ*<sub>*t*</sub> ∼ *N*(0,*σ*<sub>*ν*</sub><sup>2</sup>). For each of
the two ages (*a* = 1 and *a* = 2), we define three components of this
linear predictor: *β*<sub>1</sub>rain<sub>*t*</sub>,
*β*<sub>2</sub>temp<sub>*t*</sub>, and
*α*<sub>*a*</sub> + *ϵ*<sub>*t*</sub>. The intercept is included in the
random term since it does not vary over time and therefore doesn’t have
a contribution of its own. It could alternatively be included in either
of the covariate terms, or in a separate component with zero
contribution.

In the BUGS code, we collect all the linear predictor components for
each age and time step into a three dimensional array indexed by time,
age, and the predictor component (the dimension in this case is 26 x 2 X
3). This is done by the following lines in the BUGS model code (see file
`script.r` in Supplement 1):

      # Probability of nest (=breeding attempt) surviving until age of ringing (nu).
      for (yr in 1:nYears) {
        for (age in 1:2) {
          lpc.nu[yr, age, 1:np.ns] <- X.nu[yr, 1:np.ns] * beta.nu[1:np.ns]
          lpc.nu[yr, age, np.ns + 1] <- mu.nu[age] + s.nu * eps.nu[yr]
          logit(nu[yr, age]) <- sum(lpc.nu[yr, age, 1:(np.ns + 1)])
        }
      }

Here, the `X.nu` design matrix has the rain covariates as its first
column and the temperature covariates as its second. We are also using
matrix multiplication to simultaneously compute the linear predictor
components for both covariates.

From the `lpc.nu` array, the linear predictors for both ages and all
time points can be simply computed as the sum across the last dimension
of the array. This is then used to compute the nest survival
probabilities using the logit link function.

We use similar coding for all the time varying vital rates in the IPM,
and then make sure to monitor all the linear predictor component arrays
in the MCMC runs so that we can easily access them in the MCMC output.

This set up should make the code reasonably robust to changes as we can
change the linear predictors in the BUGS code, e.g. by adding or
removing variables, without messing up the computation of the
contributions below.

### Computing contributions from MCMC output in R

Having fitted the model using MCMC (in NIMBLE in our case, see the file
‘script.r’ for details) we collect the MCMC output in a matrix called
`samples`. In this output we have discarded initial burn-in samples, and
concatenated any multiple chains to get a single matrix (e.g. using
`rbind`). When fitting the model with NIMBLE, the matrix can be obtained
using this code:

    samples = do.call(rbind, sampleList)

where `sampleList` is the MCMC output from the NIMBLE `runMCMC`
function.

Each row of the `samples` matrix contains a draw from the posterior
distribution of the IPM. We will define an R function, `contrib` (see
file `contrib.r` in Supplement 1), that computes the contributions for a
single set of parameters, in other words for a single row of the
`samples` matrix. The first row of the function therefore is:

    contrib = function(pars) { ...

where the `pars` argument should be a single row of the `samples`
matrix. Later, we can apply the function across all the rows of the
`samples` matrix to get the posterior distribution of the contributions.

The `pars` argument is a named vector where the names represent which
parameter the corresponding value refers to. Calculations based directly
on this vector can become tricky due to the need to index it correctly.
To aid with this we first extract the linear predictor components in the
vector into R arrays. The first few lines of the `contrib` function does
exactly this to retrieve the linear predictor components of each of the
time varying vital rates.

      lpc.phi = array(pars[startsWith(parNames, 'lpc.phi[')], dim = c(nyear-1, 3, 5))
      lpc.psi = array(pars[startsWith(parNames, 'lpc.psi[')], dim = c(nyear, 3))
      lpc.nu = array(pars[startsWith(parNames, 'lpc.nu[')], dim = c(nyear, 2, 3))
      lpc.omega = array(pars[startsWith(parNames, 'lpc.omega[')], dim = c(nyear - 1, 2))
      N.struct = matrix(pars[startsWith(parNames, 'N.struct[')], ncol = 2)[1:(nyear - 1), ]
      colnames(N.struct) = c('N1.struct', 'N2.struct')

In the above linear predictor components, age is one of the dimensions
of the arrays. We next split them up into matrices so that we for each
age class and linear predictor get a matrix containing all the linear
predictor components for each time step. This is mainly for convenience
to avoid having to retrieve them by index later. At the same time we add
the covariate names to the columns of the matrices. We also drop the
last year for some of predictor components that will not be used to for
calculating growth rates, so that all matrices have the same number of
rows (years, 25 time steps in our case).

      # Names of the linear predictor components, important that they match model specification.
      phinames = c('rain.winter' , 'rain.summer', 'temperature.summer', 'random', 'density')
      lpc.phi0 = lpc.phi[,1,]; colnames(lpc.phi0) = paste0('phi0.', phinames)
      lpc.phi1 = lpc.phi[,2,]; colnames(lpc.phi1) = paste0('phi1.', phinames)
      lpc.phi2 = lpc.phi[,3,]; colnames(lpc.phi2) = paste0('phi2.', phinames)
      
      # Sanity check of dimensions
      stopifnot(sum(grepl('lpc.phi', parNames)) ==  
                  prod(dim(lpc.phi0)) +prod(dim(lpc.phi1)) + prod(dim(lpc.phi2)))
      
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

The formula for computing the contributions in the main text is:
contribution(*θ*) := ∇*f*<sub>*θ̄*</sub> ∘ *Σ*<sub>*θ*</sub> ∇*f*<sub>*θ̄*</sub>
where (∇*f*)<sub>*θ̄*</sub> is the gradient of the realized growth rate
as a function of all the linear predictor components, evaluated at the
means of those components, and *Σ*<sub>*θ*</sub> is the covariance
matrix (over time) for those components. We therefore need to compute
these two quantities.

To compute the gradient, we first evaluate all of the vital rates at the
means (over time) of the linear predictors. In other words, we compute
point estimates for each of the vital rates by using the means of the
covariate effects. This is done in the code below.

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

Next, we compute the gradient using the above point estimates. To do
this we need to take the derivative of the realized growth rate with
respect to the linear predictors for each vital rate, including the
inverse of the link functions. This is a straightforward application of
derivatives, but it’s easy to do mistakes so one should do it carefully.
We have done it by hand, but it’s also possible to use symbolic
differentiation software, e.g. using the function `D` in R.

To illustrate, we show how to compute the derivative of the realized
growth rate with respect to the linear predictor of first year survival.
The realized growth rate of the model in the main text, omitting
subscripts for time is:
*λ* = (*ν*<sub>1</sub>*r**ψ**ζ**ϕ*<sub>0</sub>+*ω*<sub>1</sub>)*ñ*<sub>1</sub> + (*ν*<sub>2</sub>*r**ψ**ζ**ϕ*<sub>0</sub>+*ω*<sub>1</sub>)(1−*ñ*<sub>1</sub>) − (*ϕ*<sub>1</sub>+*ω*<sub>2</sub>)*ñ*<sub>1</sub> − (*ϕ*<sub>2</sub>+*ω*<sub>2</sub>)(1−*ñ*<sub>1</sub>)
To take the derivative of this with respect to the linear predictor of
*ϕ*<sub>0</sub> we can use the chain rule to first take the derivative
of *λ* with respect to *ϕ*<sub>0</sub>, and then multiply by the
derivative of *ϕ*<sub>0</sub> with respect to its linear predictor. The
first derivative is:
$$
\frac{\partial \lambda}{\partial \phi\_0} = \nu\_{1} r \psi \zeta  \tilde{n}\_{1} + \nu\_{2} r \psi \zeta (1- \tilde{n}\_{1}),
$$
and the derivative of $\phi\_0 = \frac{\exp(x)}{1+\exp(x)}$ with respect
to *x* (representing the linear predictor) is
*ϕ*<sub>0</sub>(1−*ϕ*<sub>0</sub>). By multiplying these two we get the
derivative of first year survival we required.

To do the same with symbolic differentiation, we could write an R
expression for *λ* as a function of *x* and compute the derivative:

    lambda.expr = expression((nu1 * r * psi * zeta * phi0 + omega1) * n1 + 
                              (nu2 * r * psi * zeta * phi0 + omega1) * (1- n1) - 
                                (phi1  + omega2) * n1  - 
                                  (phi2 + omega2) * (1 - n1))
    phi0.expr = expression(1/(1+exp(-x)))

    # Using chain rule
    D(lambda.expr, 'phi0')

    ## nu1 * r * psi * zeta * n1 + nu2 * r * psi * zeta * (1 - n1)

    D(phi0.expr, 'x')

    ## exp(-x)/(1 + exp(-x))^2

This is identical to the derivatives above. The expression for the
derivative of *ϕ*<sub>0</sub> looks a little different at first sight.
It is in fact equivalent to *ϕ*<sub>0</sub>(1−*ϕ*<sub>0</sub>).

Alternatively we can substitute
$\phi\_0 = \frac{\exp(x)}{1+\exp(x)} = \frac{1}{1 + \exp(-x)}$ directly
in the expression for *λ*:

    D(expression( (nu1 * r * psi * zeta * 1 / (1 + exp(-x)) + omega1) * n1 + 
                    (nu2 * r * psi * 1 / (1 + exp(-x)) + omega1) * (1- n1) - 
                    (phi1  + omega2) * n1  - (phi2 + omega2) * (1 - n1)), 'x')

    ## nu1 * r * psi * zeta * 1 * exp(-x)/(1 + exp(-x))^2 * n1 + nu2 * 
    ##     r * psi * 1 * exp(-x)/(1 + exp(-x))^2 * (1 - n1)

Computing the derivatives, by hand or with symbolic differentiation, for
each of the vital rates plus age structure (represented by the
proportion of young males), and collecting them in a vector gives:

    grad = c(
        phi0 = r * zeta * (psi * nu1 * N.struct1 + psi * nu2 * N.struct2) * phi0 * (1 - phi0), 
        phi1 = N.struct1 * phi1 * (1 - phi1),
        phi2 = N.struct2 * phi2 * (1 - phi2),
        psi = r * zeta * phi0 * (nu1 * N.struct1 + nu2 * N.struct2) * psi * (1 - psi),
        nu1 = r * zeta * phi0 * psi * N.struct1 * nu1 * (1 - nu1),
        nu2 = r * zeta * phi0 * psi * N.struct2 * nu2 * (1 - nu2),
        omega1 = 1 * omega1,
        omega2 = 1 * omega2,
        N1 = r * zeta * phi0 * psi * nu1 + phi1 - 
        (r * zeta * phi0 * psi * nu2 + phi2)
      )

The next thing we need is the covariance matrix of all the components of
all linear predictors, and age structure. To compute them, we combine
all the lpc matrices in a single matrix (with years as rows) and compute
the empirical covariance matrix:

      lpc.mat = cbind(lpc.phi0, lpc.phi1, lpc.phi2, 
                      lpc.psi, 
                      lpc.nu1, lpc.nu2, 
                      lpc.omega1, lpc.omega2, 
                      N1 = N.struct[, 1, drop = FALSE])
      
      Sigma = cov(lpc.mat)

Now the last thing we have to deal with is that the gradient above is
computed over the linear predictors, plus age structure, and is a vector
of length 9, while we computed the empirical covariance for each linear
predictor component, plus age structure, so that `Sigma` is a matrix of
size 27 x 27. We therefore need to compute the gradient over all the
linear predictor components to get a vector of length 27. However, the
derivatives for the linear predictor components belonging to the same
linear predictor are identical. We therefore just need to map the
gradient we computed above to the gradient for the linear predictor
components. This is done using a matrix `M`:

    fac = sub('\\..*$', '', colnames(lpc.mat))
    fac = factor(fac, levels = unique(fac))
    M = model.matrix(~ f-1, data = data.frame(f = fac))
    rownames(M) = rownames(Sigma)

We can now compute the contributions from each of the linear predictor
components using

    contributions = (M %*% grad) * (Sigma %*% (M %*% grad))

which corresponds to the formula for *f* above and in the main text.

The relative contributions are computed by dividing by the total
contribution:

    rel.contributions = contributions / sum(contributions)

### Assessing the Taylor approximation

To assess how well the Taylor approximation that we used above to
compute the contributions is performing, we compare the variance of this
environmental growth rate to the sum of the variance contributions. To
do this, we first need to compute the exact (i.e. approximation free)
annual environmental growth rate. This is done in the BUGS code (file
`script.r`).

    lam.ES[yr] <-
          (
            sum(nu[yr, 1:2] * r * zeta * psi[yr] * phi[yr, 1] * N[yr, 1:2]) +
              sum(phi[yr, 2:3] * N[yr, 1:2]) +
              sum(omega[yr, 1:2] * sum(N[yr, 1:2]))
          ) / N.tot[yr]

We use this exact growth rate in the file `contrib.r` to compare the
variances.

    lam.ES = pars[startsWith(parNames, 'lam.ES[')]
    taylor.diagnostic = sum(contributions) / var(lam.ES)
    rel.error = 1 - taylor.diagnostic

## Estimating contributions from demographic stochasticity

### Incorporating demographic stochasticity in the IPM

Demographic stochasticity is included in the census component of the IPM
(in addition to in the demographic data components). There, starting
with `N[yr - 1, 1]` young males and `N[yr - 1, 2]` older males in year
`yr - 1`, the number of males in year `yr` is simulated by, in sequence,
drawing the number of individuals resulting from each vital rate.

First, the number of nests surviving until the time of ringing for each
age class is drawn from a binomial distribution with rate equal to the
annual probability of nest survival:

    # Number of nests surviving until the time of ringing
    N.nu[yr - 1, age] ~ dbin(prob = nu[yr - 1, age], size = N[yr - 1, age])

Next, the total number of chicks at time of ringing is drawn as the sum
of draws from the categorical distribution for the number of chicks in a
single nest across all nests that survive until ringing (given by `N.nu`
above). For this, we use a normal approximation of the sum, as a direct
coding using a categorical distribution led to convergence problems with
the MCMC. The approximation should be accurate due to the central limit
theorem as we are summing over a large number of nests in each year (it
could also be used for other demographic rates).

    # Normal approximation for total number of chicks   
    N.r[yr-1, age] ~ dnorm(N.nu[yr - 1, age] * r, sd = sqrt(N.nu[yr - 1, age] * r.V))

The number of recruits into the young age class from inside the
population (i.e. excluding immigrants) can then be drawn from a binomial
distribution with size equal to the number of chicks and using the
annual probabilities of nestling survival and first year survival after
fledging.

    # Number of recruits (one year old males)
    N.REC[yr, age] ~ dbin(prob = psi[yr - 1] * zeta * phi[yr - 1, 1], size = round(N.r[yr - 1, age]))

In a similar way, the number of breeding males surviving from year
`yr - 1` to year `yr` is drawn for each age class from binomial
distributions with rates equal to the annual survival probabilities,

    # Number of surviving breeding males
    N.SURV[yr, age] ~ dbin(prob = phi[yr - 1, age + 1], size = N[yr - 1, age])

and the number of immigrants is drawn from poisson distributions

    # Number of immigrants
    N.IMM[yr, age] ~ dpois(omega[yr - 1, age] * (N[yr - 1, 1] + N[yr -1 , 2]))

The number of young males in year `yr` is then the sum of the number of
internals recruits and the number of immigrant young, and the number of
old males is the sum of the number of surviving breeding males and the
number of immigrant old males.

    N[yr, 1] <- sum(N.REC[yr, 1:2]) + N.IMM[yr, 1]
    N[yr, 2] <- sum(N.SURV[yr, 1:2]) + N.IMM[yr, 2]

Finally, the total population size is computed as the sum of young and
old individuals

    N.tot[yr] <- sum(N[yr, 1:2])

### Computing *λ*<sup>*D*</sup> and relative contributions of demographic and environmental stochasticity

To compute *λ*<sup>*D*</sup>, we first compute the total growth rate
including both environmental and demographic stochasticity.

    lam[yr] <- N.tot[yr + 1] / N.tot[yr]

We can now compute the demographic growth rate, simply defined as the
total growth rate minus the environmental growth rate:

    lam.D[yr] <- lam[yr] - lam.ES[yr]

We have done these computations within the nimble code (in file
`script.r`), but they could also have been done on the R side.

From these different growth rates, we may compute the relative
contributions of environmental variation to variance in the realized
growth. This is done in R (file `contrib.r`) by first extracting `lam`
and `lam.ES` and then computing

    rel.environmental.var = (var(lam.ES)  + cov(lam.ES, lam.D)) / var(lam)

Similarly, the relative contribution of demographic variation is:

    rel.demographic.var = (var(lam.D)  + cov(lam.ES, lam.D)) / var(lam)

which is equal to `1 - rel.environmental.var`.

### Computing contributions of environmental variables to the total growth rate

Above we estimated contributions of environmental variables to the
variance of the environmental growth rate *λ*<sup>*E**S*</sup>. We can
also estimate contributions to the variance of the total growth rate *λ*
that includes both environmental and demographic stochasticity. The
calculation is essentially the same, but we add covariance terms between
each environmental variable and *λ*<sup>*D*</sup>, and then compute
relative contribution to the variance in *λ*.

    # Compute contributions to overall (both demographic and environmental) variance in growth 
    contributions.overall = contributions + grad * apply(lpc.mat, 2, cov, y = lam.D)
    rel.contrib.overall = contributions.overall / var(lam)

The covariance terms are necessary to make the sum of the relative
contributions (approximately) equal to `rel.environmental.var`.

### Decomposing *λ*<sup>*D*</sup>

To decompose *λ*<sup>*D*</sup> we consider all (six) additive components
of its nominator
||*n*<sub>*t* + 1</sub>|| − ||*A*<sub>*t*</sub>*n*<sub>*t*</sub>||.
These consist of the deviations for the number of juvenile recruits from
both age classes, the number of survivors from both age classes, and the
number of immigrants into both age classes. The deviations are the
differences between the number of individuals in the components at time
*t* + 1 (i.e. the terms in the sum ||*n*<sub>*t* + 1</sub>||), and the
expected number of individuals in these components under environmental
variation given the population vector at time *t* (i.e. the terms of
||*A*<sub>*t*</sub>*n*<sub>*t*</sub>||). To compute the contributions we
first extract the number of individuals in each of the components at the
end of each time step, as well as the number of individuals in each
stage at the start of each time step, from the MCMC output (this is done
in the file `contrib.r`):

    N.REC = array(pars[grepl('N.REC\\[', parNames)], dim = c(26,2))[-1,]
    N.SURV = array(pars[grepl('N.SURV\\[', parNames)], dim = c(26,2))[-1,]
    N.IMM = array(pars[grepl('N.IMM\\[', parNames)], dim = c(26,2))[-1,]
    N0 = array(pars[grepl('N\\[', parNames)], dim = c(26,2))[-26,]

We then compute the deviations for each of the components by subtracting
the expected number of individuals.

    D.REC = N.REC- 
            cbind(plogis(rowSums(lpc.nu1)) * r * plogis(rowSums(lpc.psi)) * zeta * plogis(rowSums(lpc.phi0)) * N0[, 1], 
                  plogis(rowSums(lpc.nu2)) * r * plogis(rowSums(lpc.psi)) * zeta * plogis(rowSums(lpc.phi0)) * N0[, 2])
    D.SURV = N.SURV - cbind(plogis(rowSums(lpc.phi1)) * N0[,1], plogis(rowSums(lpc.phi2)) * N0[,2])
    D.IMM = N.IMM - cbind(exp(lpc.omega1) * rowSums(N0), exp(lpc.omega2) * rowSums(N0))

The deviations are then collected into a matrix and we divide by the
number of individuals at the start of each time step to get rates.
Finally, we compute contributions to *λ*<sup>*D*</sup> from the
empirical covariance matrix, and then relative contributions by dividing
with the variance.

    D = cbind(D.REC / rowSums(N0), D.SURV / rowSums(N0), D.IMM / rowSums(N0))
    colnames(D) = paste0(rep(c('REC', 'SURV', 'IMM'), each = 2), 1:2)
      
    contributions.d = as.matrix(rowSums(cov(D)))
    rel.contrib.d = contributions.d / var(lam.D) 
