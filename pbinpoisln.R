pbinpoisln <- function(n, w, pd, lam=list(meanlog, sdlog))
{
  ## Purpose:
  ##
  ## Calculate the probability of acceptance under the Binomial-Poisson-Lognormal
  ## distribution assuming a zero tolerance approach.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## n: number of samples drawn from the lot
  ## w: weight/surface area of each sample
  ## pd: proportion of the lot that is contaminated
  ## lam: is the rate of contamination (per g / per cm2 basis). It consists
  ##      of the two parameters from the Lognormal distribution, namely meanlog
  ##      and sdlog - note that these are on the ln rather than log10 scale!
  ##      Follows the approach of the various lnorm functions.
  ## ----------------------------------------------------------------------

  f1 <- function(x, mu, sd, n, pd, w){
    dlnorm(x, meanlog=mu, sdlog=sd)*(1-pd+pd*exp(-w*x))^n
  }

  ## Bounds for the integration
  l.lo <- qlnorm(0.001, meanlog=lam$meanlog, sdlog=lam$sdlog)
  l.hi <- qlnorm(0.999, meanlog=lam$meanlog, sdlog=lam$sdlog)
  
  fx0 <- try(integrate(f1, lower=l.lo, upper=l.hi, mu=lam$meanlog, sd=lam$sdlog,
                   n=n, pd=pd, w=w, subdivisions=1000), silent=TRUE)

  ## If this doesn't work then the upper limit should be set to infinity.
  if(class(fx0)=="try-error")
    fx0 <- try(integrate(f1, lower=l.lo, upper=Inf, mu=lam$meanlog,
                         sd=lam$sdlog, n=n, pd=pd, w=w, subdivisions=1000)) 

  ## If this doesn't work then the lower limit should be set to zero.
  if(class(fx0)=="try-error")
    fx0 <- try(integrate(f1, lower=0, upper=Inf, mu=lam$meanlog,
                         sd=lam$sdlog, n=n, pd=pd, w=w, subdivisions=1000)) 

  ## If this still doesn't work then return NA to signal a problem.
  if(class(fx0)=="try-error")
    fx0 <- list(value=NA)

  ## Return the solution of the integral and the arithmetic mean
  return(list(fx0=fx0$value, mean.lambda=exp(lam$meanlog+0.5*lam$sdlog^2)))
}
