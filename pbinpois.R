pbinpois <- function(n, w, pd, lam)
{
  ## Purpose:
  ## 
  ## Calculate the probability of acceptance under the Binomial-Poisson
  ## (Habraken) approach assuming that zero tolerance is used.
  ## 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## 
  ## n: number of samples drawn from the lot
  ## w: weight/surface area of each sample
  ## pd: proportion of the lot that is contaminated
  ## lam: is the rate of contamination (per g / per cm2 basis) - fixed
  ## 
  ## ----------------------------------------------------------------------
  ## Author: Andreas Kiermeier, Date: 22 Feb 2012, 13:06
  return ( (1 - pd + pd*exp(-w*lam))^n )
}
