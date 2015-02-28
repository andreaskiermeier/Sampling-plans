phyppois <- function(n, w, N, D, lam)
{
  ## Purpose:
  ## 
  ## Calculate the probability of acceptance under the Hypergeometric-Poisson
  ## (Habraken) approach assuming that zero tolerance is used.
  ## 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## 
  ## n: number of samples drawn from the lot
  ## w: weight/surface area of each sample
  ## N: Number of primary sampling units in the lot (consists of units which are sub-sampled)
  ## D: Number of 'defective'/contaminated primary sampling units in the lot.
  ## lam: is the rate of contamination (per g / per cm2 basis) - fixed across
  ##      primary sampling units
  ## 
  ## ----------------------------------------------------------------------
  ## Author: Andreas Kiermeier, Date: 30 May 2012, 12:08
  i <- 0:n
  return(sum(dhyper(i, D, N-D, n)*exp(-i*w*lam)))
}
