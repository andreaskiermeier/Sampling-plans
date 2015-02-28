pbinpoisg <- function(n, w, pd, lam=list(a, b))
{
  ## Purpose:
  ##
  ## Calculate the probability of acceptance under the Binomial-Poisson-Gamma
  ## distribution assuming a zero tolerance approach.
  ## 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## 
  ## n: number of samples drawn from the lot
  ## w: weight/surface area of each sample
  ## pd: proportion of the lot that is contaminated
  ## lam: is the rate of contamination (per g / per cm2 basis). It consists
  ##      of the two parameters from the Gamma distribution, namely a=shape and
  ##      b=rate.
  ## 
  ## ----------------------------------------------------------------------
  ## Author: Andreas Kiermeier, Date: 22 Feb 2012, 13:12
  y <- 0:n
  fx0 <- (lam$b^lam$a)*sum(dbinom(y, n, pd)*(w*y+lam$b)^(-lam$a))
  return(list(fx0=fx0, mean.lambda=lam$a/lam$b))
}
