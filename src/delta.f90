FUNCTION delta(x, sigma)
  !
  USE constants, only : dp, sqrtpi
  !
  IMPLICIT NONE
  !
  real(dp) delta
  real(dp) x, sigma
  !
  delta=exp(-(x/sigma)**2)/(sigma*sqrtpi)
  !
END FUNCTION
