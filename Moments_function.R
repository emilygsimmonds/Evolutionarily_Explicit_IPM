#### Script containing function to extract four moments from model simulation outputs

# calculates moments of the results
moments.fun <- function(x,z){
  mu <- sum(z*x)/sum(x) # calculates mean
  va <- sum(z*z*x)/sum(x)-mu^2 # calculates variance
  sk <- (sum(z*z*z*x)/sum(x)-3*mu*va-mu^3)/(sqrt(va)^3) # calculates skew
  m4 <- sum(z*z*z*z*x)/sum(x)
  m3 <- sum(z*z*z*x)/sum(x)
  m2 <- sum(z*z*x)/sum(x)
  ku <- (m4-4*mu*m3+6*(mu^2)*m2-3*(mu^4))/(va^2)-3
  res <- c(mu,va,sk,ku,m2,m3,m4)
  names(res) <- c('Mean','Variance','Skew','Kurtosis','Moment 2','Moment 3','Moment 4')
  return(res)
}
