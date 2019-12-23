#### Script containing set up for all demographic functional forms

# COLLAPSE ALL FOR EASIER VIEWING

#### FUNDAMENTAL FUNCTIONS ####

# Survival - minimal form so most versatile
# covariates all included as 'extra'
S.fun <- function(z, int, slope.z, slope.PS, extra, Ns){ # extra = sum of all other slopes*est
  u <- exp(int+(z*slope.z)+(Ns*slope.PS)+extra)
  u/(1+u)
}

# Development
D.fun <- function(z, zz, mu.int, mu.zs, mu.ps, extra, sigma.d, Ns) {
  mu.z <- mu.int+(mu.zs*z)+(mu.ps*Ns)+extra
  sigma.d2 <- sigma.d^2
  temp1 <- sqrt(2*pi)*sigma.d
  temp2 <- ((zz-mu.z)^2)/(2*sigma.d2)
  return(exp(-temp2)/temp1)
}

# Recruitment
R.fun <- function(z, int, slope.z, slope.ps, extra, Ns){
  u <- exp(int+(z*slope.z)+(Ns*slope.ps)+extra)
  u
}

# Inheritance
H.fun <- function(z, zz, mu.int, mu.N, sigma.h, extra, Ns) {
  mu.z <- mu.int+(mu.N*Ns)+extra
  sigma.h2 <- sqrt(sigma.h)
  temp1 <- sqrt(2*pi)*sigma.h2
  temp2 <- ((zz-mu.z)^2)/(2*sigma.h)
  return(exp(-temp2)/temp1)
}
# standard IPM below, quantgen above
H.fun.PO <- function(z, zz, mu.int, mu.N, mu.mZ, sigma.h, extra, Ns) {
  mu.z <- mu.int+(mu.N*Ns)+extra
  sigma.h2 <- sigma.h^2
  temp1 <- sqrt(2*pi)*sigma.h
  temp2 <- ((zz-mu.z)^2)/(2*sigma.h2)
  return(exp(-temp2)/temp1)
}


# Half fall (caterpillar timing)
HF.fun <- function(int, slope.ST, slope.WT, slope.SP, slope.WP, ST, WT, WP, SP, B) {
  u <- int+(ST*slope.ST)+(WT*slope.WT)+(WP*slope.WP)+(SP*slope.SP)+B
  u
}
