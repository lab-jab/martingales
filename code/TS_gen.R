library(fracdiff)

epsilon <- function(size, st, df) {
  set.seed(NULL)
  if(st == FALSE) {
    y <- rnorm(size, 0, 1)
  } else {
    y <- rt(size, df)
  }
  return(y)
}

GARCH <- function(size, omega, alpha, beta, st) {
  set.seed(NULL)
  e <- epsilon(size, st, 3)
  h <- 1
  y <- rep(0,size)
  for(t in 2:size) {
    h <- omega + alpha * e[t-1]**2 + beta * h 
    y[t] <- sqrt(h) * e[t]
  }
  return(y)
}

SV <- function(size, st) {
  set.seed(NULL)
  e <- epsilon(size, st, 3)
  h <- rnorm(1, 0, 1)
  y <- rep(0,size)
  for(t in 1:size) {
    y[t] <- exp(0.5 * h) * e[t]
    h <- 0.95 * h + rnorm(1,0,0.1)
  }
  return(y)
}

AR <- function(size, phi, omega, alpha, beta, st, inov) {
  set.seed(NULL)
  if(inov == 2) {
    u <- SV(size, st)
  } else {
    u <- GARCH(size, omega, alpha, beta, st)
  }
  y <- rep(0, size)
  for(t in 2:size) {
    y[t] <- phi * y[t-1] + u[t]
  }
  return(y)
}

ARFIMA <- function(size, omega, alpha, beta, st, inov) {
  set.seed(NULL)
  if(inov == 2) {
    u <- SV(size, st)
  } else {
    u <- GARCH(size, omega, alpha, beta, st)
  }
  return(fracdiff.sim(size, ar = NULL, ma = NULL, 0.1, innov = u)$series)
}
