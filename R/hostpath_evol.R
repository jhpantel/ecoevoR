hostpath_evol <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  alpha <- x[3]
  ## now extract the parameters
  B <- params["B"]
  mu <- params["mu"]
  V <- params["V"]
  gam <- params["gam"]
  c <- params["c"]
  beta <- c*(alpha^(1/gam))

  ## now code the model equations
  dSdt <- B - beta*S*I - mu*S
  dIdt <- beta*S*I - (mu + alpha)*I
  dalphadt <- V*((S*(c/gam)*alpha^((1/gam)-1))-1)
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dalphadt)
  ## return result as a list
  list(dxdt)
}
