#' Function with ODE system for host-pathogen eco-evolutionary model
#'
#' @param t A vector of time steps for which the differential equations will be evaluated across
#' @param x A vector of starting values for S, I, and alpha
#' @param params A vector of values for the model parameters B, μ, V, γ, and c
#'
#' @return The values of the derivatives in the ODE system at time t.
#' @export
#'
#' @examples
#' parms <- c(B = (10^7), mu = 1, V = 0, gam = 2, c = 3*10^-7)
#' times <- seq(from=0,to=40,by=1)
#' xstart <- c(S=10^7,I=1,alpha=4)
#' out <- hostpath_evol(times,xstart,parms)
#'
lk_evol <- function (t, x, params) {
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
