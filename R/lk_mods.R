#' @name lk_mods
#' @aliases lk_ecoXevo
#' @aliases lk_eco_evo
#' @title Functions with ODE system for host-pathogen eco-evolutionary model
#'
#' @description lk_ecoXevo and lk_eco_evo are functions with ordinary differential equations for an eco-evolutionary host pathogen model. They differ in whether virulence evolution evolves and impacts host population size S - in lk_ecoXevo these are coupled and in lk_eco_evo virulence evolution is not dependent on S.
#'
#' @param t A vector of time steps for which the differential equations will be evaluated across
#' @param x A vector of starting values for S, I, and alpha
#' @param params A vector of values for the model parameters B, μ, V, γ, and c
#'
#' @return The values of the derivatives in the ODE system at time t.
#'
#' @examples
#' parms <- c(B = (10^7), mu = 1, V = 0, gam = 2, c = 3*10^-7)
#' times <- seq(from=0,to=40,by=1)
#' xstart <- c(S=10^7,I=1,alpha=4)
#' out <- hostpath_evol(times,xstart,parms,method="ecoXevo")
#'
#' @rdname lk_mods
#' @export
lk_ecoXevo <- function (t, x, params) {
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

#' @rdname lk_mods
#' @export
lk_eco_evo <- function (t, x, params) {
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
  dalphadt <- V*(((5*10^6)*(c/gam)*alpha^((1/gam)-1))-1)
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dalphadt)
  ## return result as a list!
  list(dxdt)
}
