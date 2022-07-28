#' Simulate host-pathogen eco-evolutionary model
#'
#' Simulate host-pathogen dynamics under the eco-evolutionary model described in Luo & Koelle (2013). In this model, *S* and *I* are the susceptible and infected host population size, *B* is a constant birth rate, μ is the per capita death rate, α is the disease-induced mortality rate (virulence), and β is the transmission rate. The transmission-virulence trade-off is modeled by $\beta(\alpha) = c\alpha^{\frac{1}{\gamma}}$, where *c* is a positive constant, γ indicates the magnitude of the transmission-virulence tradeoff, and *V* is the additive genetic variance.
#' The function calls ode (deSolve) to solve the system of ordinary differential equations given in function "lk_evol".
#'
#' @param t A vector of time steps for which the differential equations will be evaluated across
#' @param x A vector of starting values for S, I, and alpha
#' @param params A vector of values for the model parameters B, μ, V, γ, and c
#'
#' @return A dataframe with the elements:
#' \itemize{
#' \item{"time"}{Numeric, time values from parameter t}
#' \item{"S"}{Numeric, a time series of parameter S (susceptible host population size)}
#' \item{"I"}{Numeric, a time series of parameter I (infected host population size)}
#' \item{"alpha"}{Numeric, a time series of parameter α (the virulence of the pathogen)}
#' }
#' @export
#' @import deSolve
#'
#' @examples
#' parms <- c(B = (10^7), mu = 1, V = 0, gam = 2, c = 3*10^-7)
#' times <- seq(from=0,to=40,by=1)
#' xstart <- c(S=10^7,I=1,alpha=4)
#' out <- hostpath_evol(times,xstart,parms)
#'
hostpath_evol <- function(t, x, params){
  out <- as.data.frame(deSolve::ode(y=x,times=t,func=lk_evol,parms=params))
}
