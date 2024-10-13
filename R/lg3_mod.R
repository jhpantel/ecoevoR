#' Discrete-time evolutionary rescue + competition model for three species
#'
#' Equations for Leslie-gower population growth (discrete-time Lotka-Volterra) models for three species. Population growth is modelled as population mean fitness as per Gomulkiewicz & Holt 1995. These formulas calculate the population's maximum growth rate (with a trait value at the local environmental optimum), the distance from the optimal trait value and evolutionary inertia of trait given phenotypic variation, and the change in population size N from one time step to the next given intra- and interspecific interaction strength.
#'
#' The model is a discrete time, stochastic, quantitative genetic model of evolutionary rescue (Gomulkiewicz & Holt 1995), with some modifications. Populations experience logistic growth: \eqn{N_t+1 = N_t + \lambda N_t ((1-N_t)/K)} and the population growth rate \eqn{\lambda} is the absolute mean fitness of a population in generation *t*. Mean fitness is determined by a quantitative trait *x* as follows:
#' \eqn{\lambda = \hat{W}e^{w+(1-h^2)P / (P+w)(E-x_t)]^2}{2(P+w)}}
#'
#' where \eqn{\hat{W} = W_max sqrt(\frac{w}{P + w})}, *W_max* is the absolute fitness when the phenotype *x* equals the local environmental optimum trait value *E*, *w* is the width of the Gaussian fitness function (which determines the strength of selection), *P* is the width of the distribution of the phenotype *x*, $h^2$ is the heritability of the phenotype *x*, and *x_t* is the population's mean phenotype at time *t*.
#'
#' @param n A vector with the population size for all three species (N1, N2, N3) at time t
#' @param E A scalar with the local environmental optimum trait value E
#' @param tx A vector with the trait value for all three species (x1, x2, x3) at time t
#' @param parms A vector with the elements:
#' \describe{
#' \item{"a11"}{Intraspecific competition coefficient for species 1}
#' \item{"a12"}{Per-capita impact of species 2 on species 1}
#' \item{"a13"}{Per-capita impact of species 3 on species 1}
#' \item{"a21"}{Per-capita impact of species 1 on species 2}
#' \item{"a22"}{Intraspecific competition coefficient for species 2}
#' \item{"a23"}{Per-capita impact of species 3 on species 2}
#' \item{"a31"}{Per-capita impact of species 1 on species 3}
#' \item{"a32"}{Per-capita impact of species 2 on species 3}
#' \item{"a33"}{Intraspecific competition coefficient for species 1}
#' \item{"wmax"}{The maximum growth rate of a species}
#' \item{"w"}{The strength of selection, or the width of the trait-fitness Gaussian distribution}
#' \item{"P"}{The width of the phenotypic distribution}
#' \item{"h21"}{Heritability of trait x in species 1}
#' \item{"h22"}{Heritability of trait x in species 2}
#' \item{"h23"}{Heritability of trait x in species 3}
#' }
#' @return The values of population size, trait values, and population growth rate for species 1, 2, and 3.
#' @export
#'
#' @examples
#' sim_num <- 10
#' m1 <- rand_sim(n=sim_num,x0=c(0.4,0.5,0.6))
#' ## Visualize the results
#' matplot(1:300, m1$pop_record[1,,7:9], type = "l", xlab="time",ylab="N")
lg3_mod <- function (n, E, tx, parms)
{
  with(as.list(parms), {
    ## Eqns 4, 5, & 6 from Gomulkiewicz & Holt 1995.
    ## Calculates popultion growth rate at env. optima and evolutionary inertia of trait given phenotypic variation.
    ## For 1 time step t
    what <- wmax*sqrt(w/(P+w))
    k1 <- (w + (1-h21)*P) / (P+w)
    dx1 <- k1 * (E-tx[1])
    tx1 <- E-dx1
    l1 <- what*exp((-dx1^2) / (2*(P+w)))
    y1 <- (l1 * n[1]) / (1 + (a11 * n[1]) + (a12 * n[2]) + (a13 * n[3]))
    dn1dt <- y1

    k2 <- (w + (1-h22)*P) / (P+w)
    dx2 <- k2 * (E-tx[2])
    tx2 <- E-dx2
    l2 <- what*exp((-dx2^2) / (2*(P+w)))
    y2 <- (l2 * n[2]) / (1 + (a22 * n[2]) + (a21 * n[1]) + (a23 * n[3]))
    dn2dt <- y2

    k3 <- (w + (1-h23)*P) / (P+w)
    dx3 <- k3 * (E-tx[3])
    tx3 <- E-dx3
    l3 <- what*exp((-dx3^2) / (2*(P+w)))
    y3 <- (l3 * n[3]) / (1 + (a33 * n[3]) + (a31 * n[1]) + (a32 * n[2]))
    dn3dt <- y3

    list(c(round(dn1dt), round(dn2dt), round(dn3dt), tx1, tx2, tx3, l1, l2, l3))
  })
}
