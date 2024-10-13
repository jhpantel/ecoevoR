#' Simulate three-species discrete-time growth with competition and evolving growth rate
#'
#' Simulate competition, growth, and trait evolution in an eco-evolutionary version of a discrete-time Leslie-Gower model. The function calls the formulas in lg3_mod over the time intervals for which population growth is calculated.
#'
#' @param N0 A vector with the initial population size for all three species (\eqn{N_{1,0}$, $N_{2,0}$, $N_{3,0}})
#' @param x0 A vector with the initial trait value *x* for all three species (\eqn{x_{1,0}$, $x_{2,0}$, $x_{3,0}})
#' @param E A scalar with the local environmental optimum trait value E
#' @param time_points A numeric value of the number of time steps for the simulation
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
#'
#' @return A dataframe with the elements:
#' \describe{
#' \item{"x1,x2,x3"}{The population size of all three species at time t}
#' \item{"tx1,tx2,tx3"}{The trait value of all three species at time t}
#' \item{"y1,y2,y3"}{The population size of all three species at time t+1}
#' \item{"yt1,yt2,yt3"}{The trait value of all three species at time t+1}
#' \item{"E"}{The local environmental optimum trait value at time t+1}
#' \item{"l1,l2,l3"}{The per-capita  growth rate of all three species from time t to t+1}
#' }
#'
#' @export
#'
#' @examples
#' sim_num <- 10
#' m1 <- rand_sim(n=sim_num,x0=c(0.4,0.5,0.6))
#' ## Visualize the results
#' matplot(1:300, m1$pop_record[1,,7:9], type = "l", xlab="time",ylab="N")
LV_evol <- function(N0,x0,E,time_points,parms)
{
  ## Set up data object structure
  pop <- matrix(NA,nrow=time_points,ncol=17)
  colnames(pop) <- c("x1","x2","x3","tx1","tx2","tx3","y1","y2","y3","yt1","yt2","yt3","E","time","l1","l2","l3")
  ## Initial values
  pop[1,1:3] <- N0
  pop[1,4:6] <- x0
  pop[1,13] <- E
  pop[1,14] <- 1
  pop[1,c(7:12,15:17)] <- lg3_mod(pop[1,1:3],pop[1,13],pop[1,4:6],parms)[[1]]

  ## Time series
  for(i in 2:time_points){
    pop[i,1:3] <- pop[i-1,7:9]
    pop[i,4:6] <- pop[i-1,10:12]
    pop[i,13] <- E
    pop[i,14] <- i
    pop[i,c(7:12,15:17)] <- lg3_mod(pop[i,1:3],pop[i,13],pop[i,4:6],parms)[[1]]
  }
  return(pop)
}
