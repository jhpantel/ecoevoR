#' Run a simulation with parameter values randomly chosen from prior distributions (currently for 3-species Leslie-Gower model)
#'
#' @param N0 A vector with the initial population size for all three species ($N_{1,0}$, $N_{2,0}$, $N_{3,0}$)
#' @param x0 A vector with the initial trait value *x* for all three species ($x_{1,0}$, $x_{2,0}$, $x_{3,0}$)
#' @param E A scalar with the local environmental optimum trait value E
#' @param time_points A numeric value of the number of time steps for the simulation
#' @param n A numeric value of the number of simulations to run
#' @param evol A boolean where FALSE sets the heritability for all species equal to 0 and TRUE draws a heritability value for all species from a uniform distribution between 0-1.
#' @return A list with the elements:
#' \itemize{
#' \item{"pop_record"}{A three-dimensional array that returns the following for all simulation iterations (n) at all time points:
#' \describe{
#'   \item{"x1,x2,x3"}{The population size of all three species at time t}
#'   \item{"tx1,tx2,tx3"}{The trait value of all three species at time t}
#'   \item{"y1,y2,y3"}{The population size of all three species at time t+1}
#'   \item{"yt1,yt2,yt3"}{The trait value of all three species at time t+1}
#'   \item{"E"}{The local environmental optimum trait value at time t+1}
#'   \item{"l1,l2,l3"}{The per-capita  growth rate of all three species from time t to t+1}
#' }}
#' \item{"parm_record"}{An n-by-18 matrix with the values for all model parameters used for each of the n simulation iterations.}
#' @export
#'
#' @examples
#' sim_num <- 10
#' m1 <- rand_sim(n=sim_num,x0=c(0.4,0.5,0.6))
#' ## Visualize the results
#' matplot(1:300, m1$pop_record[1,,7:9], type = "l", xlab="time",ylab="N")
rand_sim <- function(N0=c(10, 10, 10),x0=NULL,E=.1,time_points=300,n,evol=FALSE){
  # Record of all n simulations
  pop_record <- array(NA,dim=c(n,time_points,17))
  parm_record <- matrix(NA,nrow=n,ncol=15,dimnames=list(NULL,c("a11", "a12", "a13", "a21", "a22", "a23", "a31", "a32", "a33", "wmax", "w", "P", "h21", "h22", "h23")))

  ## Provide fixed values and draw random values from priors
  # Fixed values
  wmax <- 2
  # N0, x0, E, time_points
  if(is.null(x0)) {
    x0 <- runif(3)
  }

  for(g in 1:n){
    print(paste("Running sim ",g," out of ",n,sep=""))
    # Draws from priors
    x0 <- runif(3)
    a11 <- rbeta(1,0.25,10)
    a12 <- rbeta(1,0.25,10)
    a13 <- rbeta(1,0.25,10)
    a21 <- rbeta(1,0.25,10)
    a22 <- rbeta(1,0.25,10)
    a23 <- rbeta(1,0.25,10)
    a31 <- rbeta(1,0.25,10)
    a32 <- rbeta(1,0.25,10)
    a33 <- rbeta(1,0.25,10)

    w <- rgamma(1,.25,.25)
    #w <- 0.5
    P <- 0.05*w
    if (evol==FALSE){
      h21 = h22 = h23 = 0
    } else {
      h21 <- runif(1,0,1)
      h22 <- runif(1,0,1)
      h23 <- runif(1,0,1)
    }

    parms <- c(a11 = a11, a12 = a12, a13 = a13,
               a21 = a21, a22 = a22, a23 = a23,
               a31 = a31, a32 = a32, a33 = a33,
               wmax = wmax, w = w, P = P,
               h21 = h21, h22 = h22, h23 = h23)

    pop_record[g,,] <- LV_evol(N0,x0,E,time_points,parms)
    parm_record[g,] <- parms
  }
  return(list(pop_record=pop_record,parm_record=parm_record))
}
