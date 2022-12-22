#' Simulations of 3-species Lotka-Volterra model with data prepared for ABC analysis
#'
#' A dataset containing the summary statistics and parameter values for two sets of models (without and with trait evolution, using the beginning- and end-point trait data or only the end-point data)
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{abc_sim_parm}{An nx2 by 15 matrix, with the number of simulations for Model 1 (no evolution) and Model 2 (with evolution) as rows and the 15 model parameter values as columns. Each row gives the parameter values chosen in *rand_sim* for that simulation iteration.}
#'   \item{abc_sim_parm_x0}{Same as abc_sim_parm but for the model with only the end-point trait data}
#'   \item{abc_sim_summ}{An nx2 by 69 matrix, with the number of simulations for Model 1 and Model 2 as rows and the 69 summary statistics as columns. The initial trait value was included when rand_sim was run, so the only included summary statistic is the end-point trait value. The other summary statistics are the population size values for each species at each of the intervals recorded (in this example days <- seq(1,300,by=14)).}
#'   \item{abc_sim_summ_x0}{Same as for abc_sim_summ. In this instance the initial trait value was randomly chosen in the simulation model during rand_sim, and again only the end-point trait value is included as a summary statistic.}
#'   \item{abc_sim_mods}{A length nx2 vector with elements that give the model (no_evol or evol) that was used to generate the data for all rows in the list elements for each simulation.}
#' }
#' @source See ecoevo-ABC section '2 Simulation'
"abc_lv"

#' Random forest feature selection, parameter tuning, and model classification for simulations of 3-species Lotka-Volterra model
#'
#' A dataset containing the results of a boruta feature selection from simulated data of 3 species population size and trait values over time (Boruta::Boruta), parameter tuning (randomForest:tuneRF), and model classification (randomForest:randomForest). Simulated data comes from a 3-species Lotka-Volterra model with a quantitative genetic model for potential evolution of population growth rate (see vignette ecoevo-ABC). Hypotheses for model classification of simulated observed data are No Evolution(h2 = 0) and Evolution (h2 > 0).
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{bor}{A list with the elements resulting from feature selection algorithm Boruta::Boruta}
#'   \item{tune}{A list with the elements resulting from parameter tuning randomForest:tuneRF}
#'   \item{rf}{A list with the elements resulting from parameter tuning randomForest:randomForest}
#' }
#' @source See ecoevo-ABC section '2 Simulation'
"rf_lv"
