######################################
###            METADATA            ###
######################################
# Version: 1.0
#
# Author: Oskar Hagen
#
# Date: 26.10.2020
#
# Landscape: SouthAmerica
#
# Publications: R-package gen3sis
#
# Description: Example config used at the introduction vignette and similar to case study global configs in Hagen et al. 2020.
# O. Hagen, B. Flück, F. Fopp, J.S. Cabral, F. Hartig, M. Pontarp, T.F. Rangel, L. Pellissier. gen3sis: The GENeral Engine for Eco-Evolutionary SImulationS on the origins of biodiversity.
######################################


######################################
###         General settings       ###
######################################

# set the random seed for the simulation
random_seed = 24

# set the starting time step or leave NA to use the earliest/highest time-step
start_time = 40

# set the end time step or leave as NA to use the latest/lowest time-step (0)
end_time = NA

# maximum total number of species in the simulation before it is aborted
max_number_of_species = 50000

# maximum number of species within one cell before the simulation is aborted
max_number_of_coexisting_species = 10000

# a list of traits to include with each species
trait_names = c("temp",  "dispersal")

# ranges to scale the input environments with:
environmental_ranges = list("temp" = c(-45, 55), "area"=c(2361.5, 12923.4), "arid"=c(1,0.5))

######################################
###            Observer            ###
######################################

# a place to inspect the internal state of the simulation and collect additional information if desired
end_of_timestep_observer = function(data, vars, config){
  save_species()
  #plot_richness(data$all_species, data$landscape)
  # example 1 plot over simulation
    # par(mfrow=c(2,3))
    # plot_raster_single(data$landscape$environment[,"temp"], data$landscape, "temp", NA)
    # plot_raster_single(data$landscape$environment[,"arid"], data$landscape, "arid", NA)
    # plot_raster_single(data$landscape$environment[,"area"], data$landscape, "area", NA)
    # plot_richness(data$all_species, data$landscape)
    # plot_species_presence(data$all_species[[1]], data$landscape)
    # plot(0,type='n',axes=FALSE,ann=FALSE)
    # mtext("STATUS",1)
  # example 2 plot over simulations saving plots
    # plot_richness(data$all_species, data$landscape)
    # plot_landscape(data$landscape)

}

######################################
###         Initialization         ###
######################################

# the initial abundance of a newly colonized cell, both during setup and later when colonizing a cell during the dispersal
initial_abundance = 1

#defines the initial species traits and ranges
#place species within rectangle, our case entire globe

create_ancestor_species <- function(landscape, config) {
  range <- c(-95, -24, -68, 13)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
    co[, "x"] <= range[2] &
    co[, "y"] >= range[3] &
    co[, "y"] <= range[4]

  new_species <- list()
  for(i in 1:10){
    initial_cells <- rownames(co)[selection]
    initial_cells <- sample(initial_cells, 1)
    new_species[[i]] <- create_species(initial_cells, config)
    #set local adaptation to max optimal temp equals local temp
    new_species[[i]]$traits[ , "temp"] <- landscape$environment[initial_cells,"temp"]
    new_species[[i]]$traits[ , "dispersal"] <- 1
    #plot_species_presence(landscape, species=new_species[[i]])
  }

  return(new_species)
}



######################################
###             Dispersal          ###
######################################

# returns n dispersal values
get_dispersal_values <- function(n, species, landscape, config) {
  traits <- species[["traits"]]
  values <- traits[, "dispersal"]
  valmat <- matrix(NA,nrow=length(values),ncol=(n / length(values)))

  for(i in 1:length(values)){
    valmat[i,] <- rweibull(n=(n / length(values)), shape = 1.5, scale = 133 + values[i])
  }
  return(as.vector(valmat))
}

######################################
###          Speciation            ###
######################################

# threshold for genetic distance after which a speciation event takes place
divergence_threshold = 2 #this is 1Myrs

# factor by which the divergence is increased between geographically isolated population
# can also be a matrix between the different population clusters
get_divergence_factor <- function(species, cluster_indices, landscape, config) {

  return(1)
}


######################################
###            Evolution           ###
######################################

# mutate the traits of a species and return the new traits matrix
apply_evolution <- function(species, cluster_indices, landscape, config) {

  trait_evolutionary_power <- 0.001
  trait_evolutionary_power_fast <- 0.01
  traits <- species[["traits"]]
  cells <- rownames(traits)
  #homogenize trait based on abundance
  for(cluster_index in unique(cluster_indices)){
    cells_cluster <- cells[which(cluster_indices == cluster_index)]
    mean_abd <- mean(species$abundance[cells_cluster])
    weight_abd <- species$abundance[cells_cluster]/mean_abd
    traits[cells_cluster, "temp"] <- mean(traits[cells_cluster, "temp"]*weight_abd)
    traits[cells_cluster, "dispersal"] <- mean(traits[cells_cluster, "dispersal"]*weight_abd)
  }
  #mutations
  # vcv is the variance-covariance matrix
  tcor <- -.25
  vcv <- matrix(c((trait_evolutionary_power^2),(tcor*(trait_evolutionary_power)*(trait_evolutionary_power_fast)),(tcor*(trait_evolutionary_power)*(trait_evolutionary_power_fast)),(trait_evolutionary_power_fast^2)),2,2)

  if(length(traits) == 2){
    mutation_deltas <- mvrnorm(length(traits[1]),mu=c(0,0),Sigma=vcv)
    traits[, "temp"] <- traits[, "temp"] + mutation_deltas[1]
    traits[, "dispersal"] <- traits[, "dispersal"] + mutation_deltas[2]
  } else {
    mutation_deltas <- mvrnorm(length(traits[, "temp"]),mu=c(0,0),Sigma=vcv)
    traits[, "temp"] <- traits[, "temp"] + mutation_deltas[,1]
    traits[, "dispersal"] <- traits[, "dispersal"] + mutation_deltas[,2]
  }
  return(traits)
}


######################################
###             Ecology            ###
######################################

# called for every cell with all occurring species, this function calculates the who survives in the current cells
# returns a vector of abundances
# set the abundance to 0 for every species supposed to die

apply_ecology <- function(abundance, traits, landscape, config) {
  abundance_scale = 10
  abundance_threshold = 1
  #abundance threshold
  survive <- abundance>=abundance_threshold
  abundance[!survive] <- 0
  abundance <- (( 1-abs( traits[, "temp"] - landscape[, "temp"]))*abundance_scale)*as.numeric(survive)
  #abundance threshold
  abundance[abundance<abundance_threshold] <- 0
  k <- ((landscape[,"area"]*(landscape[,"arid"]+0.1)*(landscape[,"temp"]+0.1))*abundance_scale^2)
  total_ab <- sum(abundance)
  subtract <- total_ab-k
  if (subtract > 0) {
    # print(paste("should:", k, "is:", total_ab, "DIFF:", round(subtract,0) ))
    while (total_ab>k){
      alive <- abundance>0
      loose <- sample(1:length(abundance[alive]),1)
      abundance[alive][loose] <- abundance[alive][loose]-1
      total_ab <- sum(abundance)
    }
    #set negative abundances to zero
    abundance[!alive] <- 0
  }

  return(abundance)
}
