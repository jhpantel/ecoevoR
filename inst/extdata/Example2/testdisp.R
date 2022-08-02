### Recreate dispersal process for some point in time
### JH Pantel, 02/08/2022
presence_spi_ti <- names(species_t_40[[5]]["abundance"]$abundance)
all_cells <- rownames(landscape_t_40$coordinates)
free_cells <- all_cells[!(all_cells %in% presence_spi_ti)]
num_draws <- length(free_cells) * length(presence_spi_ti)
# To get dispersal values
traits <- species_t_40[[5]]["traits"]$traits
values <- traits[, "dispersal"]
# To get distance matrix
neighbour_file <- paste0(datapath,"/created_landscapes/distances_local/distances_local_65.rds")
con <- gzfile(neighbour_file)
distance_neighbours <- readRDS(file=neighbour_file)
habitable_cells <- as.integer(rownames(landscape_t_40$coordinates))
num_cells <- nrow(distance_neighbours)

sourceCpp(paste0(datapath,'/distances.cpp'))
distance_matrix <- get_distance_matrix(habitable_cells,num_cells,distance_neighbours@p,distance_neighbours@i,distance_neighbours@x,Inf)

# To use dispersal trait values
#get_dispersal_values <- function(n, species, landscape, config) {
#  traits <- species[["traits"]]
#  values <- traits[, "dispersal"]
  valmat <- matrix(NA,nrow=length(values),ncol=(num_draws / length(values)))

  for(i in 1:length(values)){
    valmat[i,] <- rweibull(n=(num_draws / length(values)), shape = 1.5, scale = 133 + values[i])
  }
#  return(as.vector(valmat))
#}
r_disp <- as.vector(valmat)

geo_disp <- distance_matrix[presence_spi_ti, free_cells, drop=FALSE] #lines mark where they are present, cols the possible suitable sites
# Short test to be sure about comparison
#blah1 <- matrix(data=c(rep(1,6),rep(2,6),rep(3,6),rep(4,6)),nrow = 4,ncol = 6,byrow = TRUE)
#blah2 <- c(rep(c(1:4),times=6))
#blah3 <- blah1 <= blah2

geo_disp <- geo_disp <= r_disp

## What is needed for the config dispersal function
# returns n dispersal values
get_dispersal_values <- function(n, species, landscape, config) {
  traits <- species[["traits"]]
  values <- traits[, "dispersal"]
  valmat <- matrix(NA,nrow=length(values),ncol=(num_draws / length(values)))

  for(i in 1:length(values)){
    valmat[i,] <- rweibull(n=(num_draws / length(values)), shape = 1.5, scale = 133 + values[i])
  }
  return(as.vector(valmat))
}

r_disp <- get_dispersal_values(num_draws, species, landscape, config)
