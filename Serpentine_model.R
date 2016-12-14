# "Serpentine Model" - basic start.

# initialize grids for two different plants

#create the distance matrix and the coordinate matrix...

{  
  coord = matrix(0,nrow = 2, ncol = mn) # coord is where coordinate information for calculating distances will be stored. 
  coord[1,] = rep(1:m, each = n) # rep functions generate the pattern that describes the coordinates in the appropriate order for the coordinate matrix
  coord[2,] = rep(1:n, times = m)
  dimnames(coord) = list(c("row", "column"))
  dist = matrix(0, nrow = mn, ncol= mn) # establish a matrix for storing distance values - each row represents an individual with columns listing distance from that individual to other individuals by index.

  # rows are all seeds from a certain point
  # columns are are the points to which they distribute
  for(j in  1:mn){   # from every point
    for (i in 1:mn){ # to every other point
      dist[j,i] = sqrt((coord["row",i]-coord["row",j])^2+(coord["column",i]-coord["column",j])^2) # distance = sqrt(a^2 + b^2), where a and b are differences in x and y direction respectively from a point of interest
   }}
}

# create all the other grids which will store information on the plants.
{
grids <- initialize_grids(rows = m, cols = n, positionx = xpos, positiony = ypos, 
                          obspos = obspos_coords, obsintens = .5) # create plant landscape according to the parameters above.
plant_grid <- grids[[1]]
plant_age <- grids[[2]]
new_plants <- grids[[3]]
obs <- grids[[4]]
nspecies <- grids[[5]]
# initial plant ages
plant_age = plant_age + 1 * (plant_grid > 0)
population = matrix(0, nrow = generations, ncol = nspecies)
}
    
for (y in 1:generations){
  
  # run reproduction and veggrowth for each plant
  # output: locations of "new plants" for each species

  occupied_cells = apply(plant_grid, MARGIN = c(1, 2), FUN = sum)
  
for (s in 1:nspecies){

  plant_grid[,,s] <- obs * plant_grid[,,s]
  
  matrix_reprod <- matrix_reproduction_mp(plant_cells = plant_grid[,,s], dist_cells = dist, 
                                          slope = slope[s], b = b[s])
  matrix_veg <- veggrowth_mp(plant_cells = plant_grid[,,s], dist_cells = dist,
                             vegchance = vegchance[s])
      # now we have 2 grids, one with new plants developed through veg, the other through seeds.

  new_plants[,,s] = pmin(1, (matrix_reprod + matrix_veg)) # adds the vegetative and seed plants together
                        
  new_plants[,,s] }

# propegules compete
  
  all_new = apply(new_plants, MARGIN = c(1,2), FUN = sum) 
  points = which(all_new > 1, arr.ind = TRUE)
      # points is a matrix with two columns - rows and collumns
      # gives locations of points in the array where there are overlaps.
  if (nrow(points) > 0){

  for (i in 1:nrow(points)){
    location = new_plants[points[i,1], points[i,2],] #specifies location of cell, spits out all values in those cells
    if (obs[points[i,1], points[i,2]] < 1){
      weight = tolerence * location # weights tolerence based on occupied/unoccupied - unoccupied  cells have zero chance of being selected for   
    }else{
      weight = comp * location #weights competition based on a value rating of competative ability
    }
    winner = sample(nspecies, size = 1, prob = weight) # winner chosen among all species - unoccupied cells weighted zero, all others weighted by tolerence.
    new_plants[points[i,1], points[i,2],] = 0 # removes all plants at that location
    new_plants[points[i,1], points[i,2],winner] = 1 # establishes winner plant species at that location
    
   }}
  
# verify that new_plants is all what is expected

new_plants_test = apply(new_plants, MARGIN = c(1,2), FUN = sum)
errors = which(new_plants_test > 1)
stopifnot(length(errors) == 0)

# revise plant grid to include the newest plants. 
for (s in 1:nspecies){
  new_plants[,,s] = new_plants[,,s]*(1 - occupied_cells) # new plants minus new plants in places where other plants win minus places which are occupied
  plant_grid[,,s] = plant_grid[,,s] + new_plants[,,s] # add new plant to existing plant grids for each species
  }

# all plants age up one year -- new plant age = 1  

plant_age = plant_age + 1 * (plant_grid > 0) 
plant_grid = (plant_grid > 0) * 1     # resets all plants to 1 or zero, no in-betweens.

dead = mortality(plant_age, nspecies = nspecies)
plant_age[dead == 1] = 0
plant_grid = plant_grid - dead * plant_grid

# additional statistics on results:
population[y,] <- apply(plant_grid, MARGIN = 3, FUN = sum) 
    
}