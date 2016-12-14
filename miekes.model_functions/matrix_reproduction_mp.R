# matrix reproduction recovery

matrix_reproduction_mp = function(plant_cells, dist_cells, slope = (-3/5), b = 3)  
  # slope is slope in equation describing the relationship between seeds distributed and distance
  # y-intercept in equation describing the relationship between seeds distributed and distance
{ #how many seeds go to that distance?
  rows = dim(plant_cells)[1]
  cols = dim(plant_cells)[2]
  # calculates the number of seeds distrubuted from each cell to every other cell
  seed_dist = dist_cells * (slope) + b  # calculates seeds distributed to a cell from every other cell based on distance value              # this equation expresses the number N of new seeds based on distance from propagule
  seed_dist[which(seed_dist < 0)] = 0   # this turns any negative seed values to zero (more realistic)
  if (sum(plant_cells) == 0){
    new_seeds = matrix(0, nrow = rows, ncol = cols, byrow = TRUE) 
  } else {
    propegules = as.vector(plant_cells) %*% seed_dist# this calculates the number of seeds going to each cell, regardless of current occupents
    
    max_propegules = max(propegules)
    chance = propegules/max_propegules
    test_vector = runif(length(propegules))
    new_seeds = 1 * (chance > test_vector)
    
    #finally, new seeds are added to the plant grid - they remain ageless until the next time reproduction occurs though.   
    #new_seeds = vector(mode = "integer", length = (rows * cols))
    #new_seeds[which(propegules >  0)] = 1               # adds a seed to the seed grid when there is a propegule calculated as being present there.
    new_seeds = matrix(new_seeds, nrow = rows, ncol = cols, byrow = TRUE)
  }  
  return(new_seeds)
  
}