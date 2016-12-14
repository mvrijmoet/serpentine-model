initialize_grids = function(rows , cols, positionx, positiony, obspos, obsintens = 0.5){ 
  numspecies <- length(positionx)
  
  plant_grid = array(data = 0, dim = c(rows, cols, numspecies))
  plant_age = array(data = 0, dim = c(rows, cols, numspecies))
  
  obs = matrix(data = 1, nrow = m, ncol = n) 

for (items in 1:length(obspos)){
  obs[obspos[[items]][1,], obspos[[items]][2,]] = obsintens # adjust values of obspos to reflect coordinates from the obspos_coords list
}
  
    new_seeds = array(data = 0, dim = c(rows, cols, numspecies))
  for (i in 1:numspecies){
    plant_grid[positionx[i], positiony[i], i] = 1
    plant_age[positionx[i], positiony[i], i] = 1
  }
  
  return(list(plant_grid = plant_grid, 
              plant_age = plant_age, 
              new_seeds = new_seeds, 
              obs = obs, 
              numspecies = numspecies))
}
