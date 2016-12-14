veggrowth_mp = function(plant_cells, dist_cells, vegchance){
  rows = dim(plant_cells)[1]
  cols = dim(plant_cells)[2]
  
  rowscols = rows * cols
  veg_cells = matrix((vegchance), nrow = rowscols, ncol = rowscols) # assigns all cells in this dummy matrix to a quarter chance 
  veg_cells[which(dist_cells > 1.5)] = 0 # plants placed in cells further than 1.5 cells are reduced to  a zero chance of growth
  
  propegules = as.vector((plant_cells)) %*% veg_cells # calculates the "quarters" of vegetative growth for each cell, based on where plants are currently located. 
  new_veg = matrix(0, nrow = rows, ncol = cols) # creates a matrix for storing places where there will be new vegetative growth
  
  # each value in new_veg is the chance of there being veg growth
  # plants are added to new_plants if they pass the test against a randomly generated number.
  random_matrix = runif(rowscols) # assigns a random number between zero and 1 
  
  for (i in 1:length(propegules)){ # compare random matrix to propegules
    new_veg[i] = ifelse(propegules[i] > random_matrix[i], 1, 0)
  } # ifelse(test, yes, no)
  
  return(new_veg = new_veg)}
