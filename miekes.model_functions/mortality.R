mortality = function(ages, chance_factor = 1/10, nspecies){
  death = array(0, dim = dim(ages))
  for (n in 1:nspecies){
    mortality_chance = ages[,,n] * chance_factor
    test = matrix(runif(length(ages[,,n])), nrow = dim(ages)[1], ncol = dim(ages)[2])
    death[,,n] = mortality_chance > test * 1
  }
  return(death = death)}
