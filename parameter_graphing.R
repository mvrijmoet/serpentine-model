# this script will parameterize and run the model multiple times. 
# outputs: population data, shannon index, winner, for each model run. 

# next:
# fix shannon index - currently bugs when it hits zero (I think) (4-18)
# outputs: text file containing population matrix, paramter values, and final outputs for  each  script run. 

# REMEMBER: setwd to miekes.model_additional before  running this script!

library('lhs') # for Latin Hypercube Sampling

# universal parameters 
m = 24 # number of rows in plant_grid
n = 24 # number of columns in plant_grid
stopifnot(m%%8 == 0)

mn = m * n
generations = 3 * 10 
nspecies = 6
color = rainbow(nspecies, start = 0, end = max(1, nspecies - 1)/nspecies)
runs = 64 * 3
reps = 3

xpos = c()
ypos = c()

# establish empty vectors for parameters
{vegchance = c()
slope = c()
b = c()
tolerence = c()
comp = c()
vegchance_raw = c()
slope_raw = c()
b_raw = c()
tolerence_raw = c()
comp_raw = c()
params = list()
}


# parameters are chosen - lets choose 10 at a time for one plant, and one each for the other two.

# establish a dataframe to store outputs/parameter values
{ plant_bin <- array(data = 0, dim = c(m, n, nspecies, reps * runs))
  population_bin <- array(data = 0, dim = c(generations, nspecies, reps * runs))
  community_info = c("dominent", "shannon", "spdens_toxic", "spdens_norm", "sprich_toxic", 
                     "sprich_norm", "generations", "nspecies", "toxic_area", "total_area", 
                     "frac_dim1", "frac_dim2", "frac_dim3")
  param_names = c("b", "slope", "vegchance", "pop", "tol", "comp")
  param_col = c()
  param_raw_col = c()
  for (sp in 1:nspecies){ # generate the list of parameter columns for each plant
    param_raw_col = c(param_raw_col, paste("raw", param_names, sp, sep = "_"))
    param_col = c(param_col, paste(param_names, sp, sep = "_"))
    
  }
  param_col = c(community_info, param_col, param_raw_col)
  
# all randomized parameters chosen off a latin hypercube
  param_raw = randomLHS(runs, length(param_names) * nspecies) # generates values for all parameters and species, though some of these are dependent within the model. 
  # fractal dimensions generated
  dim1 = rep(c(rep(.25, 16), rep(.5, 16), rep(.75, 16), rep(1.0, 16)), 3)#length = runs 
  dim2 = rep(rep(c(rep(.25, 4), rep(.5, 4), rep(.75, 4), rep(1.0, 4)), 4),3)
  dim3 = rep(c(rep(c(.25, .5, .75, 1.0), 16)),3)
    
  run_stats = data.frame(matrix(0, nrow = reps * runs, ncol = length(param_col)))
  colnames(run_stats) <- param_col 
}
  

# model is run for every set of parameters. 

for (g in 1:runs){

obspos_coords = obspos_generator(m = m, n = n, dim1 = dim1[g], dim2 = dim2[g], dim3 = dim3[g])

xpos[1:nspecies] = round(runif(nspecies, 1, m))
ypos[1:nspecies] = round(runif(nspecies,1,n))

vegchance[1:nspecies] =  param_raw[g, 1:nspecies] #1
vegchance_raw[1:nspecies] = param_raw[g, 1:nspecies]
slope[1:nspecies] =  param_raw[g, 2 * (1:nspecies)] * -1
slope_raw[1:nspecies] = param_raw[g, 2 * (1:nspecies)]
b[1:nspecies] =  param_raw[g, 3 * (1:nspecies)] * 10
b_raw[1:nspecies] = param_raw[g, 3 * (1:nspecies)]
tolerence[1:nspecies] = param_raw[g, 4 * (1:nspecies)]
tolerence_raw[1:nspecies] = param_raw[g, 4 * (1:nspecies)]
comp[1:nspecies] = param_raw[g, 5 * (1:nspecies)]
comp_raw[1:nspecies] = param_raw[g, 5 * (1:nspecies)]

# stores all these parameters in a list, coded by names (for use later, but most of script relys on raw values above)
params[[g]] = list(b, slope, vegchance, tolerence, comp, xpos, ypos, b_raw, slope_raw, vegchance_raw, tolerence_raw, comp_raw)  
names(params[[g]]) <- c("b", "slope", "vegchance", "tol", "comp", "xpos", "ypos", "raw_b", "raw_slope", "raw_vegchance", "raw_tol", "raw_comp")

# the serpentine script is run using these parameters (generic 1 and 2, particular 3)
for (rp in 1:reps){
  sequ = (g-1) * length(1:reps) + rp # names the sequence this falls in in all runs.
  
  source("Serpentine_model.R")
  
  plant_bin[,,,sequ] <- plant_grid[,,] # store plant grids generated

  # plot
  alltogether = matrix(data = 0, nrow = m, ncol = n)
  for (s in 1:nspecies){
    alltogether = alltogether + s * plant_grid[,,s]
  }
  image(t(alltogether), col = color)
  
  population_bin[1:generations,,sequ] <- population
  params[[g]][["pop"]] <- population[generations,]
  
  run_stats$dominent[sequ] = which.max(population_bin[generations,,sequ])
  
  #calculate the shannon index
  p = c(rep(0,nspecies))
  vals = c(rep(0, nspecies))
  for (h in 1:nspecies){
    p[h] = population_bin[generations,h,sequ]/sum(population_bin[generations,,sequ]) #computes relative proportion of species h in final generation for model run f
    if (is.na(p[h])){
      p[h] = 0
    }
    vals[h] = p[h] * log(p[h])
  } #compute relative proportion of each species
  run_stats$shannon[sequ] = -sum(vals) #compute shannon value for assemblage
  if (is.nan(run_stats$shannon[sequ])){
    run_stats$shannon[sequ] = 0
  }
  
  # compute species richness in toxic patches
  zone_toxic = dim1[g] * dim2[g] * dim3[g] # calculate area of toxicity
  plant_bin_tox = plant_bin[,,,sequ]
  for (la in 1:nspecies){
    plant_bin_tox[,,la] = plant_bin[,,la,sequ] * (obs < 1)
  }
  population_toxic = apply(plant_bin_tox, MARGIN = c(3), sum) #calculates population in toxic zones for each plant species
  richness_toxic = length(which(population_toxic > 0))
  run_stats$spdens_toxic[sequ] = richness_toxic / (zone_toxic * mn) 
  run_stats$sprich_toxic[sequ] = richness_toxic
  
  # compute species richness in non-toxic patches
  population_norm = apply(plant_bin[,,,sequ], MARGIN = 3, sum) - population_toxic #calculates population in normal soils for each plant
  richness_norm = length(which(population_norm > 0))
  run_stats$spdens_norm[sequ] = richness_norm / ((mn) - zone_toxic * mn) 
  run_stats$sprich_norm[sequ] = richness_norm
  
  # store additional run stats:
  run_stats$generations[sequ] = generations
  run_stats$nspecies[sequ] = nspecies
  run_stats$toxic_area[sequ] = zone_toxic * mn
  run_stats$total_area[sequ] = mn
  run_stats$frac_dim1[sequ] = dim1[g]
  run_stats$frac_dim2[sequ] = dim2[g]
  run_stats$frac_dim3[sequ] = dim3[g]  
  # store parameter values used
  
  {for (sp in 1:nspecies){
    sp_params <- paste(param_names, sp, sep = "_") # generates secific species by species name for columns in run_stats
    sp_params_raw <- paste("raw", param_names, sp, sep = "_") # generate names for the raw values of these params (randomized proportions from latin hypercube)
    
    for (pp in 1:length(sp_params)){
      this_param = param_names[pp]
      this_param_raw = c(paste("raw", param_names, sep = "_"))[pp]
      run_stats[sequ,sp_params[pp]]  <- params[[g]][[this_param]][sp] 
    if (this_param_raw != "raw_pop"){
      run_stats[sequ,sp_params_raw[pp]] <- params[[g]][[this_param_raw]][sp]
    }
      # pulls this species value for this parameter from the list of vectors of parameters for this run. 
      # stores them by specific name in run_stats
    }}}
  
  # store corresponding raw parameter value for each run for each parameter (for use in multiple regression)
  
  }
  
}  

tally(run_stats$dominent)
View(run_stats)
# store: the population matrix as a layer of an array (population_bin)
# store: the parameters for that run and some calculated statistics - in a data-frame (run_stats) 

