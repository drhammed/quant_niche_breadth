#########################################################################################
#################  Asymmetric Gaussian to simulate species distribution ################## 
############## methods and functions developed by Pedro Peres-Neto, May 2023 ############
#########################################################################################

########### utility functions

set_max_min <- function(x, new_max, new_min) {
  # Find the current maximum and minimum values of x
  current_max <- max(x)
  current_min <- min(x)
  
  # Calculate the range of x
  current_range <- current_max - current_min
  
  # Calculate the range of the new values
  new_range <- new_max - new_min
  
  # Scale x to the new range
  x_scaled <- (x - current_min) * new_range / current_range + new_min
  
  # Replace the maximum and minimum values with the desired values
  x_scaled[x == current_max] <- new_max
  x_scaled[x == current_min] <- new_min
  
  return(x_scaled)
}

weighted.variance <- function(x, w) {
  # Check if x is a matrix or data.frame
  if (!is.matrix(x) && !is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  # Initialize an empty vector to store the weighted variances
  w.var <- numeric(ncol(x))
  
  # Loop through each column of x
  for (i in 1:ncol(x)) {
    # Calculate the weighted mean
    w.mean <- sum(w * x[, i]) / sum(w)
    
    # Calculate the weighted sum of squares
    w.ss <- sum(w * (x[, i] - w.mean)^2)
    
    # Calculate the weighted variance and store it in w.var vector
    w.var[i] <- w.ss / sum(w)
  }
  return(w.var)
}


########### function to generate abundances for a single species and for a single environmental
# gradient (i.e., x is a vector)
generate_abundances.Asymetric.Gaussian <- function(x, optimum, peak, left_breadth, right_breadth){
  # simulate a single species abundances as a function of x using asymmetric Gaussian function
  left_side <- ifelse(x < optimum, peak * exp(-((x - optimum) / left_breadth)^2), 0)
  right_side <- ifelse(x >= optimum, peak * exp(-((x - optimum) / right_breadth)^2), 0)
  Abundance.oracle <- left_side + right_side
  Abundance.sim <- rpois(length(Abundance.oracle),lambda=Abundance.oracle) 
  output <- list(Abundance.oracle=Abundance.oracle,Abundance.sim=Abundance.sim)
  return(output)
  # NOTE - function generate_abundances.Asymmetric.Gaussian.multi.x below
  # can also be used for a single predictor but easy to understand 
  # how a non-symmetric Gaussian works based on the implementation here.
}

########### functions that generalizes the above for multiple species

# function to generate abundances for a single species but multiple predictors
# the function simulates abundances as in a Poisson regression under an asymmetric 
# Gaussian. Variable strength (akin to slopes for each environment; x column) is
# somewhat not as easy as the Poisson model is set a linear. This is further
# discussed below and the solution used works well under certain constraints (also discussed below)

generate_abundances.Asymetric.Gaussian.multi.x <- function(x, optima, peak, left_breadth, right_breadth,strength.x){
  # simulate single species abundances as a function of multiple x (a single x can be used too) using an asymmetric Gaussian function
  # it generates for a single species
  n.sites <- NROW(x)
  n.x <- NCOL(x)
  
  left_matrix <- matrix(0,n.sites,n.x)
  right_matrix <- matrix(0,n.sites,n.x)
  Abundance.oracle <- matrix(0,n.sites,n.x)
  
  niche.breadth.Oracle <- 0
  
  for(i in 1:n.x){
    # the two lines of code below are based on the fact that 50 * exp(10) = exp(log(50) + 10) &
    # that 50 * exp(-10) = exp(log(50) - 10)
    # note also that the reciprocal (1/strength.x[i]) make it more intuitive for determining
    # predictor strength of non-linear poisson
    left_matrix[,i]  <- ifelse(x[,i] <  optima[i], exp((log(peak) - (1/strength.x[i]) * ((x[,i] - optima[i]) / left_breadth[i])^2)), 0)
    right_matrix[,i] <- ifelse(x[,i] >= optima[i], exp((log(peak) - (1/strength.x[i]) * ((x[,i] - optima[i]) / right_breadth[i])^2)), 0)
    #left_breadth.Oracle <- left_breadth.Oracle + sqrt(weighted.variance(x[,i], left_matrix[,i]))
    #right_breadth.Oracle <- right_breadth.Oracle + sqrt(weighted.variance(x[,i], right_matrix[,i]))
    Abundance.oracle[,i] <- left_matrix[,i] + right_matrix[,i]
    niche.breadth.Oracle <- niche.breadth.Oracle + sqrt(weighted.variance(x[,i], Abundance.oracle[,i]))
    
    #niche.breadth.Oracle <- niche.breadth.Oracle + 
    #                  (sqrt(weighted.variance(x[,i], left_matrix[,i])) + sqrt(weighted.variance(x[,i], right_matrix[,i])))
  }
  
  Abundance.oracle <- apply(Abundance.oracle,1,sum) / n.x
  Abundance.sim <- rpois(length(Abundance.oracle),Abundance.oracle)
  #left_breadth.Oracle <- left_breadth.Oracle/n.x
  #right_breadth.Oracle <- right_breadth.Oracle/n.x
  breadth.Oracle <- niche.breadth.Oracle/n.x
  
  output <- list(Abundance.oracle=Abundance.oracle,
                 Abundance.sim=Abundance.sim,
                 niche.breadth.Oracle=niche.breadth.Oracle)
  return(output)
}

generate_community.Asymetric.Gaussian.multi.x <- function(n.species,n.sites,n.x,shape,rate){
  # generates multiple species abundances for multiple environmental variables
  # number of sites and number of communities mean the same here
  
  community_abundance_matrix <- matrix(0,n.sites,n.species)
  Abundance.oracle <- matrix(0,n.sites,n.species)
  repeat {
    x <- matrix(rnorm(n.x * n.sites),n.sites,n.x) 
    optima <- matrix(runif(n.x * n.species,-2,2),n.x,n.species) 
    peak <- round(runif(n.species,10,500))
    # strength.x <- matrix(runif(n.x * n.species,1),n.x,n.species)
    strength.x <- matrix(runif(n.x * n.species,1,1),n.x,n.species) # fixed every variable for every species at 1
    
    # Uniform distribution
    left_breadth <- matrix(runif(n.x * n.species,0.1,1.2),n.x,n.species)  
    right_breadth <- matrix(runif(n.x * n.species,0.1,1.2),n.x,n.species)
    
    
    # makes all x variables have the same breadth for any given species
    #left_breadth <- asymmetric_data <- rgamma(n.species, shape = shape, rate = rate)
    #left_breadth <- set_max_min(left_breadth, new_max=5, new_min=0.05)
    #left_breadth <- t(matrix(rep(left_breadth,n.x),n.species,n.x))
    #right_breadth <- asymmetric_data <- rgamma(n.species, shape = shape, rate = rate)
    #right_breadth <- set_max_min(right_breadth, new_max=5, new_min=0.05)
    #right_breadth <- t(matrix(rep(right_breadth,n.x),n.species,n.x))
    
    breadth.Oracle <- numeric(n.species)
    for(species in 1:n.species){
      sim <- generate_abundances.Asymetric.Gaussian.multi.x(x, optima[,species], peak[species], left_breadth[,species], right_breadth[,species],strength.x[,species])
      community_abundance_matrix[,species] <- sim$Abundance.sim
      breadth.Oracle[species] <- sim$niche.breadth.Oracle
      Abundance.oracle[,species] <- sim$Abundance.oracle
    }
    n_species_c <- sum(colSums(community_abundance_matrix)!=0) # _c for check
    n_communities_c <- sum(rowSums(community_abundance_matrix)!=0)
    if ((n_species_c == n.species) & (n_communities_c==n.sites)){break}
  }
  colnames(community_abundance_matrix) <- paste0("sp", 1:n.species)
  # a given species has the same niche breadth for all x variables as set above; so, just pick the 1st x 
  # both left_breadth/right_breadth gives the same results as left_breadth.from.oracle/right_breadth.from.oracle; 
  # just a way to demonstrate that generating data using the asymmetric Gaussian was properly coded
  left_breadth <- left_breadth[1,]
  right_breadth <- right_breadth[1,] 
  names(left_breadth) <- paste0("sp", 1:n.species)
  names(right_breadth) <- paste0("sp", 1:n.species)
  
  #names(left_breadth.from.oracle) <- paste0("sp", 1:n.species)
  #names(right_breadth.from.oracle) <- paste0("sp", 1:n.species)
  names(breadth.Oracle) <- paste0("sp", 1:n.species)
  #breadth.nicheOVER <- run.nicheOVER(Abundance.oracle, x)
  breadth.nicheOVER <- c()
  output <- list(community_abundance_matrix=community_abundance_matrix,left_breadth=left_breadth,
                 right_breadth=right_breadth,breadth.Oracle=breadth.Oracle,x=x,
                 breadth.nicheOVER=breadth.nicheOVER,Abundance.oracle=Abundance.oracle)
  return(output)  
}



