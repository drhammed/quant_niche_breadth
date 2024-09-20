
generate_abundances.Symmetric.Gaussian.multi.x <- function(x, optima, breadth, h) {
  n.sites <- nrow(x)
  n.x <- ncol(x)
  
  
  Abundance.oracle <- matrix(0, n.sites, n.x)
  
  niche.breadth.Oracle <- numeric(n.x)
  
  for (i in 1:n.x) {
    Abundance.oracle[, i] <- rpois(n.sites, 30 * h[i] * sqrt(exp(-(x[, i] - optima[i])^2 / (2 * breadth[i]^2))))
    
    # Calculate the weighted variance
    w.mean <- sum(Abundance.oracle[, i] * x[, i]) / sum(Abundance.oracle[, i])
    
    # Check if the denominator is zero
    if (sum(Abundance.oracle[, i]) == 0) {
      niche.breadth.Oracle[i] <- 0
    } else {
      w.ss <- sum(Abundance.oracle[, i] * (x[, i] - w.mean)^2)
      w.var <- w.ss / sum(Abundance.oracle[, i])
      
      niche.breadth.Oracle[i] <- sqrt(w.var)
    }
  }
  
  
  
  Abundance.oracle <- apply(Abundance.oracle, 1, sum) / n.x
  Abundance.sim <- rpois(length(Abundance.oracle), Abundance.oracle)
  
  output <- list(Abundance.oracle = Abundance.oracle,
                 Abundance.sim = Abundance.sim,
                 niche.breadth.Oracle = niche.breadth.Oracle)
  
  return(output)
}


generate_community.Symmetric.Gaussian.multi.x <- function(n.species, n.sites, n.x) {
  community_abundance_matrix <- matrix(0, n.sites, n.species)
  Abundance.oracle <- matrix(0, n.sites, n.species)
  x <- NULL  # Initialize x
  
  repeat {
    #x <- NULL  # Reset x to NULL at the beginning of each iteration
    x <- matrix(rnorm(n.x * n.sites), n.sites, n.x) 
    optima <- matrix(runif(n.x * n.species, -2, 2), n.x, n.species) 
    h <- matrix(runif(n.x*n.species, min=0.3, max=1), n.x, n.species)
    breadth <- matrix(runif(n.x*n.species, 0.1, 1.2), n.x, n.species)  # Uniform distribution for breadth
    
    breadth.Oracle <- numeric(n.species)
    for (species in 1:n.species) {
      sim <- generate_abundances.Symmetric.Gaussian.multi.x(x, optima[, species], breadth[,species], h[, species])
      community_abundance_matrix[,species] <- sim$Abundance.sim
      breadth.Oracle[species] <- mean(sim$niche.breadth.Oracle)
      Abundance.oracle[,species] <- sim$Abundance.oracle
    }
    n_species_c <- sum(colSums(community_abundance_matrix)!=0) # _c for check
    n_communities_c <- sum(rowSums(community_abundance_matrix)!=0)
    if ((n_species_c == n.species) & (n_communities_c==n.sites)){break}
    
  }
  colnames(community_abundance_matrix) <- paste0("sp", 1:n.species)
  names(breadth.Oracle) <- paste0("sp", 1:n.species)
  
  output <- list(community_abundance_matrix=community_abundance_matrix,
                 breadth.Oracle=breadth.Oracle, x=x,
                 Abundance.oracle=Abundance.oracle)
  return(output)  
}
