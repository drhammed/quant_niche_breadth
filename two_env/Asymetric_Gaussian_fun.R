source("AsymetricGaussian_functions.R")
source("MyNewNicheBreadthFunctions.R")
source("DavidZelenyThetaFunctions.R")
#source("asymetric_gaus_multi_run.R")



#########################################################################################
#########################################################################################
########### demonstration for multiple predictors for multiple species
#########################################################################################


# Number of simulations
#num_simulations <- 100
# Create a list to store the results of each simulation
#simulated_data <- list()
# Run the simulations
#for (i in 1:num_simulations) {
#  simu_data <- generate_community.Asymetric.Gaussian.multi.x(n.species = 30, n.sites = 500, n.x = 2, shape = 1, rate = 1)
#  simulated_data[[i]] <- simu_data
#}

#sim_1 = simulated_data[[1]]

#simulated.com <- generate_community.Asymetric.Gaussian.multi.x_run(n.species=30,n.sites=500,n.x=2,
#                                                               shape = 0.1,rate=0.5,num_runs = 100)

######## lots of rare species (i.e., narrows niches), shape = 0.1, rate=0.5
n.species <- 30
set.seed(99)
simulated.com <- generate_community.Asymetric.Gaussian.multi.x(n.species=30,n.sites=500,n.x=2,
                                                               shape = 1,rate=1)

hist(simulated.com$breadth.Oracle)
hist(simulated.com$left_breadth + simulated.com$right_breadth)
sim.breadth <- simulated.com$left_breadth + simulated.com$right_breadth


min(simulated.com$breadth.Oracle)

PA <- simulated.com$community_abundance_matrix
PA[PA > 0] <- 1 
hist(apply(PA,2,sum))


nicheOVER.res=run.nicheROVER(PA, simulated.com$x)
additive.theta <- calculate.theta(as.data.frame(PA), method = 'add')
nicheBreadth_Gam <- estimate_nicheBreadth_Gam(PA, simulated.com$x)
nicheBreadth_Gam.latents <- estimate_nicheBreadth_Latents(PA, simulated.com$x,nlv=5)
nicheBreadth_avg.Dist <- estimate_nicheBreadth_avg.Dist(PA, simulated.com$x)

# niche.breadth.simulated2 <- simulated.com$left_breadth + simulated.com$right_breadth
cor(cbind(simulated.com$breadth.Oracle,nicheOVER.res,additive.theta$theta,
          nicheBreadth_Gam$Niche_Breadth,nicheBreadth_Gam.latents,
          nicheBreadth_avg.Dist))


######## lots of common species (i.e., broad niches), shape = 5,rate=0.5
simulated.com <- generate_community.Asymetric.Gaussian.multi.x(n.species=30,n.sites=500,n.x=4,
                                                               shape = 5,rate=0.5)

sim.breadth <- simulated.com$left_breadth + simulated.com$right_breadth
PA <- simulated.com$community_abundance_matrix
PA[PA > 0] <- 1 
hist(apply(PA,2,sum))
