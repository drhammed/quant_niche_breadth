#########################################################################################
########################## All methods to estimate niche breadth- Hammed Akande ########################
############## this include the new ones by Pedro Peres-Neto ############
#########################################################################################


#index 1- co-occurence based estimates

co_occur <- function(GOSmat, reps, psample, psample2) {
  
  #create function to compute similarity index (sim)
  
  "sim" <- 
    function(x, coord=NULL, method="soer", dn=NULL, normalize = FALSE, listin = FALSE, listout = FALSE, ...)
    {	
      if (!is.na(pmatch(method, "jaccard"))) 
        method <- "jaccard"
      METHODS <- c("soerensen", "jaccard", "simpson")
      method <- pmatch(method, METHODS)
      if (is.na(method)){
        stop("invalid similarity method")
      }
      if (method == -1){
        stop("ambiguous similarity method")
      }
      if (listin) {
        x <- mama(x)
        x <- as.matrix(x)
      }
      x <- x > 0
      df <- as.matrix(x)
      zeina <- row.names(df)
      anz <- nrow(df)
      a <- df %*% t(df) 
      b <- df %*% (1 - t(df)) 
      c <- (1 - df) %*% t(df) 
      d <- ncol(df) - a - b - c
      if (normalize) {
        an <- a/(a+b+c)
        bn <- b/(a+b+c)
        cn <- c/(a+b+c)
        a <- an
        b <- bn
        c <- cn
      }
      if (method == 1) {
        dis <- (2*a)/((2*a) + b + c)
      }
      else if (method == 2) {
        dis <- a / (a + b + c)
      }
      else if (method == 3) {
        dis <- pmin(b,c) / (pmin(b,c) + a) 
      }
      
      
      dis <- as.dist(dis)
      attr(dis, "Size") <- anz
      attr(dis, "Labels") <- zeina
      attr(dis, "method") <- METHODS[method]
      attr(dis, "call") <- match.call()
      class(dis) <- "dist"
      if (listout) {
        dis <- liste(dis, entry=METHODS[method])
        dis$a <- a[row(a) > col(a)]
        dis$b <- b[row(b) > col(b)]
        dis$c <- c[row(c) > col(c)]
        dis$d <- d[row(d) > col(d)]
      }
      if (!is.null(coord)){
        xydist <- liste(dist(coord), entry="distance")
        dis <- cbind(xydist, as.vector(dis))
        names(dis)[4] <- METHODS[method]
        X <- (outer(coord[,1], coord[,1], FUN="+"))*0.5
        Y <- (outer(coord[,2], coord[,2], FUN="+"))*0.5	   
        dis$X <- X[row(X) > col(X)]
        dis$Y <- Y[row(Y) > col(Y)]
        dis$xdist <- dist(coord[,1])
        dis$ydist <- dist(coord[,2])
        dis$a <- a[row(a) > col(a)]
        dis$b <- b[row(b) > col(b)]
        dis$c <- c[row(c) > col(c)]
        dis$d <- d[row(d) > col(d)]
        if (!is.null(dn)) {
          if(length(dn)==1){
            dis <- dis[(dis$distance <= dn), ]
          }
          else{
            dis <- dis[((dis$distance >= min(dn)) & (dis$distance <= max(dn))), ]
          }
        }
      }
      return(dis)
    }
  
  
  
  
  # Create function for multi-site simpson based on the formula presented in Balsega 2007 
  #(also, add multi-site Soerensen index)
  
  
  minbibj<-function(matrix) {
    ## nr is the number of rows
    nr<-dim(matrix)[1];
    ## This variable contains the sum of 'min(bi,bj)' on row pairs.
    sumbibj<-0.;
    ## We 'loop' on every different row pairs.
    for(i in 1:nr-1) {
      for(j in (i+1):nr) { 
        ## 'bi' and 'bj' are respectively the number of species appearing
        ## only in the ith and the jth site (row).
        
        bi<-sum(matrix[i,]&(!matrix[j,]));
        bj<-sum(matrix[j,]&(!matrix[i,]));
        
        ## We sum up 'min(bi,bj)' on every different pair of rows.
        sumbibj<-sumbibj+min(bi,bj);
      }
    }
    
    ## We return the sum
    sumbibj;
  }
  zeroColumn<-function(matrix) {
    sum<-0;
    nc<-dim(matrix)[2];
    for(i in 1:nc) if(!sum(matrix[,i])) sum<-sum+1;
    #if(sum!=0)
    #  warning(sum," ",call.=FALSE,immediate.=TRUE);
    sum
  }
  
  Simpson.multi<-function(x) {
    ## x must be a data frame
    
    matrix<-as.matrix(x);
    
    ## Si is the number of species present in the ith site. We sum all
    ## these values.
    
    sumSi<-sum(matrix);
    
    ## St = total number of species in all sites, thus the number of columns, excepting if a given
    ## species is not present at all (in any site), so we must not count this column
    ## out.
    
    St<-ncol(matrix)-zeroColumn(matrix);
    
    ## 'a' is the number of species shared for at least two sites
    
    a<-sumSi-St;
    index<-a/(minbibj(matrix)+a);
    result<-index
  }
  
  Soerensen.multi<-function(x) {
    
    ## x must be a data frame
    
    matrix<-as.matrix(x);
    
    ## Si is the number of species present in the ith site. We sum all
    ## these values.
    
    sumSi<-sum(matrix);
    
    ## St is the total number of species in all sites, thus the number of columns, excepting if a given
    ## species is not present at all (in any site), so we must not count this column
    ## out.
    
    St<-ncol(matrix)-zeroColumn(matrix);
    
    ## T is the number of sites in the matrix
    
    T<-nrow(matrix)
    
    index<-(T/(T-1))*(1-(St/sumSi))
    result<-index
  }    
  
  
  
  #convert the species presence and absence to a 2 column list for the generalist-specialist algorithm
  sim.com
  spp.vec<-NULL
  plot.vec<-NULL
  n_sites_sel <- nrow(sim.com)
  for(i in 1:n_sites_sel) {
    vec.true<-as.logical(sim.com[i,])
    #print(sim.com[i,])
    #print(sum(sim.com[i,]))
    plot.vec<-c(plot.vec,rep(i,length=sum(sim.com[i,])))
    spp.vec<-c(spp.vec,c(1:n.species)[vec.true])
  }
  
  #output
  out.simul <- data.frame(plot.vec,spp.vec)
  
  
  # set up loop format for generalist-specialist calculation
  
  c.range<-list(NULL)
  sp.loop <- 1
  GOSmat<- out.simul
  
  #spp matrix- col 1 is numeric SppID, col 2 is species label (name, "sp" in this case)
  SppMat<-data.frame(sort(unique(GOSmat[,2])),paste("sp",sort(unique(GOSmat[,2])),sep="")) 
  
  #factorized plot vector
  plotID<-factor(GOSmat[,1])		
  
  #species per plot as numbers or codes
  SppID<-GOSmat[,2]			
  
  #number of plots
  Nplots<-length(levels(plotID))		
  #vector of number of species in each plot
  richness<-tapply(SppID,plotID,length)	
  max.rich<-max(richness)			#maximum local richness value
  metacom<-table(plotID,SppID)		#plot x species matrix
  
  #Select subset of species for analysis that occur in "plot.cut" # of plots or more
  plots.per.spp<-tapply(plotID,SppID,length)			#vector of number plot occurrences for each species
  Species<-sort(unique(GOSmat[,2]))[plots.per.spp>=psample]	#vector of selected species
  Nspp<-length(Species)						#number of selected species
  
  
  # loop through all species and calculate each metric
  
  #SPECIES LOOP
  sci.name<-rep(0,Nspp)
  Beta.a<-rep(0,Nspp)
  Beta.w<-rep(0,Nspp)
  multi.sim<-rep(0,Nspp)
  
  #set an empty data frame for the specific stat you want to output
  #you can add more or reduce here
  Theta.out <- data.frame(
    sci.name = character(Nspp),
    multi.sim = numeric(Nspp),
    Beta.a = numeric(Nspp),
    Beta.w = numeric(Nspp)
  )
  
  for(sp in 1:Nspp) {
    
    #print(sp)
    
    #Plot selection
    lab<-as.numeric(labels(metacom)[2][[1]])
    xlab<-c(1:dim(metacom)[2])
    metacol<-xlab[lab==Species[sp]]
    sp.plots<-as.logical(metacom[,metacol])
    sp.metacom<-metacom[sp.plots,]
    Np<-dim(sp.metacom)[1]
    wide<-length(xlab)
    
    
    ##Loop to calculate different similarity measures of Beta diversity from each of "reps" random selections of "psample2" plots
    
    multi.sim.rep<-rep(0,reps)
    
    
    for(k in 1:reps) {
      
      ###selects psample2 plots randomly from sp.metacom
      rowselect<-sample(row.names(sp.metacom),psample2)    #selects random plot-IDs from sp.metacom
      trueRow<-is.element(row.names(sp.metacom),rowselect) #logic vector of true and false rows
      rand.mat<-sp.metacom[trueRow,]                       #gives the smaller subset of psample2 plots
      
      
      #Calculate dissimilarity measures
      
      multi.sim.rep[k]<-1-(Simpson.multi(rand.mat))
      
    }
    
    #calculate the mean for simpson
    multi.sim[sp]<-mean(multi.sim.rep)
  
    #Monte Carlo procedure
    
    rpmat<-matrix(c(1:Np),reps,Np,byrow=T)					#"reps" rows of plot sequences
    rpmat<-t(apply(rpmat,1,function(x)sample(x,psample))	)		#randomize plot sequence orders, taking "psample" plots
    mc.mat<-array(0,dim=c(psample,wide,reps))				#monte carlo matrix: psamples x #allspecies x #reps
    for(i in 1:reps) {
      mc.mat[,,i]<-sp.metacom[rpmat[i,],]
    }
    
    
    #-----------
    colsum<-apply(mc.mat,c(2,3),sum)		#sum columns of each rep, reps become columns
    colsum[colsum>0]<-1				#convert >0 values to ones
    rich.vec<-colSums(colsum)-1			#vector of # cooccurrences for each rep
    rich.vec2<-colSums(colsum)      #vector of # species totals for each rep
    mc.mat[mc.mat>0]<-1				#convert species numbers to ones
    rowsum<-apply(mc.mat,c(1,3),sum)		#sum rows of each rep, reps become columns
    Walpha.vec<-colMeans(rowsum)			#vector of "avg local richness" (Whittaker's alpha) for each rep
    Wbeta.vec<-rich.vec-Walpha.vec
    
    
    
    Beta.a.vec<-rich.vec-Walpha.vec        #Landes additive partitioning for each rep
    Beta.w.vec<-rich.vec2/Walpha.vec       #Whittaker's beta calculation for each rep
    Beta.a[sp]<-mean(Beta.a.vec)       #mean additive partitioning value for all reps 
    Beta.w[sp]<-mean(Beta.w.vec)    #mean Whittakers beta value for all reps
    
    
    sci.name[sp]<-as.character(SppMat[,2][SppMat[,1]==Species[sp]])		#scientific name
    
    # Store the results for the current species in the result data frame
    Theta.out$sci.name[sp] <- sci.name[sp]
    Theta.out$multi.sim[sp] <- multi.sim[sp]
    Theta.out$Beta.a[sp] <- Beta.a[sp]
    Theta.out$Beta.w[sp] <- Beta.w[sp]
  }
  
  # Return the result data frame
  return(Theta.out)
  
}  




##niche breadth based on the omi index- index 2

omi_params <- function(env_vars, PA, nf = 2) {
  # Perform PCA on environmental variables
  pca.Env.virt <- ade4::dudi.pca(env_vars, center = TRUE, scale = TRUE, scannf = FALSE, nf = nf)
  
  # Calculate niche parameters using niche() and niche.param()
  omi <- niche(pca.Env.virt, as.data.frame(PA), scannf = FALSE)
  omi.param <- niche.param(omi)
  
  # Extract the tolerance values
  om_tol <- as.data.frame(omi.param)$Tol
  om_df <- as.data.frame(cbind(niche.breadth, om_tol))
  colnames(om_df) <- c("niche_b", "sci.name", "om_tol")
  # Combine the true niche breadth (nb) and omi tolerance values into a data frame
  result_df <- om_df
  
  return(result_df)
}

#for when using replications, use this for omi_fun 

omi_params_fun <- function(env_vars, PA, nf = 2) {
  # Perform PCA on environmental variables
  pca.Env.virt <- ade4::dudi.pca(env_vars, center = TRUE, scale = TRUE, scannf = FALSE, nf = nf)
  
  # Calculate niche parameters using niche() and niche.param()
  omi <- niche(pca.Env.virt, as.data.frame(PA), scannf = FALSE)
  omi.param <- niche.param(omi)
  
  # Extract the tolerance values
  om_tol <- as.data.frame(omi.param)$Tol
  nb <- as.data.frame(niche.breadth)
  nb$sci.name <- row.names(nb)
  om_df <- as.data.frame(cbind(nb, om_tol))
  colnames(om_df) <- c("niche_b", "sci.name" ,"om_tol")
  # Combine the true niche breadth (nb) and omi tolerance values into a data frame
  result_df <- om_df
  
  return(result_df)
}



##additive theta - index 3


add_theta <- function(sim_com, m_sim_df, method = 'add') {
  # Calculate additive theta using calculate.theta()
  additive.theta <- calculate.theta(sim_com, method = method)
  
  # Combine the niche_b and theta values into a data frame
  add_t <- as.data.frame(cbind(m_sim_df$niche_b, additive.theta$theta))
  colnames(add_t) <- c("niche_b", "theta")
  
  return(add_t)
}

#when reps, using this add_theta_fun
add_theta_fun <- function(sim_com, niche_breadth, method = 'add') {
  # Calculate additive theta using calculate.theta()
  additive.theta <- calculate.theta(sim_com, method = method)
  
  # Combine the niche_b and theta values into a data frame
  add_t <- as.data.frame(cbind(niche_breadth, additive.theta$theta))
  colnames(add_t) <- c("niche.breadth", "theta")
  
  return(add_t)
}


#Levins Index - index 4


#####################################
#Levins Index
####################################


Levins <- function(data, matrix = NULL, community = FALSE, abundance = FALSE) {
  
  ##### Species
  # check of input data
  if(nrow(data[!complete.cases(data),]) > 0) stop('NAs present in the data file.') # NAs
  if(missing(matrix)) { # Species
    d.prop <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data))) # proportion
    for (i in 1:nrow(d.prop)) {d.prop[i, ] <- data[i, ] / rowSums(data[i, ])}
    d.sq <- d.prop^2 # square
    d.ba <- data.frame(Levins = matrix(nrow = nrow(d.sq))) # Levins
    for (i in 1:nrow(d.ba)) {d.ba[i, 1] <- 1 - (1 / sum(d.sq[i, ]) - 1) / (ncol(d.sq) - 1)}
    rownames(d.ba) <- rownames(data)
    invisible(d.ba)
  } # Species
  ##### Assemblages
  else { # else Species
    # check of input data
    if(nrow(data[!complete.cases(data),]) > 0) stop('NAs present in the data file.') # NAs
    if(nrow(matrix[!complete.cases(matrix),]) > 0) stop('NAs present in the community matrix.') # NAs
    if(nrow(data) != ncol(matrix)) stop('Number of rows in data file does not match number of columns in the community matrix.')
    spp <- data.frame(data = rownames(data), matrix = colnames(matrix)) ; spp$F = spp[, 1] == spp[, 2] # species order
    4
    if (nrow(spp[spp[, 3] == F, ]) > 0) stop('Order of rownames in data file does not match the order of colnames in the community matrix, or the names differ.')
    matrix.pa <- data.frame(ifelse(matrix == 0, yes = 0, no = 1)) # p/a matrix
    ### (1) Specialization as a mean of present species
    if(community == F) { # 1
      d.prop <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data))) # proportion
      for (i in 1:nrow(d.prop)) {d.prop[i, ] <- data[i, ] / rowSums(data[i, ])}
      d.sq <- d.prop^2 # square
      d.ba <- data.frame(Levins = matrix(nrow = nrow(d.sq))) # Levins
      for (i in 1:nrow(d.ba)) {d.ba[i, 1] <- 1 - (1 / sum(d.sq[i, ]) - 1) / (ncol(d.sq) - 1)}
      pa.sp <- matrix.pa
      for(i in 1:ncol(pa.sp)) {pa.sp[, i] <- ifelse(pa.sp[, i] == 1, yes = d.ba[i, 1], no = NA)}
      # (1a) without abundance
      if (abundance == F) { # 1a
        loc.sp.pa <- data.frame(SR = rowSums(matrix.pa), Levins = rowMeans(pa.sp, na.rm = T)) # Levins
        invisible(loc.sp.pa)
      } # 1a
      # (1b) with abundance
      else { # 1b
        loc.sp.ab <- data.frame(SR = rowSums(matrix.pa), Individuals = rowSums(matrix), Levins = NA) # empty data frame
        for (i in 1:nrow(loc.sp.ab)) {loc.sp.ab[i, 3] <- weighted.mean(pa.sp[i, ], matrix[i, ], na.rm = T)} # Levins
        invisible(loc.sp.ab)
      } #
    } # 
    ### Specialization of the whole assemblage
    else { 
      #  without abundance
      if (abundance == F) { # 2a
        locs.pa <- lapply(1:nrow(matrix.pa), function(x) list()) # empty list
        for (i in 1:nrow(matrix.pa)) {
          locs.pa[[i]] <- data[colnames(matrix.pa[i, matrix.pa[i, ] == 1]), ] # for each assemblage take data of only the present species
          locs.pa[[i]] <- colSums(locs.pa[[i]])} # sum of categories of all species present in an assemblage
        loc <- data.frame(matrix(unlist(locs.pa), nrow = nrow(matrix.pa), byrow = T))
        colnames(loc) <- colnames(data) ; rownames(loc) <- rownames(matrix.pa)
        loc.prop <- loc # proporce
        for (i in 1:nrow(loc.prop)) {loc.prop[i, ] <- loc[i, ] / rowSums(loc[i, ])}
        loc.sq <- loc.prop^2 # square
        loc.com.pa <- data.frame(SR = rowSums(matrix.pa), Levins = NA) # empty data frame
        for (i in 1:nrow(loc.com.pa)) {loc.com.pa[i, 2] <- 1 - (1 / sum(loc.sq[i, ]) - 1) / (ncol(loc.sq) - 1)} # Levins
        invisible(loc.com.pa)
      } 
      # with abundance
      else { 
        locs.ab <- lapply(1:nrow(matrix), function(x) list()) # list
        for (i in 1:nrow(matrix)) {
          locs.ab[[i]] <- data[colnames(matrix[i, matrix[i, ] > 0]), ] # for each assemblage take data of only the present species
          ab <- matrix[i, colnames(matrix[i, matrix[i, ] > 0])] # species abundance for each assemblage
          for (j in 1:nrow(locs.ab[[i]])) {locs.ab[[i]][j, ] <- locs.ab[[i]][j, ] * ab[, j]} # multiply kategories with abundance
          locs.ab[[i]] <- colSums(locs.ab[[i]]) # sum of all species in every category for each locality
        }
        loc <- data.frame(matrix(unlist(locs.ab), nrow = nrow(matrix), byrow = T))
        colnames(loc) <- colnames(data) ; rownames(loc) <- rownames(matrix)
        loc.prop <- loc # proportion
        for (i in 1:nrow(loc.prop)) {loc.prop[i, ] <- loc[i, ] / rowSums(loc[i, ])}
        loc.sq <- loc.prop^2 # square
        loc.com.ab <- data.frame(SR = rowSums(matrix.pa), Individuals = rowSums(matrix), Levins = NA) # prazdny data frame
        for (i in 1:nrow(loc.com.ab)) {loc.com.ab[i, 3] <- 1 - (1 / sum(loc.sq[i, ]) - 1) / (ncol(loc.sq) - 1)} # Levins
        invisible(loc.com.ab)
      } 
    } 
  } # else Species
} # Levins


#melt the data frame


df_melt_fun <- function(data, env_vars) {
  # Combine the sim.com and environmental variables
  df_melt <- cbind(data, env_vars)
  
  # Melt the data frame
  df_melted <- gather(df_melt, key = "species", value = "PA", sp1:sp30)
  
  # Select the species, presence-absence, and environmental variables
  df_melted <- df_melted[, c("species" , "PA", "env.1", "env.2")]
  
  # Delete rows where PA == 0
  df_melted <- df_melted[df_melted$PA != 0, ]
  
  # Remove the PA column
  df_melted <- df_melted[, -2]
  
  return(df_melted)
}




#Bayesian approach (Swandon et al. 2015) - index 5

#Calculate hypervolume niche size estimates

nr_hypervolume <- function(df_melted, nsamples = 1000) {
  # Calculate niche size estimates using niw.post() for each species
  sim.par <- tapply(1:nrow(df_melted), df_melted$species,
                    function(ii) niw.post(nsamples = nsamples, X = df_melted[ii, 2:3]))
  
  # Calculate the posterior distribution of niche size by species
  sim.size <- sapply(sim.par, function(spec) {
    apply(spec$Sigma, 3, niche.size, alpha = 0.95)
  })
  
  # Calculate point estimate and standard error
  rover_nz <- rbind(est = colMeans(sim.size),
                    se = apply(sim.size, 2, sd))
  

  # Extract the hypervolume size estimated from the niche rover
  
  
  nicheRover_hyp <- as.data.frame(t(t(rover_nz[-2, ])))
  colnames(nicheRover_hyp) <- "nr_hypervolume"
  nicheRover_hyp$sci.name <- rownames(nicheRover_hyp)
  
  # Select matching rows between niche breadth and hypervolume estimates
  nb_sim <- nb[nb$sci.name %in% nicheRover_hyp$sci.name, ]
  nb_sim <- nb_sim[, c("sci.name", "niche.breadth")]
  nicheRover_hyp_m <- merge(nb_sim, nicheRover_hyp, by = "sci.name", sort = FALSE)
  colnames(nicheRover_hyp_m) <- c("sci.name", "niche_b", "nr_hypervolume")
  
  return(nicheRover_hyp_m)
}



#in case of reps, using this nr_hypervolume_fun


nr_hypervolume_fun <- function(data, env_vars, nsamples = 1000) {
  # Combine the sim.com and environmental variables
  df_melted <- cbind(data, env_vars)
  
  # Melt the data frame
  df_melted <- gather(df_melted, key = "species", value = "PA", sp1:sp30)
  
  # Select the species, presence-absence, and environmental variables
  df_melted <- df_melted[, c("species" , "PA", "env.1", "env.2")]
  
  # Delete rows where PA == 0
  df_melted <- df_melted[df_melted$PA != 0, ]
  
  # Remove the PA column
  df_melted <- df_melted[, -2]
  
  # Calculate niche size estimates using niw.post() for each species
  sim.par <- tapply(1:nrow(df_melted), df_melted$species,
                    function(ii) niw.post(nsamples = nsamples, X = df_melted[ii, 2:3]))
  
  
  # Calculate the posterior distribution of niche size by species
  sim.size <- sapply(sim.par, function(spec) {
    apply(spec$Sigma, 3, niche.size, alpha = 0.95)
  })
  
  # Calculate point estimate and standard error
  rover_nz <- rbind(est = colMeans(sim.size),
                    se = apply(sim.size, 2, sd))
  
  
  # Extract the hypervolume size estimated from the niche rover
  
  
  nicheRover_hyp <- as.data.frame(t(t(rover_nz[-2, ])))
  colnames(nicheRover_hyp) <- "nr_hypervolume"
  nicheRover_hyp$sci.name <- rownames(nicheRover_hyp)
  
  nb = as.data.frame(niche.breadth)
  nb$sci.name <- rownames(nb)
  
  # Select matching rows between niche breadth and hypervolume estimates
  nb_sim <- nb[nb$sci.name %in% nicheRover_hyp$sci.name, ]
  nb_sim <- nb_sim[, c("sci.name", "niche.breadth")]
  nicheRover_hyp_m <- merge(nb_sim, nicheRover_hyp, by = "sci.name", sort = FALSE)
  colnames(nicheRover_hyp_m) <- c("sci.name", "niche_b", "nr_hypervolume")
  
  return(nicheRover_hyp_m)
}




run.nicheROVER_p <- function(data, env_vars){
  
  data <- as.matrix(data)
  n.env.variables <- NCOL(env_vars)
  n.species <- ncol(data)
  n.sites <- nrow(data)
  
  vectorized <- data.frame(as.factor(rep(colnames(data),each=n.sites)),
                           c(data),env_vars[rep(1:n.sites, times = n.species), ])
  
  eliminate.zeros <- which(vectorized[,2]==0)
  vectorized <- vectorized[-eliminate.zeros,]
  vectorized <- vectorized[,c(-2)]
  colnames(vectorized) <- c("species",paste0("env.", 1:n.env.variables))
  
  nsamples <- 1000
  vectorized.par <- tapply(1:nrow(vectorized), vectorized$species,
                           function(ii) niw.post(nsamples = nsamples, X = vectorized[ii,2:(n.env.variables+1)]))
  vectorized.size <- sapply(vectorized.par, function(spec) {
    apply(spec$Sigma, 3, niche.size, alpha = .95)
  })
  
  desired_order <- colnames(data)
  ordered_vectorized.size <- vectorized.size[, match(desired_order, colnames(vectorized.size))]
  nr_nb <- apply(ordered_vectorized.size,2,mean)
  
  #nb
  nb <- as.data.frame(niche.breadth)
  nb$sci.name <- row.names(nb)
  
  
  nr_nb = as.data.frame(nr_nb)
  nr_nb$sci.name <- rownames(nr_nb)
  nr_nb <- merge(nb, nr_nb, by= "sci.name", sort = FALSE)
  colnames(nr_nb) <- c("sci.name", "niche_b", "nicheOVER_ppn")
  
  return(nr_nb)
}








#Hypervolume (ref. Blonder et al. 2018) - index 6


hypervolume_blond <- function(df_melted) {
  
  #Extract environmental variables for each species and calculate hypervolume
  
  
  # Create an empty list to store the results for each species
  hv_list <- list()
  
  # Get unique species values from the "species" column
  species_values <- unique(df_melted$species)
  
  # Iterate over each species
  for (sp in species_values) {
    # Subset the data for the current species
    sim_data <- df_melted[df_melted$species == sp, c("env.1", "env.2")]
    # Calculate the hypervolume
    hv <- hypervolume_gaussian(sim_data, name = sp, samples.per.point = 10)
    
    # Store the hypervolume object in the list
    hv_list[[sp]] <- hv
  }
  
  # Create an empty data frame to store the results
  hv_results <- data.frame(species = character(), hypervolume = numeric(), stringsAsFactors = FALSE)
  
  # Populate the hv_results data frame using the values stored in the list
  for (sp in species_values) {
    sp_hv <- get_volume(hv_list[[sp]])
    hv_results <- rbind(hv_results, data.frame(species = sp, hypervolume = sp_hv))
  }
  
  nb <- as.data.frame(niche.breadth)
  nb$species <- row.names(nb)
  #cbind the niche breadth and the hypervolume for each species
  hyp_nb_df <- cbind(as.data.frame(nb$niche.breadth), hv_results$hypervolume)
  colnames(hyp_nb_df)<- c("niche_b", "hypervolume")
  
  return(hyp_nb_df)
}



# in case of reps, use this hypervolume_blond_fun

hypervolume_blond_fun <- function(data, env_vars) {
  # Combine the species presence-absence data and environmental variables
  df_melted <- cbind(data, env_vars)
  
  # Melt the data frame
  df_melted <- gather(df_melted, key = "species", value = "PA", sp1:sp30)
  
  # Select the species, presence-absence, and environmental variables
  df_melted <- df_melted[, c("species" , "PA", "env.1", "env.2")]
  
  # Delete rows where PA == 0
  df_melted <- df_melted[df_melted$PA != 0, ]
  
  # Remove the PA column
  df_melted <- df_melted[, -2]
  
  # Create an empty list to store the results for each species
  hv_list <- list()
  
  # Get unique species values from the "species" column
  species_values <- unique(df_melted$species)
  
  # Iterate over each species
  for (sp in species_values) {
    # Subset the data for the current species
    sim_data <- df_melted[df_melted$species == sp, c("env.1", "env.2")]
    # Calculate the hypervolume
    hv <- hypervolume_gaussian(sim_data, name = sp, samples.per.point = 10)
    
    # Store the hypervolume object in the list
    hv_list[[sp]] <- hv
  }
  
  # Create an empty data frame to store the results
  hv_results <- data.frame(species = character(), hypervolume = numeric(), stringsAsFactors = FALSE)
  
  # Populate the hv_results data frame using the values stored in the list
  for (sp in species_values) {
    sp_hv <- get_volume(hv_list[[sp]])
    hv_results <- rbind(hv_results, data.frame(species = sp, hypervolume = sp_hv))
  }
  

  
  nb <- as.data.frame(niche.breadth)
  nb$species <- row.names(nb)
  #cbind the niche breadth and the hypervolume for each species
  hyp_nb_df <- cbind(as.data.frame(nb$niche.breadth), hv_results$hypervolume)
  colnames(hyp_nb_df) <- c("niche_b", "hypervolume")
  
  return(hyp_nb_df)
}




#### New functions


eliminate_cols <- function(df, min_sum, max_sum) {
  # Get the sums of each column in the data frame
  col_sums <- colSums(df)
  
  # Find the indices of columns to keep
  keep_cols <- which(col_sums >= min_sum & col_sums <= max_sum)
  
  # Return the data frame with only the selected columns
  return(df[, keep_cols])
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

run.nicheROVER <- function(distribution_data, env_variables){
  library(nicheROVER) 
  
  distribution_data <- as.matrix(distribution_data)
  n.env.variables <- NCOL(env_variables)
  n.species <- ncol(distribution_data)
  n.sites <- nrow(distribution_data)
  
  vectorized <- data.frame(as.factor(rep(colnames(distribution_data),each=n.sites)),
                           c(distribution_data),env_variables[rep(1:n.sites, times = n.species), ])
  
  eliminate.zeros <- which(vectorized[,2]==0)
  vectorized <- vectorized[-eliminate.zeros,]
  vectorized <- vectorized[,c(-2)]
  colnames(vectorized) <- c("species",paste0("env.", 1:n.env.variables))
  
  nsamples <- 1000
  vectorized.par <- tapply(1:nrow(vectorized), vectorized$species,
                           function(ii) niw.post(nsamples = nsamples, X = vectorized[ii,2:(n.env.variables+1)]))
  vectorized.size <- sapply(vectorized.par, function(spec) {
    apply(spec$Sigma, 3, niche.size, alpha = .95)
  })
  
  desired_order <- colnames(distribution_data)
  ordered_vectorized.size <- vectorized.size[, match(desired_order, colnames(vectorized.size))]
  nr_nb <- apply(ordered_vectorized.size,2,mean)
  
  #nb
  nb <- as.data.frame(niche.breadth)
  nb$sci.name <- row.names(nb)
  
  
  nr_nb = as.data.frame(nr_nb)
  nr_nb$sci.name <- rownames(nr_nb)
  nr_nb <- merge(nb, nr_nb, by= "sci.name", sort = FALSE)
  colnames(nr_nb) <- c("sci.name", "niche_b", "nicheOVER_ppn")
  
  return(nr_nb)
}








#########################################################################################
########################## New methods to estimate niche breadth ########################
############## methods and functions developed by Pedro Peres-Neto, May 2023 ############
#########################################################################################

fit_gam_mod <- function(data.df, predictors) {
  # Construct the formula for the GAM model
  formula_str <- paste("response ~", paste(paste("s(", predictors, ")", sep = ""), collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  result <- NULL
  warningMessage <- NULL
  
  result <- tryCatch({
    withCallingHandlers({
      gam_model <- mgcv::gam(formula_obj, data = data.df, family = binomial(link = "logit"))
    }, warning = function(w) {
      warningMessage <<- w$message
      invokeRestart("muffleWarning")
    })
    gam_model
  }, error = function(e) {
    e
  })
  
  result <- list(gam_model=gam_model,warning="noWarning")
  if (!is.null(warningMessage)) {
    if (warningMessage %in% c("Iteration limit reached without full convergence - check carefully")) {
      result <- list(gam_model=gam_model,warning="Convergence Issues - check carefully")
    } 
  }
  return(result)
}


fit_gam <- function(data.df, predictors) {
  # Construct the formula for the GAM model
  formula_str <- paste("response ~", paste(paste("s(", predictors, ")", sep = ""), collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  result <- NULL
  warningMessage <- NULL
  
  tryCatch({
    gam_model <- withCallingHandlers({
      mgcv::gam(formula_obj, data = data.df, family = binomial(link = "logit"))
    }, warning = function(w) {
      warningMessage <<- w$message
      invokeRestart("muffleWarning")
    })
    
    # Check for convergence issues
    if (inherits(gam_model, "gam")) {
      if (any(gam_model$converged != 0)) {
        warningMessage <- "Convergence Issues - check carefully"
      }
    }
    
    result <- list(gam_model = gam_model, warning = warningMessage)
  }, error = function(e) {
    result <- list(gam_model = NULL, warning = paste("Error: ", conditionMessage(e)))
  })
  
  return(result)
}


estimate_nicheBreadth_Gam <- function(distribution_data, env_variables) {
  # name env_variables so that the gam function can be generalized
  Env_variables <- data.frame(env_variables)
  colnames(Env_variables) <- paste0("Env", 1:NCOL(env_variables)) # NCOL in case there is only one environmental variable
  
  # this avoids fitting models for species that are present everywhere
  # or in only a few sites
  #distribution_data.smaller <- eliminate_cols(distribution_data, min_sum=10, max_sum=nrow(distribution_data)-10)
  
  n.species <- ncol(distribution_data)
  n.sites <- nrow(distribution_data)
  niche_breadth <- numeric(n.species)
  gam_Warning <- numeric(n.species)
  predicted_values <- matrix(0,n.sites,n.species)
  
  # Loop through each species
  for (i in 1:n.species) {
    data.df <- data.frame(response = distribution_data[, i], Env_variables)
    gam_model <- fit_gam(data.df, predictors = colnames(Env_variables))
    
    gam_Warning[i] <- if (!is.null(gam_model$warning)) gam_model$warning else "noWarning"
    
    if (!is.null(gam_model$gam_model)) {
      gam_model <- gam_model$gam_model
      print(c(i, summary(gam_model)$r.sq))
      predicted_values[, i] <- predict(gam_model, type = "response")
    }
  }
  
  niche.breadth <- apply(as.matrix(dist(t(predicted_values))),2,mean)
  names(niche_breadth) <- paste0("sp", 1:n.species)
  
  results <- data.frame(Niche_Breadth = niche.breadth,
                        Gam_Warning = gam_Warning)
  
  return(results)
}


estimate_nicheBreadth_Gam_mod <- function(distribution_data, env_variables) {
  # name env_variables so that the gam function can be generalized
  Env_variables <- data.frame(env_variables)
  colnames(Env_variables) <- paste0("Env", 1:NCOL(env_variables)) # NCOL in case there is only one environmental variable
  
  # this avoids fitting models for species that are present everywhere
  # or in only a few sites
  #distribution_data.smaller <- eliminate_cols(distribution_data, min_sum=10, max_sum=nrow(distribution_data)-10)
  
  n.species <- ncol(distribution_data)
  n.sites <- nrow(distribution_data)
  niche_breadth <- numeric(n.species)
  gam_Warning <- numeric(n.species)
  predicted_values <- matrix(0,n.sites,n.species)
  # Loop through each species
  for (i in 1:n.species) {
    #print(i)
    # Fit the GAM with a Poisson distribution and smooth term for the environmental variable
    # gam_model <- gam(abundance_data[, i] ~ s(env_variable), family = poisson(link = "log"))
    # Fit the GAM with a binomial distribution and logistic link for the environmental variable
    data.df <- data.frame(response=distribution_data[, i],Env_variables)
    gam_model <- fit_gam(data.df, predictors=colnames(Env_variables))
    gam_Warning[i] <- gam_model$warning
    gam_model <- gam_model$gam_model
    print(c(i,summary(gam_model)$r.sq))
    
    predicted_values[,i] <- predict(gam_model, type = "response")
    
  }
  niche.breadth <- apply(as.matrix(dist(t(predicted_values))),2,mean)
  names(niche_breadth) <- paste0("sp", 1:n.species)
  
  results <- data.frame(Niche_Breadth = niche.breadth,
                        Gam_Warning = gam_Warning)
  
  return(results)
}

estimate_nicheBreadth_Latents <- function(distribution_data, env_variables,nlv=5){
  n.species <- ncol(distribution_data)
  distribution_data.smaller <- eliminate_cols(distribution_data, 20, 480)
  eco_PA <- ecoCopula::stackedsdm(distribution_data.smaller, ~ 1,data=distribution_data,family="binomial")
  eco_lvs <- ecoCopula::cord(eco_PA, nlv=nlv)
  eco_lvs <- matrix(eco_lvs$scores,nrow(distribution_data.smaller),nlv)
  #predicted_values <- matrix(0,n.sites,n.species)
  predicted_values <- matrix(0, nrow = nrow(distribution_data.smaller), ncol = n.species)
  # Initialize results data frame
  # Loop through each species
  for (i in 1:n.species) {
    # print(i)
    # Fit the GAM with a binomial distribution and logistic link based on the latents; non linear for now
    logistic_model <- mgcv::gam(distribution_data[, i] ~ eco_lvs, family = binomial(link = "logit"))
    
    print(c(i,summary(logistic_model)$r.sq))
    predicted_values[,i] <- predict(logistic_model, type = "response")
    
  }
  # Return results
  niche.breadth <- apply(as.matrix(dist(t(predicted_values))),2,mean)
  names(niche.breadth) <- paste0("sp", 1:n.species)
  
  
  return(niche.breadth)
}

estimate_nicheBreadth_avg.Dist <- function(distribution_data, env_variables){
  library(vegan)
  spe.Hellinger <- decostand(distribution_data, method = "hellinger")
  niche.breadth <- apply(as.matrix(dist(t(spe.Hellinger))),2,mean)
  names(niche.breadth) <- paste0("sp", 1:n.species)
  
  
  return(niche.breadth)
}



