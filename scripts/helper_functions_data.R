#'############################################################################'#
#'#################### HELPER FUNCTIONS FOR DATA SAMPLING ####################'#
#'############################################################################'#

library(sf)
library(mvtnorm)
library(mgcv)
library(Directional)
library(stlnpp)
library(sfnetworks)
library(tidygraph)
library(spatstat)

### 1.1 SIMULATION STUDY ON [0,1]x[0,1] ----------------------------------------
# Get domain boundary
domain.1 <- st_sfc(st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE))))

# Get mixing weights
get.priors.1 <- c(1,1,1,1)/4

# Get mean and covariance matrix of the d-th component of the mixture
get.parameters.1 <- function(d){
  if(d == 1){
    # First Gaussian distribution
    mu = c(1, 1)/3 # Mean
    Sigma = matrix(data = c(0.8, -0.5, -0.5, 1)/100, nrow = 2, ncol = 2) # Covariance Matrix
  } else if(d == 2){
    # Second Gaussian distribution
    mu = c(2, 1)/3 # Mean
    Sigma = matrix(data = c(1.5, 0, 0, 1.5)/100, nrow = 2, ncol = 2) # Covariance Matrix
  } else if(d == 3){
    # Third Gaussian distribution
    mu = c(1, 2)/3 # Mean
    Sigma = matrix(data = c(0.6, 0, 0, 0.6)/100, nrow = 2, ncol = 2) # Covariance Matrix
  } else if(d == 4){
    # Fourth Gaussian distribution
    mu = c(2, 2)/3 # Mean
    Sigma = matrix(data = c(1, 0.8, 0.8, 1)/100, nrow = 2, ncol = 2) # Covariance Matrix
  }
  
  return (list(mu, Sigma))
}

# Generate purely spatial point patterns sampling N data
generate.spatial.data.1 <- function(N){
  
  # Create working directory
  dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
  
  if(!file.exists(paste0("data/","sim1_",N,"data.txt"))){
    
    # Data
    data = c()
    
    # Number of data sampled
    n = 1
    
    # Generate N data
    while (n <= N){
      p = get.priors.1
      d = sample(c(1,2,3,4), 1, prob = p)
      parameters = get.parameters.1(d)
      
      l = mvtnorm::rmvnorm(n = 1, mean = parameters[[1]], sigma = parameters[[2]], method = "svd")
      
      is_within = st_within(st_sfc(st_point(l)), domain.1)
      if(lengths(is_within) > 0){
        if(is_within[[1]]){
          data = rbind(data, l)
          n = n+1
        }
      }
      
    }
    
    # Export data
    write.table(data, paste0("data/","sim1_",N,"data.txt"), row.names = F, col.names = F)
    # Import data
    return(data)
    
  } else {
    
    # Import data
    data = read.table(file = paste0("data/","sim1_",N,"data.txt"))
    return(data)
    
  }
}

# Density function
dens.func.1 <- function(data){
  p = get.priors.1
  dens = p[1]*mvtnorm::dmvnorm(data, mean = get.parameters.1(1)[[1]], sigma = get.parameters.1(1)[[2]]) +
    p[2]*mvtnorm::dmvnorm(data, mean = get.parameters.1(2)[[1]], sigma = get.parameters.1(2)[[2]]) +
    p[3]*mvtnorm::dmvnorm(data, mean = get.parameters.1(3)[[1]], sigma = get.parameters.1(3)[[2]]) +
    p[4]*mvtnorm::dmvnorm(data, mean = get.parameters.1(4)[[1]], sigma = get.parameters.1(4)[[2]])
  
  return (dens)
}


### 1.2 SIMULATION STUDY ON A HORSESHOE-SHAPED DOMAIN --------------------------
# Get domain boundary
boundary.domain.2 <- function(r0 = 0.0875, r = 0.50625, l = 3, n.theta = 20){
  rr = r + (r - r0)
  theta = seq(pi, pi/2, length = n.theta/2)
  x = rr * cos(theta)
  y = rr * sin(theta)
  x = c(x, seq(0.1, 2.9, length = 10))
  y = c(y, rep(y[length(y)],10))
  theta = seq(pi/2, -pi/2, length = n.theta/2)
  x = c(x, (r-r0) * cos(theta) + l)
  y = c(y, (r-r0) * sin(theta) + r)
  x = c(x, seq(2.9, 0.1, length = 10))
  y = c(y, rep(y[length(y)],10))
  theta = seq(pi/2, pi, length = n.theta/5)
  x = c(x, r0 * cos(theta))
  y = c(y, r0 * sin(theta))
  n = length(x)
  x = c(x, x[n:1])
  y = c(y, -y[n:1])
  return (unique(cbind(x, y)))
}

domain.2 <- st_sfc(st_polygon(list(rbind(boundary.domain.2(), boundary.domain.2()[1,]))))

# Sample N points uniformly
sample.uniform.2 <- function(N){
  data = c()
  n = 1
  while (n <= N) {
    x = runif(n = 1, min = min(boundary.domain.2()[,1]), max = max(boundary.domain.2()[,1]))
    y = runif(n = 1, min = min(boundary.domain.2()[,2]), max = max(boundary.domain.2()[,2]))
    if (mgcv::inSide(bnd = list(x = boundary.domain.2()[,1], y = boundary.domain.2()[,2]),
                     x = x, y = y)) {
      data = rbind(data, c(x,y))
      n = n+1
    }
  }
  
  return(data)
}

# Generate purely spatial point patterns sampling N data
generate.spatial.data.2 <- function(N, mesh){
  
  # Create working directory
  dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
  
  if(!file.exists(paste0("data/","sim2_",N,"data.txt"))){
    
    # Data
    data = c()
    
    # Number of data sampled
    n = 1
    
    data = c()
    
    # Generate N data
    while (n <= N){
      
      # Sample uniformly
      data_unif = sample.uniform.2(N)
      
      # Evaluate the density function on the generated points
      f <- dens.func.2(data_unif, mesh)
      f[is.na(f)] = max(f, na.rm = TRUE)
      
      # For each generated point, decide to keep it or discard it on the base of the 
      # criterium: runif(n = 1, min = 0, max = 1) < density.value
      for(i in 1:N){
        if(runif(n = 1, min = 0, max = 1) < f[i]){
          data = rbind(data, data_unif[i,])
          n = n+1
        }
      }
    }
    
    data = data[1:N,]
    
    # Export data
    write.table(data, paste0("data/","sim2_",N,"data.txt"), row.names = F, col.names = F)
    
  } else {
    
    # Import data
    data = read.table(file = paste0("data/","sim2_",N,"data.txt"))
    
  }
  
  return(data)
}

# Density function
dens.func.2 <- function(data, mesh, time){
  
  bnd = list(x = boundary.domain.2()[,1], y = boundary.domain.2()[,2])
  
  dens = vector(mode = "numeric", length = nrow(data))
  
  for(i in 1:length(dens)){
    
    x = data[i,1]
    y = data[i,2]
    
    if (mgcv::inSide(bnd = bnd, x = x, y = y)) {
      
      dens[i] = mgcv::fs.test(x = x, y = y, b = 1) + 5
      
    } else {
      
      dens[i] = NA
      
    }
  }
  
  return (dens / integral.dens.func.2(mesh = mesh))
}

# Compute the integral of the test function (to get a proper density)
integral.dens.func.2 <- function(mesh){
  return(sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(create.FEM.basis(mesh = mesh)) %*% (mgcv::fs.test(x = mesh$nodes[,1], y = mesh$nodes[,2], b = 1) + 5), na.rm = TRUE))
}


### 2.1 SIMULATION STUDY ON THE SURFACE OF THE UNIT SPHERE ---------------------
# Parameters of the mixture
mu1 = c(-0.5, -0.5, 0.8) 
mu1 = mu1 / sqrt( sum(mu1^2) )
gamma11 = c(-0.7789378,  0.6157424,  0.1188163)
gamma12 = c(-0.5695773, -0.6154000, -0.5448528)
k1 = 18
beta1 = 0

mu2 = c(-0.3, -0.3, 0.2)
mu2 = mu2 / sqrt( sum(mu2^2) )
gamma21 = c(-0.8651146,  0.3803316, -0.3269933)
gamma22 = c( 0.1482597, -0.4288975, -0.8911038)
k2 = 15
beta2 = 7

mu3 = c(0.5, -0.5, 0.8)
mu3 = mu3 / sqrt( sum(mu3^2) )
gamma31 = c(-0.66647307, -0.74323532, -0.05843723)
gamma32 = c( 0.5753645,  -0.4629244,  -0.6742824)
k3 = 20
beta3 = 10

mu4 = c(0.2, -1, 0)
mu4 = mu4 / sqrt( sum(mu4^2) )
gamma41 = c( 0.6364658, -0.0303920, -0.7707059)
gamma42 = c(-0.7545879, -0.2314437, -0.6140285)
k4 = 20
beta4 = 7

mu5 = c(0.6, -0.5, 0.3)
mu5 = mu5 / sqrt( sum(mu5^2) )
gamma51 = c( 0.6364658, -0.0303920, -0.7707059)
gamma52 = c( 0.7545879, -0.2314437, -0.6140285)
k5 = 20
beta5 = 4

# Sample N points uniformly
sample.uniform.3 <- function(N){
  theta = 2*pi*runif(n = N, min = 0, max = 1)
  phi = acos(1-2*runif(n = N, min = 0, max = 1))
  
  x = sin(phi) * cos(theta)
  y = sin(phi) * sin(theta)
  z = cos(phi)
  
  return(cbind(x,y,z))
}

# Density function
dens.func.3 <- function(data){
  
  G1 <- cbind(mu1, gamma11, gamma12)
  G2 <- cbind(mu2, gamma21, gamma22)
  G3 <- cbind(mu3, gamma31, gamma32)
  G4 <- cbind(mu4, gamma41, gamma42)
  G5 <- cbind(mu5, gamma51, gamma52)
  
  return (dkent(data, G1, param = c(k1, beta1)) / 5 +
            dkent(data, G2, param = c(k2, beta2)) / 5 + 
            dkent(data, G3, param = c(k3, beta3)) / 5 + 
            dkent(data, G4, param = c(k4, beta4)) / 5 +
            dkent(data, G5, param = c(k5, beta5)) / 5 )
  
}

# Generate purely spatial point patterns sampling N data
generate.spatial.data.3 <- function(N, density.function = dens.func){
  
  # Create working directory
  dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
  
  if(!file.exists(paste0("data/","sim3_",N,"data.txt"))){
    
    # Data
    data = c()
    
    # Number of data sampled
    n = 1
    
    data = c()
    
    # Generate N data
    while (n <= N){
      
      # Sample uniformly
      data_unif = sample.uniform.3(N)
      
      # Evaluate the density function on the generated points
      f <- dens.func.3(data_unif)
      f[is.na(f)] = max(f, na.rm = TRUE)
      
      # For each generated point, decide to keep it or discard it on the base of the 
      # criterium: runif(n = 1, min = 0, max = 1) < density.value
      for(i in 1:N){
        if(runif(n = 1, min = 0, max = 1) < f[i]){
          data = rbind(data, data_unif[i,])
          n = n+1
        }
      }
    }
    
    data = data[1:N,]
    
    # Export data
    write.table(data, paste0("data/","sim3_",N,"data.txt"), row.names = F, col.names = F)
    
  } else {
    
    # Import data
    data = read.table(file = paste0("data/","sim3_",N,"data.txt"))
    
  }
  
  return (data)
}


### 2.2 SIMULATION STUDY ON A CURVED SURFACE -----------------------------------
# Sample N points uniformly
sample.uniform.4 <- function(N){
  load("data/sim4.fullPoints.proj.RData")
  numFullPoints = nrow(fullPoints.proj)
  unifPointsIndex = sample(1:numFullPoints, N)
  unifPoints = fullPoints.proj[unifPointsIndex,]
  
  return (unifPoints)
}

# Density function
dens.func.4 <- function(data, mesh){
  fun = ((sin(data[,1]) + cos(data[,2]) + data[,3] + 5.16)^3)
  
  return (fun / integral.dens.func.4(mesh = mesh))
}

# Compute the integral of the test function (to get a proper density)
integral.dens.func.4 <- function(mesh){
  return(sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(create.FEM.basis(mesh = mesh)) %*% (((sin(mesh$nodes[,1]) + cos(mesh$nodes[,2]) + mesh$nodes[,3] + 5.16)^3)), na.rm = TRUE))
}

# Generate purely spatial point patterns sampling N data
generate.spatial.data.4 <- function(N, mesh){
  
  # Create working directory
  dir.create(file.path(getwd(), "data"), showWarnings = FALSE)
  
  if(!file.exists(paste0("data/","sim4_",N,"data.txt"))){
    
    # Data
    data = c()
    
    # Number of data sampled
    n = 1
    
    data = c()
    
    # Generate N data
    while (n <= N){
      
      # Sample uniformly
      data_unif = sample.uniform.4(N)
      
      # Evaluate the density function on the generated points
      f <- dens.func.4(data_unif, mesh)
      f[is.na(f)] = max(f, na.rm = TRUE)
      
      # For each generated point, decide to keep it or discard it on the base of the 
      # criterium: runif(n = 1, min = 0, max = 1) < density.value
      for(i in 1:N){
        if(runif(n = 1, min = 0, max = 1) < f[i]){
          data = rbind(data, data_unif[i,])
          n = n+1
        }
      }
    }
    
    data = data[1:N,]
    
    # Export data
    write.table(data, paste0("data/","sim4_",N,"data.txt"), row.names = F, col.names = F)
    
  } else {
    
    # Import data
    data = read.table(file = paste0("data/","sim4_",N,"data.txt"))
    
  }
  
  return(data)
}


### 3.1 SIMULATION STUDY ON A LINEAR NETWORK -----------------------------------
mesh.5 <- function(){
  
  data("Eastbourne", package = "stlnpp")
  
  # Build sfnetwork
  sfnetwork <- as_sfnetwork(Eastbourne$domain, directed = FALSE, edges_as_lines = TRUE)
  
  # Simplify 
  sfnetwork <- sfnetwork %>%
    activate("edges") %>%
    filter(!edge_is_multiple()) %>%
    filter(!edge_is_loop())
  
  # Clean
  sfnetwork <- sfnetwork %>% 
    convert(to_spatial_subdivision, .clean = TRUE)
  
  # Selecting full connected graph
  sfnetwork <- sfnetwork %>% 
    convert(to_components, .clean = TRUE, .select = 1L)
  
  mesh = as.mesh.1.5D(sfnetwork)
  mesh = normalize.mesh.unit(mesh)$mesh
  
  FEMbasis = create.FEM.basis(mesh)
  
  return(mesh)
  
}

# Convert mesh.1.5D into spatstat.linnet
as.linnet.mesh.1.5D <- function(X,...){
  vertices = ppp(x = X$nodes[,1], 
                 y = X$nodes[,2], 
                 window = owin(xrange = c(min(X$nodes[,1]),max(X$nodes[,1])),
                               yrange = c(min(X$nodes[,2]),max(X$nodes[,2]))))
  spat.stat.linnet = linnet(vertices = vertices, edges = X$edges, sparse = T)
  
  return(spat.stat.linnet)
}

# Convert sfnetwork into mesh.1.5D
setGeneric("as.mesh.1.5D", function(x) standardGeneric("as.mesh.1.5D"))

setMethod("as.mesh.1.5D", signature = "sfnetwork", function(x){
  nodes = st_coordinates(x, "nodes")
  edges = st_as_sf(x, "edges")
  edges = cbind(edges$from, edges$to) 
  mesh = create.mesh.1.5D(nodes, edges)
  return(mesh)
})

# Mesh normalization
normalize.mesh.unit <- function(mesh){
  
  x.min = min(mesh$nodes[,1])
  y.min = min(mesh$nodes[,2])
  
  x.max = max(mesh$nodes[,1])
  y.max = max(mesh$nodes[,2])
  
  x.norm = (mesh$nodes[,1] - x.min)/(x.max - x.min)
  y.norm = (mesh$nodes[,2] - y.min)/(y.max - y.min)
  
  mesh.norm = create.mesh.1.5D(nodes = cbind(x.norm,y.norm), edges = mesh$edges)  
  
  ret = list(mesh = mesh.norm,
             x.min = x.min, y.min = y.min, x.max = x.max, y.max = y.max)
  
  return(ret)
}

# Data normalization
normalize.data.unit <- function(data, mesh){
  
  x.min = min(mesh$nodes[,1])
  y.min = min(mesh$nodes[,2])
  
  x.max = max(mesh$nodes[,1])
  y.max = max(mesh$nodes[,2])
  
  data.x.norm = (data[,1] - x.min)/(x.max - x.min)
  data.y.norm = (data[,2] - y.min)/(y.max - y.min)
  
  data.norm = data.frame(x = data.x.norm, y = data.y.norm)
  
  return(data.norm)
}

# Intensity function
intens.func.5 <- function(x, y, L = mesh_linnet, sigma = 0.2){
  
  current_sources <- as.matrix(50, 1, 1)
  
  nodes.lpp <- ppp(x = L$vertices$x, y = L$vertices$y, window = L$window)
  
  PP <- ppp(x = x, y = y, window = nodes.lpp$window)
  ND <- crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  coeff = 100 * 1/sqrt(2*pi*sigma^2) * exp(-ND[current_sources[,1],]^2/(2*sigma^2))
  
  return (coeff)
}

# Compute the integral of the intensity function (to get a proper density)
integral.dens.func.5 <- function(f){
  
  integral <- sum(fdaPDE:::CPP_get.FEM.Mass.Matrix(f) %*%
                    intens.func.5(x = f$mesh$nodes[,1], y = f$mesh$nodes[,2]))
  
  return (integral)
  
}
