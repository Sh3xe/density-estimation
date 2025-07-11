library(fdaPDE)
library(viridis)
rm(list = ls())
directory <- "data/infections_southampton"

# Import helper functions to sample and visualize data
source("scripts/helper_functions_plot.R")
source("scripts/helper_functions_data.R")

nodes <- read.csv(file.path(directory,"mesh_vertices.csv"))[,c(2,3)]
triangles <- read.csv(file.path(directory, "mesh_elements.csv"))[,c(2,3,4)]
nodesmarker <- read.csv(file.path(directory, "mesh_boundary.csv"))[,2]

# Create a regular mesh of the unit square
mesh <- create.mesh.2D(nodes = nodes, triangles = triangles, nodesattributes = nodesmarker)

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Load the data
f_init <- read.csv(file.path(directory, "f_init.csv"))
sample <- read.csv(file.path(directory, "sample.csv"))
log_dens <- read.csv(file.path(directory, "log_dens.csv"))

# Plot
FEMfunction <- FEM(coeff = f_init[,2], FEMbasis = FEMbasis)
evaluation <- eval.FEM(FEM = FEMfunction, locations = mesh$nodes)
estimated_density <- exp(evaluation)
plot.density.2D.map(coeff = estimated_density, mesh = mesh, colorscale = viridis)