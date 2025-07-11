library(fdaPDE)
rm(list = ls())

# This script is located under [repo/scripts] and the data is located at [repo/data] 
setwd(file.path(getwd(), ".."))

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")

# Load data
load(file = "data/infections_southampton.RData")

# Define the domain
nodes <- cbind(as.data.frame(xyt$window$bdry[[1]]$x), as.data.frame(xyt$window$bdry[[1]]$y))
segments <- cbind(1:dim(nodes)[1], c(2:dim(nodes)[1], 1))

# Define the mesh
mesh <- create.mesh.2D(nodes = nodes, segments = segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 4.25, minimum_angle = 25)

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Import the data
data <- cbind(xyt$x, xyt$y)

# Plot
plot.sample.2D.map(data = data, mesh = mesh)

# Save the mesh and samples
output_dir <- "data/infections_southampton"
write.csv(data, file.path(output_dir, "sample.csv"))
write.csv(mesh$nodes, file.path(output_dir,"mesh_vertices.csv"))
write.csv(mesh$triangles, file.path(output_dir, "mesh_elements.csv"))
write.csv(mesh$nodesmarker, file.path(output_dir, "mesh_boundary.csv"))

# Calculate f_init
lambda_proposal <- 10^seq(from = -1, to = -5, by = -0.5)

de <- DE.FEM(
	data=data,
	FEMbasis = FEMbasis,
	lambda = lambda_proposal,
	nfolds = 10,
	tol1 = 1e-4,
	nsimulations = 1000,
	preprocess_method = "RightCV",
	step_method = "Wolfe_Method",
	direction_method = "L-BFGS10"
)

best_lambda_id <- match(de$lambda, lambda_proposal)

write.csv(de$f_init[,best_lambda_id], file.path(output_dir, "f_init.csv"))